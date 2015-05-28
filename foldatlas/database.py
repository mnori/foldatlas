# from app import db

from sqlalchemy import create_engine, and_
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

import settings, os
import sys

engine = create_engine(settings.database_uri, convert_unicode=True)
db_session = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))
Base = declarative_base()
Base.query = db_session.query_property()

import models

from models import Strain, Gene, Transcript, Feature, AlignmentEntry, ReactivityMeasurement, GeneLocation

def hydrate_db():
    try:
        print("Rebuilding schema...")
        Base.metadata.drop_all(bind=engine)
        Base.metadata.create_all(bind=engine)

        # Add the annotations
        SequenceHydrator().hydrate() 

        # Add DMS reactivities
        DmsReactivityHydrator().hydrate() 

        # Do alignments so we can see polymorphism
        # This takes a pretty long time, will probably have to run on HPC when doing the real thing.
        # TranscriptAligner().align() 

        print("Hydration Complete.")

    except Exception as e: # catch the exception so we can display a nicely formatted error message
        print(str(e).replace("\\n", "\n").replace("\\t", "\t"))
        raise e

# Parses genome sequence .fa and annotation .gff3 files into the database.
class SequenceHydrator():

    # how many genes to process before committing rows to the database.
    gene_chunk_size = 2500 

    genes_to_write = []
    gene_locations_to_write = []
    transcripts_to_write = []
    features_to_write = []

    genes_seen = {}
    transcripts_seen = {}

    # for duplicate transcript ID detection
    transcript_ids_seen_this_strain = set()

    # limit on genes to process - for testing purposes
    gene_limit = None

    # limit on chromosome sequence to add, in bp - for testing
    bp_limit = None

    # max strains
    strain_limit = 1

    # Use the genome sequence and annotation files to populate the database.
    def hydrate(self):
        n_strains = 0;
        for strain in settings.strains:
            self.hydrate_strain(strain)
            n_strains += 1
            if n_strains >= self.strain_limit:
                break
        db_session.commit()

    def hydrate_strain(self, strain_config):
        self.transcript_ids_seen_this_strain = set()

        print("Hydrating strain ["+strain_config["name"]+"]")

        # add the strain
        strain = Strain(id=strain_config["name"], description=strain_config["description"])
        db_session.add(strain)
        db_session.commit()

        # add the chrosomomes
        self.hydrate_chrosomomes(strain_config)

        # add genes, transcripts, and feature annotations
        self.hydrate_genes(strain_config)

    # Adding chromosomes to the DB is a little bit tricky, since the sequences are huge.
    # Therefore a LOAD DATA INFILE strategy is used to import the data.
    def hydrate_chrosomomes(self, strain_config):

        print("Adding chromosomes...")
        filepath = settings.genomes_sauce_folder+"/"+strain_config["sequence_filename"]

        for record in SeqIO.parse(filepath, "fasta"): # loop through chromosomes

            chr_id = record.id
            if (chr_id in settings.ignored_chromosomes):
                continue

            seq_str = str(record.seq)

            temp_filepath = settings.temp_folder+"/tmp.fa"

            # Save a row of chromosome data to a text file
            temp_file = open(temp_filepath, "w")
            temp_file.write(strain_config["name"]+"\t"+chr_id+"\t"+seq_str)
            temp_file.close()

            # Import file into the DB
            sql = """
                LOAD DATA LOCAL INFILE '/tmp/tmp.fa'
                REPLACE INTO TABLE chromosome
            """
            db_session.execute(sql)
            db_session.commit()

            # Delete the file
            os.remove(temp_filepath)

            print("Added ["+chr_id+"]")

        print("Finished adding chromosomes to ["+strain_config["name"]+"]")

    def hydrate_genes(self, strain_config):
        genes_added = 0

        feature_rows = []

        # open the annotation file and go through it line by line
        with open(settings.genomes_sauce_folder+"/"+strain_config["annotation_filename"]) as gff_file:
            for gff_line in gff_file:
                if gff_line[0] == "#": # ignore comments
                    continue

                bits = gff_line.split("\t")
                feature_type = bits[2]
                if feature_type == "gene":
                    
                    if len(feature_rows) > 0: # this is needed to stop it going wrong at the beginning
                        self.hydrate_gene(feature_rows, strain_config["name"])
                        genes_added += 1
                        feature_rows = []
                        if genes_added % 100 == 0:
                            if self.gene_limit != None and genes_added >= self.gene_limit:
                                break
                            print (str(genes_added)+" genes processed")

                        if genes_added % self.gene_chunk_size == 0: # commit at regular intervals
                            self.commit_all()

                feature_rows.append(bits)

        # gotta add that last entry
        if (len(feature_rows) > 0):
            self.hydrate_gene(feature_rows, strain_config["name"])
            genes_added += 1

        self.commit_all()

        # do the sequences
        print (str(genes_added)+" genes added total")


    def commit_all(self):
        self.commit_entities_list(self.genes_to_write, "Genes")
        self.commit_entities_list(self.gene_locations_to_write, "GeneLocations")
        self.commit_entities_list(self.transcripts_to_write, "Transcripts")
        self.commit_entities_list(self.features_to_write, "Features")

        self.genes_to_write = []
        self.gene_locations_to_write = []
        self.transcripts_to_write = []
        self.features_to_write = []
        
    def commit_entities_list(self, entities, label):
        print("Committing "+label+"...")
        for entity in entities:
            db_session.add(entity)
        db_session.commit()
        print("...done.")

    def hydrate_gene(self, feature_rows, strain_id):
        features = {}
        sequence = None
        transcript = None

        gene_id = None
        min_start = None
        max_end = None

        for feature_row in feature_rows: # Loop through annotation rows in the gff file, all related to the current gene

            # keep track of start and end
            start = feature_row[3]
            end = feature_row[4]
            direction = "forward" if feature_row[6] == "+" else "reverse"
            chromosome_id = feature_row[0]

            if min_start == None or start < min_start:
                min_start = start
            if max_end == None or end > max_end:
                max_end = end

            feature_type = feature_row[2]
            attribs = feature_row[8].strip()

            if feature_type == "gene": # Handle gene entries
                gene_id = attribs.split(";")[0].split(":")[1] # grab the gene ID - we'll want this for later

                # add the Gene entry - if it hasn't been already
                if gene_id not in self.genes_seen: 
                    gene = Gene(gene_id)
                    self.genes_to_write.append(gene)
                    self.genes_seen[gene_id] = gene

            
            else: # Handle transcript entries - only add new ones
                transcript_id = self.find_attribs_value("ID=Transcript", attribs)
                if transcript_id != None: # it's a transcript entry

                    # add the Transcript entry - if it hasn't been already
                    transcript_id = self.ensure_unique_transcript_id(transcript_id)

                    if transcript_id not in self.transcripts_seen: 
                        transcript = Transcript(
                            id=transcript_id, gene_id=gene_id
                        )
                        self.transcripts_to_write.append(transcript)
                        self.transcripts_seen[transcript.id] = transcript

                else: # Handle transcript feature entries

                    # for some reason, features for a given strain/transcript 
                    # combination are not always added

                    transcript_id = self.find_attribs_value("Parent=Transcript", attribs)

                    if transcript_id != None: # it's a transcript feature entry
                        # put a filter here? some elements are not worth storing?
                        self.features_to_write.append(Feature(
                            transcript_id=transcript_id,
                            type_id=feature_row[2],
                            strain_id=strain_id,
                            chromosome_id=chromosome_id,
                            start=start,
                            end=end,
                            direction=direction
                        ))

                    else:
                        pass # this happens for pseudogenes and TEs - which we aint interested in

        # Add to the gene location cache table
        self.gene_locations_to_write.append(GeneLocation(
            gene_id=gene_id, 
            strain_id=strain_id, 
            chromosome_id=chromosome_id, 
            start=min_start, 
            end=max_end, 
            direction=direction
        ))


    def ensure_unique_transcript_id(self, transcript_id):
        version = 1
        candidate_transcript_id = transcript_id
        while True:
            if candidate_transcript_id in self.transcript_ids_seen_this_strain:
                version += 1
                candidate_transcript_id = transcript_id+"_v"+str(version)
            else:
                self.transcript_ids_seen_this_strain.add(transcript_id)
                if candidate_transcript_id != transcript_id:
                    print("Transcript ID ["+transcript_id+"] was a duplicate, renamed to ["+candidate_transcript_id+"]")
                return candidate_transcript_id
        
    # Parse out the value of a key in the attribs field
    # e.g. 
    #   find_attribs_value("Parent=Transcript", "ID=five_prime_UTR:AT5G67630.1.1;Parent=Transcript:AT5G67630.1")
    # will return
    #   AT5G67630.1
    #   
    def find_attribs_value(self, key, attribs_str):
        entries = attribs_str.split(";")
        for entry in entries:
            entry_bits = entry.split(":")
            if (entry_bits[0] == key):
                return ":".join(entry_bits[1:]) # we need all of the bits in the array
        return None

# Class for doing alignments, one run per transcript.
class TranscriptAligner():

    def align(self):
        transcript_ids = self.fetch_transcript_ids()
        for transcript_id in transcript_ids:
            self.process_transcript_id(transcript_id)
            
    def process_transcript_id(self, transcript_id):

        print("Aligning ["+transcript_id+"]...", end="")
        sys.stdout.flush()

        seqs_to_align = list(Transcript(transcript_id).get_sequences().values())

        if len(seqs_to_align) <= 1:
            print("Warning - not enough sequences to proceed with alignment")
            return

        temp_filepath = settings.temp_folder+"/tmp.fa"

        # output to a fasta file for clustalw alignment
        output_handle = open(temp_filepath, "w")
        SeqIO.write(seqs_to_align, output_handle, "fasta")
        output_handle.close()

        # run the clustalw alignment
        clustalw_cline = ClustalwCommandline("clustalw2", infile=temp_filepath, quicktree=True)
        results = clustalw_cline()

        # parse the results into the database
        entries = AlignIO.read(settings.temp_folder+"/tmp.aln", "clustal")
        for entry in entries:
            obj = AlignmentEntry(transcript_id, entry.id, str(entry.seq))
            db_session.add(obj)
            
        db_session.commit()

        print("Aligned")

    # Fetch all the transcript IDs from the database. Order them for consistency
    def fetch_transcript_ids(self):
        transcript_ids = []
        sql = "SELECT id FROM transcript ORDER BY id ASC"
        rows = engine.execute(sql)
        for row in rows:
            transcript_ids.append(row["id"])

        return transcript_ids

# Inserts DMS reactivities into the DB.
class DmsReactivityHydrator():

    def hydrate(self):

        print("Adding reactivities...")

        # fetch all of the transcript IDs from the database, store them in a set to check against.
        sql = ("SELECT id FROM transcript ORDER BY id ASC")
        results = engine.execute(sql)
        transcript_ids = set()
        for result in results:
            transcript_ids.add(result["id"])

        # Open the DMS reactivities file. These are normalised already.
        input_file = open(settings.dms_reactivities_filepath, "r")
        line = input_file.readline()
        while(True):
            transcript_id = line.strip()
            line = input_file.readline().strip()

            if line[:2] == "AT": # another transcript id
                continue;

            if line == "": # if empty, it means end of file was reached
                break

            # got this far? assume it's reactivities
            if transcript_id not in transcript_ids:
                continue;

            count_strs = line.split("\t")

            # go through reactivity entries, adding each to the database.
            position = 0;
            for count_str in count_strs:
                position += 1
                if (count_str != "NA"): # skip adding "NA" entries.
                    obj = ReactivityMeasurement(
                        strain_id=settings.reference_strain_id, 
                        transcript_id=transcript_id, 
                        position=position, 
                        reactivity=float(count_str))
                    db_session.add(obj)        
            db_session.commit() # insert all the reactivity measurement rows into the DB
            print("Added ["+transcript_id+"] ("+str(position)+" positions)")

        input_file.close()

# experiment
#     id 
#     description

# experiment_transcript 
#     experiment_id
#     strain_id
#     transcript_id

# experiment_transcript_reactivity
#     - basically the current table
