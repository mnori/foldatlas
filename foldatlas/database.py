# from app import db

from sqlalchemy import create_engine, and_
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from sklearn import decomposition

import settings, os
import sys
import re

engine = create_engine(settings.database_uri, convert_unicode=True)
db_session = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))
Base = declarative_base()
Base.query = db_session.query_property()

import models

from models import Strain, Gene, Transcript, Feature, AlignmentEntry, NucleotideMeasurement, \
    GeneLocation, Experiment, TranscriptCoverage, Structure, StructurePosition

def hydrate_db():
    try:

        # print("Rebuilding schema...")
        # Base.metadata.drop_all(bind=engine)
        # Base.metadata.create_all(bind=engine)

        # # Add the annotations
        # SequenceHydrator().hydrate() 

        # # Add DMS reactivities
        # NucleotideMeasurementHydrator().hydrate(settings.dms_reactivities_experiment)
        # CoverageHydrator().hydrate(settings.dms_reactivities_experiment)

        # # Add ribosome profiling
        # NucleotideMeasurementHydrator().hydrate(settings.ribosome_profile_experiment)
        # CoverageHydrator().hydrate(settings.ribosome_profile_experiment)

        # # Do alignments so we can see polymorphism
        # # This takes a pretty long time, will probably have to run on HPC when doing the real thing.
        # TranscriptAligner().align() 

        # # Import RNA structures
        # StructureHydrator().hydrate(settings.structures_in_silico)
        # StructureHydrator().hydrate(settings.structures_in_vivo)

        PcaHydrator().hydrate(settings.structures_in_silico)

        print("Hydration Complete.")

    except Exception as e: # catch the exception so we can display a nicely formatted error message
        print(str(e).replace("\\n", "\n").replace("\\t", "\t"))
        raise e

# Parses genome sequence .fa and annotation .gff3 files into the database.
class SequenceHydrator():

    # how many genes to process before committing rows to the database.
    gene_chunk_size = 2500 

    genes_to_write = []
    transcripts_to_write = []
    features_to_write = []

    genes_seen = {}
    transcripts_seen = {}

    # for duplicate transcript ID detection
    transcript_ids_seen_this_strain = set()

    # limit on genes to process - for testing purposes
    gene_limit = 10

    # limit on chromosome sequence to add, in bp - for testing
    bp_limit = None

    gene_location_chunk_size = 1000

    # max strains
    strain_limit = None

    # Use the genome sequence and annotation files to populate the database.
    def hydrate(self):
        n_strains = 0
        for strain in settings.strains:
            self.hydrate_strain(strain)
            n_strains += 1
            if self.strain_limit != None and n_strains >= self.strain_limit:
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

        self.cache_gene_locations(strain_config)

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
                feature_type = bits[2] # BUG
                if feature_type == "gene": 
                    
                    if len(feature_rows) > 0: # this is needed to stop it going wrong at the beginning
                        self.hydrate_gene(feature_rows, strain_config["name"])
                        genes_added += 1
                        feature_rows = []
                        if self.gene_limit != None and genes_added >= self.gene_limit:
                            break
                        if genes_added % 100 == 0:
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

    # Cache gene locations in a redundant table by looking at the feature locations.
    def cache_gene_locations(self, strain_config):
        print("Caching gene locations...")
        start = 0
        while(True):

            sql = ( "SELECT "
                    "   transcript.gene_id, "
                    "   feature.chromosome_id, "
                    "   feature.direction, "
                    "   MIN(start) min_start, "
                    "   MAX(end) max_end "
                    "FROM feature, transcript "
                    "WHERE feature.transcript_id = transcript.id "
                    "AND feature.strain_id =  'Col_0' "
                    "AND chromosome_id =  'Chr1' "
                    "GROUP BY transcript.gene_id "
                    "LIMIT "+str(start)+", "+str(self.gene_location_chunk_size))

            results = engine.execute(sql)
            if results.rowcount == 0:
                break
            for row in results:
                
                db_session.add(GeneLocation(
                    gene_id=row["gene_id"], 
                    strain_id=strain_config["name"], 
                    chromosome_id=row["chromosome_id"], 
                    start=row["min_start"], 
                    end=row["max_end"], 
                    direction=row["direction"]
                ))

            start += self.gene_location_chunk_size

        db_session.commit()

    def commit_all(self):
        self.commit_entities_list(self.genes_to_write, "Genes")
        self.commit_entities_list(self.transcripts_to_write, "Transcripts")
        self.commit_entities_list(self.features_to_write, "Features")

        self.genes_to_write = []
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

            feature_type = feature_row[2]
            attribs = feature_row[8].strip()

            # This causes bugs.
            # if feature_type == "gene": # Handle gene entries
                # gene_id = attribs.split(";")[0].split(":")[1] # grab the gene ID - we'll want this for later

            new_gene_id = self.find_attribs_value("ID=Gene", attribs)
            if new_gene_id != None:

                # only deal with proper genes. setting gene_id to None means nothing else will be processed.
                # so it will essentially skip non-"gene" entries.
                if feature_type != "gene":
                    gene_id = None
                    continue

                gene_id = new_gene_id

                # add the Gene entry - if it hasn't been already
                if gene_id not in self.genes_seen: 
                    gene = Gene(gene_id)
                    self.genes_to_write.append(gene)
                    self.genes_seen[gene_id] = gene
            
            elif gene_id != None : # Handle transcript entries - if the gene is legit
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
class NucleotideMeasurementHydrator():

    def hydrate(self, experiment_config):

        # Add the experiment
        experiment = Experiment(
            id=experiment_config["experiment_id"],
            type=experiment_config["type"],
            description=experiment_config["description"]
        )
        db_session.add(experiment)
        db_session.commit() # insert the experiment into the DB.

        print("Inserting data from ["+experiment_config["nucleotides_filepath"]+"] ...")

        transcript_ids = get_inserted_transcript_ids()

        # Open the DMS reactivities file. These are normalised already.
        with open(experiment_config["nucleotides_filepath"], "r") as input_file:
            for line in input_file:
        
                bits = line.strip().split("\t")
                transcript_id = bits[0]

                # skip transcripts not already in DB
                if transcript_id not in transcript_ids:
                    continue

                if len(bits) <= 1: # no measurements present
                    continue

                count_strs = bits[1:]

                # go through reactivity entries, adding each to the database.
                position = 0
                for count_str in count_strs:
                    position += 1
                    if (count_str != "NA"): # skip adding "NA" entries.
                        obj = NucleotideMeasurement(
                            experiment_id=experiment_config["experiment_id"],
                            strain_id=settings.reference_strain_id, 
                            transcript_id=transcript_id, 
                            position=position, 
                            measurement=float(count_str)
                        )
                        db_session.add(obj)        
                db_session.commit() # insert all the reactivity measurement rows into the DB
                print("Added ["+transcript_id+"] ("+str(position)+" positions)")

        input_file.close()

# fetch all of the transcript IDs from the database, store them in a set to check against.
def get_inserted_transcript_ids(): 
    
    sql = ("SELECT id FROM transcript ORDER BY id ASC")
    results = engine.execute(sql)
    transcript_ids = set()
    for result in results:
        transcript_ids.add(result["id"])

    return transcript_ids

# Loads coverage data from a single file into the database.
class CoverageHydrator():
    def hydrate(self, experiment_config):

        transcript_ids = get_inserted_transcript_ids()

        coverage_filepath = experiment_config["coverage_filepath"]

        if not os.path.isfile(coverage_filepath):
            print("WARNING: skipped import of missing ["+coverage_filepath+"]")
            return

        with open(coverage_filepath) as coverage_file:
            for coverage_line in coverage_file:
                (transcript_id, coverage) = coverage_line.strip().split("\t")

                # skip transcripts not already in DB
                if transcript_id not in transcript_ids:
                    continue

                obj = TranscriptCoverage(
                    experiment_id=experiment_config["experiment_id"],
                    strain_id=settings.reference_strain_id,
                    transcript_id=transcript_id,
                    measurement=coverage
                )
                db_session.add(obj)

        db_session.commit()

class StructureHydrator():

    def hydrate(self, experiment_config):

        # Add the new experiment row to the DB
        experiment = Experiment(
            id=experiment_config["experiment_id"],
            type=experiment_config["type"],
            description=experiment_config["description"]
        )
        db_session.add(experiment)
        db_session.commit() # insert the experiment into the DB.

        print("Importing structures for ["+experiment_config["description"]+"]")

        transcript_ids = get_inserted_transcript_ids()
        for transcript_id in transcript_ids:

            structure_filepath = \
                experiment_config["sauce_filepath"] + \
                "/"+transcript_id+experiment_config["sauce_ext"]

            if not os.path.isfile(structure_filepath):
                print("["+structure_filepath+"] skipped")
            else:
                print("["+structure_filepath+"] found")
                self.parse_ct(structure_filepath, transcript_id, experiment_config)

    def parse_ct(self, ct_filepath, transcript_id, experiment_config):
        n_structs = 0
        with open(ct_filepath) as ct_file:
            for line in ct_file:
                # if it's an energy line, we're looking at a brand new structure
                if "ENERGY" in line:
                    # Parse the energy out using regex
                    search = re.search('ENERGY = (-[0-9\.]+)', line)
                    energy = search.group(1)

                    # Insert the new structure row
                    structure = Structure(
                        experiment_id=experiment_config["experiment_id"],
                        strain_id=experiment_config["strain_id"],
                        transcript_id=transcript_id,
                        energy=energy
                    )
                    db_session.add(structure)

                    # insert the experiment into the DB. can now access ID
                    db_session.commit() 
                    n_structs += 1

                else:

                    # the .ct format is a bit annoying because it's not tab delimited.
                    # instead it's delimited by variable numbers of spaces.

                    # calling split() with no parameter makes it split on any length 
                    # of whitespace - i.e. so that each element is 1 word
                    bits = line.strip().split()
                    from_pos = bits[0]
                    to_pos = bits[4]
                    letter = bits[1]

                    structure_position = StructurePosition(
                        structure_id=structure.id,
                        position=from_pos,
                        paired_to_position=to_pos,
                        letter=letter
                    )
                    db_session.add(structure_position)

        db_session.commit() # insert remaining data into DB
        print ("["+str(n_structs)+"] structures added")

# Carries out PCA using structures
class PcaHydrator():
    def hydrate(self, experiment_config):

        transcript_structures = {}

        # Get all transcript IDs for which there are structures
        results = db_session \
            .query(Structure.transcript_id) \
            .filter(
                Structure.experiment_id==experiment_config["experiment_id"],
                Structure.strain_id==experiment_config["strain_id"]
            ) \
            .distinct() \
            .all()

        for result in results:
            
            transcript_id = result[0]
            self.process_transcript_id(experiment_config, transcript_id)
    
    def process_transcript_id(self, experiment_config, transcript_id):

        # this is an implicit join - no need to use join() here.
        results = db_session \
            .query(Structure, StructurePosition) \
            .filter(
                StructurePosition.structure_id==Structure.id,
                Structure.experiment_id==experiment_config["experiment_id"],
                Structure.strain_id==experiment_config["strain_id"],
                Structure.transcript_id==transcript_id
            ) \
            .order_by(Structure.id, StructurePosition.position) \
            .all()

        # Get all the structures and structure_position elements for this transcript_id
        # Also keep a store of the structure objects
        structure_vecs = {}
        structures = {}
        for structure, structure_position in results:
            if structure.id not in structure_vecs:
                structure_vecs[structure.id] = []
            structure_vecs[structure.id].append(1 if structure_position.paired_to_position != 0 else 0)
            structures[structure.id] = structure

        # Do PCA using structure vectors
        pca_results = self.do_pca(structure_vecs)

        # Add the PC data to the DB
        for structure_id in structures:
            structure = structures[structure_id]
            structure.pc_1 = float(pca_results[structure.id][0])
            structure.pc_2 = float(pca_results[structure.id][1])
            db_session.add(structure)
        db_session.commit()

        
    def do_pca(self, structure_vecs):

        data = list(structure_vecs.values())

        # Do PCA.
        # Results always listed in the order that they were added.
        pca = decomposition.PCA(n_components=2)
        pca.fit(data)
        results = pca.transform(data)

        # Rearrange the data so that it is keyed by structure ID
        out = {}
        i = 0
        for structure_id in structure_vecs:
            out[structure_id] = list(results[i])
            i += 1

        return out







