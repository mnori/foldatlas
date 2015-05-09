# from app import db

# from flask.ext.sqlalchemy import SQLAlchemy
# app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:vagrant@127.0.0.1/rnabrowser?charset=utf8&use_unicode=0'
# db = SQLAlchemy(app)

from sqlalchemy import create_engine, and_
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from Bio import SeqIO
from Bio.Seq import Seq
import settings, os

engine = create_engine(settings.database_uri, convert_unicode=True)
db_session = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))
Base = declarative_base()
Base.query = db_session.query_property()

import models

from models import Strain, Gene, Transcript, Feature

def hydrate_db():
    try:
        print("Rebuilding schema...")
        Base.metadata.drop_all(bind=engine)
        Base.metadata.create_all(bind=engine)
        print("...done.")
        SequenceHydrator().hydrate() # add the annotations
        TranscriptAligner().align() # make the alignments

    except Exception as e: # catch the exception so we can display a nicely formatted error message
        print(str(e).replace("\\n", "\n").replace("\\t", "\t"))
        raise e

# Parses genome sequence .fa and annotation .gff3 files into the database.
class SequenceHydrator():

    # how many genes to process before committing rows to the database.
    gene_chunk_size = 2500 

    # how many bytes of sequence data to use per chromosome writing chunk
    seq_chunk_size = 1048576 # i.e. 1 MB

    genes_to_write = []
    transcripts_to_write = []
    features_to_write = []

    genes_seen = {}
    transcripts_seen = {}

    # for duplicate transcript ID detection
    transcript_ids_seen_this_strain = set()

    sequence_id = 0

    # how many genes to process
    gene_limit = 100

    # how much chromosome sequence to add, in bp
    bp_limit = None

    # Use the genome sequence and annotation files to populate the database.
    def hydrate(self):
        for strain in settings.strains:
            self.hydrate_strain(strain)

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
    # Therefore a LOAD DATA INFILE strategy is used.
    def hydrate_chrosomomes(self, strain_config):

        print("Adding chromosomes...")
        filepath = settings.genomes_sauce_folder+"/"+strain_config["sequence_filename"]

        for record in SeqIO.parse(filepath, "fasta"): # loop through chromosomes

            # Save the chromosome data to a text file
            chr_id = record.id

            if (chr_id in settings.ignored_chromosomes):
                continue

            seq_str = str(record.seq)
            temp_file = open(settings.temp_db_file, "w")
            temp_file.write(strain_config["name"]+"\t"+chr_id+"\t"+seq_str)
            temp_file.close()

            # Load the file into the DB
            sql = """
                LOAD DATA LOCAL INFILE '/tmp/tmp.txt'
                REPLACE INTO TABLE chromosome
            """
            db_session.execute(sql)
            db_session.commit()

            # Delete the file
            os.remove(settings.temp_db_file)

            print("Added ["+chr_id+"]")


            # chr_id = record.id

            

            # seq_str = str(record.seq)
            # len_seq_str = self.bp_limit if self.bp_limit != None else len(seq_str)

            # chunk = seq_str[0:len_seq_str if len_seq_str <= self.seq_chunk_size else self.seq_chunk_size]

            # # Insert the Chromosome row
            # db_session.execute(
            #     "INSERT INTO chromosome SET"
            #     "   strain_id = '"+strain_config["name"]+"',"
            #     "   chromosome_id = '"+chr_id+"',"
            #     "   sequence = '"+chunk+"'"
            # )
            # total_test += chunk

            # pos = self.seq_chunk_size

            # # Append chunks of sequence onto the row's sequence field.
            # while True: 
            #     end = pos + self.seq_chunk_size
            #     if (end > len_seq_str):
            #         end = len_seq_str

            #     chunk = seq_str[pos:end]

            #     db_session.execute(
            #         "UPDATE chromosome SET"
            #         "   sequence = CONCAT(sequence, '"+chunk+"') "
            #         "WHERE strain_id = '"+strain_config["name"]+"' "
            #         "AND chromosome_id = '"+chr_id+"'"
            #     )
            #     total_test += chunk

            #     if end == len_seq_str:
            #         break
            #     pos += self.seq_chunk_size

            # db_session.commit()

            # print("Added ["+chr_id+"]")

            # # f = open("test.fa", "w")
            # # f.write(total_test)
            # # f.close()

            # exit()

        print("Finished adding chromosomes to ["+strain_config["name"]+"]")
        # exit()

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

        for feature_row in feature_rows: # Loop through annotation rows in the gff file, all related to the current gene

            feature_type = feature_row[2]
            attribs = feature_row[8].strip()

            if feature_type == "gene": # Handle gene entries
                gene_id = attribs.split(";")[0].split(":")[1] # grab the gene ID - we'll want this for later

                # add the Gene entry - if it hasn't been already
                if gene_id not in self.genes_seen: 
                    gene = Gene(gene_id)
                    self.genes_to_write.append(gene)
                    self.genes_seen[gene_id] = gene

            
            else: # Handle transcript entries
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
                    transcript_id = self.find_attribs_value("Parent=Transcript", attribs)
                    if transcript_id != None: # it's a transcript feature entry
                        # put a filter here? some elements are not worth storing?
                        self.features_to_write.append(Feature(
                            transcript_id=transcript_id,
                            type_id=feature_row[2],
                            strain_id=strain_id,
                            chromosome_id=feature_row[0],
                            start=feature_row[3],
                            end=feature_row[4],
                            direction="forward" if feature_row[6] == "+" else "reverse"
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
        # given the transcript ID, fetch the feature sequences in the correct order.
        sql = """
            SELECT 
                feature.strain_id, 
                SUBSTR(chromosome.sequence, feature.start, feature.end - feature.start) seq
            FROM chromosome, feature
            WHERE 
                feature.strain_id = chromosome.strain_id AND
                feature.chromosome_id = chromosome.chromosome_id AND
                feature.type_id = 'exon' AND
                feature.transcript_id = '{0}'
            ORDER BY feature.strain_id, start
        """
        sql = sql.format(transcript_id)

        # print(sql)
        # exit()

        results = engine.execute(sql)

        # let's create BioPython sequence objects here
        transcript_seqs = {}
        for row in results:
            strain_id = row["strain_id"]
            if strain_id  not in transcript_seqs:
                transcript_seqs[strain_id] = Seq("")

            transcript_seqs[strain_id] += row["seq"]

        # Does 
        # if strain_id != "Col_0":

        # if (transcript_id == "AT1G01010"):
        for strain_id in transcript_seqs:
            print("transcript_id: ["+transcript_id+"]")
            print("strain_id: ["+strain_id+"]")
            print("seq: ["+transcript_seqs[strain_id]+"]")
        

        exit()

        # concatenate the feature sequences together


        # if direction is reverse, do reverse-complement


    # Fetch all the transcript IDs from the database.
    def fetch_transcript_ids(self):
        transcript_ids = set()
        sql = "SELECT id FROM transcript"
        rows = engine.execute(sql)
        for row in rows:
            transcript_ids.add(row["id"])

        return transcript_ids


        
