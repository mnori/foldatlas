# from app import db

# from flask.ext.sqlalchemy import SQLAlchemy
# app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:vagrant@127.0.0.1/rnabrowser?charset=utf8&use_unicode=0'
# db = SQLAlchemy(app)

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from Bio import SeqIO

import settings

engine = create_engine(settings.database_uri, convert_unicode=True)
db_session = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))
Base = declarative_base()
Base.query = db_session.query_property()

import models

def reset_db():
    try:
        Base.metadata.drop_all(bind=engine)
        Base.metadata.create_all(bind=engine)
        DBHydrator().hydrate()

    except Exception as e: # catch the exception so we can display a nicely formatted error message
        print(str(e).replace("\\n", "\n").replace("\\t", "\t"))
        raise e

# Fills the database with sequence data, by parsing genome sequence .fa and annotation .gff files
class DBHydrator():

    # how many genes to process before committing rows to the database.
    chunk_size = 2500 

    genes_to_write = []
    transcripts_to_write = []
    sequences_to_write = []
    transcript_sequences_to_write = []
    transcript_sequence_features_to_write = []

    # for duplicate transcript ID detection
    transcript_ids_seen_this_strain = set()

    genes_seen = {}
    transcripts_seen = {}

    sequence_id = 0
    gene_limit = None

    # Use the genome sequence and annotation files to populate the database.
    def hydrate(self):
        for strain in settings.strains:
            self.hydrate_strain(strain)

    def hydrate_strain(self, strain_config):
        from models import Strain

        self.sequences_to_write = []
        self.transcript_sequences_to_write = []
        feature_rows = []

        chr_data = None

        # first, create the strain's entry
        strain = Strain(id=strain_config["name"], description=strain_config["description"])
        db_session.add(strain)
        db_session.commit()

        genes_added = 0

        # open the annotation file and go through it line by line
        with open(settings.genomes_sauce_folder+"/"+strain_config["annotation_filename"]) as gff_file:
            for gff_line in gff_file:
                if gff_line[0] == "#": # ignore comments
                    continue

                bits = gff_line.split("\t")
                feature_type = bits[2]
                if feature_type == "gene":
                    
                    if len(feature_rows) > 0: # this is needed to stop it going wrong at the beginning
                        self.hydrate_gene(strain, feature_rows, chr_data)
                        genes_added += 1
                        feature_rows = []
                        if genes_added % 100 == 0:
                            if self.gene_limit != None and genes_added >= self.gene_limit:
                                break
                            print (str(genes_added)+" genes processed")

                        if genes_added % self.chunk_size == 0: # commit at regular intervals
                            self.commit_all()
                    
                    chr_id_in = bits[0]
                    if chr_data == None or chr_id_in != chr_data["id"]: # detected a new chromsome - let's commit everything
                        chr_data = {
                            "id": chr_id_in,
                            "seq": self.fetch_chr_seq(strain_config, chr_id_in)
                        }

                feature_rows.append(bits)

        # gotta add that last entry
        if (len(feature_rows) > 0):
            self.hydrate_gene(strain, feature_rows, chr_data)
            self.commit_all()
            genes_added += 1

        # self.commit_all()

        # best to commit everything right at the end - it's a lot faster that way

        # first commit the genes

        # for gene_id in self.genes_to_write:
        #     db_session.add(self.genes_to_write[gene_id])
        # db_session.commit()

        # # now do the transcripts - we need the genes there already because foreign key constraints
        # for transcript_id in self.transcripts_to_write:
        #     db_session.add(self.transcripts_to_write[transcript_id])
        # db_session.commit()

        # # add strain-specific Sequence objects
        # for sequence_id in self.sequences_to_write:
        #     db_session.add(self.sequences_to_write[sequence_id])
        # db_session.commit()

        # do the sequences
        print (str(genes_added)+" genes added total")

    def commit_all(self):
        self.commit_entities_list(self.genes_to_write, "Genes")
        self.commit_entities_list(self.transcripts_to_write, "Transcripts")
        self.commit_entities_list(self.sequences_to_write, "Sequences")
        self.commit_entities_list(self.transcript_sequences_to_write, "TranscriptSequences")
        self.commit_entities_list(self.transcript_sequence_features_to_write, "TranscriptSequenceFeatures")

        self.genes_to_write = []
        self.transcripts_to_write = []
        self.sequences_to_write = []
        self.transcript_sequences_to_write = []
        self.transcript_sequence_features_to_write = []


    def commit_entities_list(self, entities, label):
        print("Committing "+label+"...")
        for entity in entities:
            db_session.add(entity)
        db_session.commit()
        print("...done.")

    def hydrate_gene(self, strain, feature_rows, chr_data):

        from models import Gene, Sequence, Transcript, TranscriptSequence, TranscriptSequenceFeature

        features = {}
        sequence = None
        transcript = None

        # first: try to find CDS related features (UTRs and CDS annotations)
        for feature_row in feature_rows:

            feature_type = feature_row[2]
            start = int(feature_row[3])
            end = int(feature_row[4])
            attribs = feature_row[8].strip()

            if feature_type == "gene": # Handle gene entries

                gene_id = attribs.split(";")[0].split(":")[1] # grab the gene ID - we'll want this for later
                chr_id = feature_row[0]
                
                sequence = str(chr_data["seq"].seq[start - 1:end])

                # create sequence object with the right coords and metadata.
                # TODO also create an intergenic sequence, concatenating gene IDs together.
                # Do that somewhere else - at the end most likely

                self.sequence_id += 1
                sequence = Sequence(
                    id=self.sequence_id,
                    strain_id=strain.id,
                    gene_id=gene_id,
                    sequence=sequence,
                    chromosome=chr_data["id"],
                    start=start,
                    end=end
                )
                self.sequences_to_write.append(sequence)

                # add the Gene entry - if it hasn't been already
                if sequence.gene_id not in self.genes_seen: 
                    gene = Gene(sequence.gene_id)
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

                    # add the TranscriptSequence entry
                    transcript_sequence = TranscriptSequence(
                        transcript_id=transcript_id,
                        sequence_id=self.sequence_id,
                        start=(start - sequence.start),
                        end=(end - sequence.start),
                        direction=("forward" if feature_row[6] == "+" else "reverse")
                    )
                    self.transcript_sequences_to_write.append(transcript_sequence)

                else: # Handle transcript feature entries
                    transcript_id = self.find_attribs_value("Parent=Transcript", attribs)

                    if transcript_id != None: # it's a transcript feature entry
                        # put a filter here? some elements are not worth storing?
                        self.transcript_sequence_features_to_write.append(TranscriptSequenceFeature(
                            sequence_id=sequence.id,
                            transcript_id=transcript.id,
                            start=(start - transcript_sequence.start),
                            end=(end - transcript_sequence.start),
                            type_id=feature_row[2]
                        ))

                    else:
                        pass # this happens for pseudogenes and TEs - do we care about those?

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

    def fetch_chr_seq(self, strain_config, chr_id):
        print("Fetching ["+chr_id+"]...")
        filepath = settings.genomes_sauce_folder+"/"+strain_config["sequence_filename"]
        for record in SeqIO.parse(filepath, "fasta"):
            if record.id == chr_id:
                print("...done.")
                return record
