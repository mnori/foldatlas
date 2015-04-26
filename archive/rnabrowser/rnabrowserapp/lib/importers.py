# Classes for adding sequence data to the database
# @author Matthew Norris

from Bio import SeqIO
from django.apps import apps
from django.db import connection
from rnabrowserapp.models import Strain, Transcript, TranscriptSequence;

# Imports sequences from a FASTA file to MySQL DB
class Importer():

    def __init__(self, settings):
        self.settings = settings
        
    def do_import(self): 
        self.__truncate_tables()
        self.__import_reference()

    def __truncate_tables(self):
        print("Truncating tables...")
        cursor = connection.cursor()
        cursor.execute('SET FOREIGN_KEY_CHECKS = 0')
        cursor.execute('TRUNCATE TABLE '+str(TranscriptSequence._meta.db_table))
        cursor.execute('TRUNCATE TABLE '+str(Strain._meta.db_table))
        cursor.execute('TRUNCATE TABLE '+str(Transcript._meta.db_table))
        cursor.execute('SET FOREIGN_KEY_CHECKS = 1')
        print("...Done.")

    def __import_reference(self):
        refseq_details = self.settings.REFERENCE_SEQ

        # save the reference strain details
        strain = self.__save_strain_details(refseq_details)

        # go through sequences, save each one
        print("Saving sequences to database...")
        path = refseq_details["path"]
        handle = open(path, "rU")
        n = 0
        for record in SeqIO.parse(handle, "fasta"):
            bits = record.id.split("|")
            accession_number = bits[3]
            gi_number = bits[1]

            transcript = Transcript()
            transcript.id = accession_number
            transcript.gi_number = gi_number
            transcript.save()

            transcript_sequence = TranscriptSequence()
            transcript_sequence.strain = strain
            transcript_sequence.transcript = transcript
            transcript_sequence.sequence = str(record.seq)
            transcript_sequence.save()

            n += 1
            if n % 100 == 0:
                print("Saved "+str(n)+" sequences")
            
        handle.close()
        print("...Done. Total of "+str(n)+" sequences processed")
        # sequences = utils.load_sequences_fasta(refseq_details["path"])
        # print("n: "+len(sequences))

    def __save_strain_details(self, strain_details):
        strain = Strain()
        strain.name = strain_details["strain_name"]
        strain.description = strain_details["strain_description"]
        strain.save()
        return strain

        



