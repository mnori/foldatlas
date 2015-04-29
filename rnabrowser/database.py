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
	    hydrate_db()

    except Exception as e: # catch the exception so we can display a nicely formatted error message
    	print(str(e).replace("\\n", "\n").replace("\\t", "\t"))
    	raise e


# Use the genome sequence and annotation files to populate the database.
def hydrate_db():
	for strain in settings.strains:
		hydrate_strain(strain)

def hydrate_strain(strain_config):

	from models import Strain, Transcript

	chromosome_seq = None

	# first, create the strain's entry
	strain = Strain(name=strain_config["name"], description=strain_config["description"])
	db_session.add(strain)
	db_session.commit()

	# open the annotation file and go through it line by line
	with open(settings.genomes_sauce_folder+"/"+strain_config["annotation_filename"]) as gff_file:
		for gff_line in gff_file:
			if gff_line[0] == "#": # ignore comments
				continue

			# if we got this far, it's a line of annotation info
			bits = gff_line.split("\t")

			chr_id_in = bits[0]
			feature_type = bits[2]
			start = bits[3]
			end = bits[4]
			direction = bits[5]
			attribs = bits[7]

			# fetch the chromosome sequence if requiried
			if chromosome_seq == None or chr_id_in != chromosome_seq.id:
				chromosome_seq = fetch_chromosome_seq(strain_config, chr_id_in)

			

			# # if it's a gene, add new gene entry
			# if feature_type == "gene":
			# 	# parse out the gene_id
			# 	# do we need a new Transcript for this gene?
			# 	Transcript.query.filter_by(id='peter').first()
			# 	sequence = models.Sequence()





def fetch_chromosome_seq(strain_config, chr_id):
	filepath = settings.genomes_sauce_folder+"/"+strain_config["sequence_filename"]
	for record in SeqIO.parse(filepath, "fasta"):
		if record.id == chr_id:
			return record

# def reset_db():
# 	db.drop_all()
# 	db.create_all()

# create a new user class for testing purposes
# class User(db.Model):
#     id = db.Column(db.Integer, primary_key=True)
#     name = db.Column(db.String(80), unique=True)

#     # Constructor
#     def __init__(self, name):
#         self.name = name

#     # Return string representing the object
#     def __repr__(self):
#         return '<Test %r>' % self.username

