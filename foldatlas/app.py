from flask import Flask, render_template, request, Response, send_from_directory

from flask_sqlalchemy import SQLAlchemy
from sys import argv
import settings
from utils import FastaExporter, FastaSplitter

app = Flask(__name__)

# these commands enable migrations, see https://flask-migrate.readthedocs.org/en/latest/
# print("Creating db object")


# def load_models():
# 	from models import Gene
# 	print("Models loaded")

# from controllers import GenomeBrowser, TranscriptView, TranscriptSearcher, CoverageSearcher, \
# 	StructureDiagramView, StructureCirclePlotView, StructureDownloader, \
# 	NucleotideMeasurementDownloader

# db_session
@app.teardown_appcontext
def close_db(error):
	from database import db_session
	db_session.close()

@app.route('/static/<path:path>/<filename>')
def send_static(path, filename):
	# probably not very secure but it'll do for now
	return send_from_directory(settings.static_path, path+"/"+filename)

@app.after_request
def after_request(response):
	response.headers.add('Access-Control-Allow-Origin', '*')
	return response

# The index should show the genome browser and also a search box.
# Maybe also some introductory text.
@app.route("/")
def index():
	from controllers import GenomeBrowser
	return render_template("index.html", settings=settings, genome_browser=GenomeBrowser(), page="home")

@app.route("/search")
def search():
	from controllers import GenomeBrowser
	return render_template("index.html", settings=settings, genome_browser=GenomeBrowser(), page="search")

@app.route("/help")
def help():
	from controllers import GenomeBrowser
	return render_template("index.html", settings=settings, genome_browser=GenomeBrowser(), page="help")

# Transcript - initialise the genome browser with custom parameters to center on the gene of interest.
# Also show the transcript's details
@app.route("/transcript/<transcript_id>")
def view_transcript(transcript_id):
	from controllers import GenomeBrowser, TranscriptView
	return render_template("index.html", 
		settings=settings, 
		genome_browser=GenomeBrowser(), 
		transcript_view=TranscriptView(transcript_id),
		page="transcript")

@app.route("/ajax/genome-browser/genes")
def get_genes_ajax():
	from controllers import GenomeBrowser
	return GenomeBrowser().get_genes(request)

@app.route("/ajax/help")
def get_help_ajax():
	return render_template("help-view.html", settings=settings)

@app.route("/ajax/genome-browser/transcripts")
def get_transcripts_ajax():
	from controllers import GenomeBrowser
	return GenomeBrowser().get_transcripts(request)

@app.route("/ajax/search-transcript/<search_string>")
def search_transcripts_ajax(search_string):
	from controllers import TranscriptSearcher
	return TranscriptSearcher().search(search_string)

@app.route("/ajax/search-coverage/<page_num>")
def search_coverage_ajax(page_num):
	from controllers import CoverageSearcher
	return render_template(
		"coverage-search.html", 
		transcript_data=CoverageSearcher().fetch_transcript_data(page_num)
	)

@app.route("/ajax/get-coverage-page-count")
def get_coverage_page_count():
	from controllers import CoverageSearcher
	return str(CoverageSearcher().fetch_page_count())

@app.route("/ajax/transcript/<transcript_id>")
def view_transcript_ajax(transcript_id):
	from controllers import TranscriptView
	return render_template("transcript-view.html", transcript_view=TranscriptView(transcript_id))

@app.route("/ajax/structure-diagram/<structure_id>")
def structure_diagram_ajax(structure_id):
	from controllers import StructureDiagramView
	return StructureDiagramView(structure_id).data_json

@app.route("/ajax/structure-circle-plot/<structure_id>")
def structure_circle_plot_ajax(structure_id):
	from controllers import StructureCirclePlotView
	return StructureCirclePlotView(structure_id).data_json

# strain ID .. should really be experiment ID
# and strain ID should only be associated with experiment ID.
# can then just use experiment IDs for everything.
@app.route("/download/structure/<transcript_id>")
def download_structure(transcript_id):
	from controllers import StructureDownloader
	buf = StructureDownloader([1, 2], transcript_id).generate()
	return Response(buf, mimetype='text/plain')

@app.route("/download/bppm/<transcript_id>")
def download_bppm(transcript_id):
	from controllers import BppmDownloader
	return Response(BppmDownloader().fetch(transcript_id), mimetype='text/plain')

@app.route("/download/measurements/<experiment_id>/<transcript_id>")
def download_measurements(experiment_id, transcript_id):
	from controllers import NucleotideMeasurementDownloader
	buf = NucleotideMeasurementDownloader(experiment_id, transcript_id).get_normalised()
	return Response(buf, mimetype='text/plain')

@app.route("/download/raw_measurements/<experiment_id>/<transcript_id>")
def download_raw_measurements(experiment_id, transcript_id):
	from controllers import NucleotideMeasurementDownloader
	buf = NucleotideMeasurementDownloader(experiment_id, transcript_id).get_raw()
	return Response(buf, mimetype='text/plain')

@app.route("/download/all")
def download_all():
	return "All data here"

if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`
	if len(argv) > 1:  
		cmd = argv[1]

		# custom commands
		if argv[1] == "grabstructures":
			import utils
			utils.grab_structures()

		elif argv[1] == "hydratedb":
			# reset the database
			import database
			database.import_db(2)

		elif argv[1] == "exportfasta":
			# export sequences into a big fasta file
			FastaExporter().export()

		elif argv[1] == "splitfasta":
			# split fasta sequences into individual ones
			FastaSplitter().split()

		elif argv[1] == "export_structure_tids":
			from importers import StructureTidsExporter
			StructureTidsExporter().export()

		else: # some other command
			from models import manager
			manager.run()
	
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)
