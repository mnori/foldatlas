from flask import Flask, render_template, request
from sys import argv
from controllers import GenomeBrowser

import settings
import database

app = Flask(__name__)

@app.after_request
def after_request(response):
	response.headers.add('Access-Control-Allow-Origin', '*')
	return response

@app.route("/")
def index():
	browser = GenomeBrowser()
	chromosomes = browser.get_chromosomes()
	return render_template("index.html", settings=settings, chromosomes=chromosomes)

@app.route("/ajax/genome_browser/genes")
def get_genes():
	browser = GenomeBrowser()
	return browser.get_genes(request)

@app.route("/ajax/genome-browser/transcripts")
def get_transcripts():
	browser = GenomeBrowser()
	return browser.get_transcripts(request)

@app.route("/ajax/transcript/<transcript_id>") # add transcript ID to the URL
def get_transcript(transcript_id):
	buf = "transcript_id: ["+transcript_id+"]"
	return buf
	

if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`
	if (len(argv) > 1 and argv[1] == "hydratedb"):
		# reset the database
		database.hydrate_db()
		
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)
