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

@app.route("/browser-data/genes")
def get_genes():
	browser = GenomeBrowser()
	return browser.get_genes(request)

@app.route("/browser-data/transcripts")
def get_transcripts():
	browser = GenomeBrowser()
	return browser.get_transcripts(request)

if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`
	if (len(argv) > 1 and argv[1] == "resetdb"):
		# reset the database
		database.reset_db()
		
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)
