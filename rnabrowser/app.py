from flask import Flask
from flask import render_template
from sys import argv

import settings
import database

app = Flask(__name__)

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route("/")
def hello():
	print(settings.base_path+"/templates/hello.html")
	return render_template("hello.html", message="Hello world!")

@app.route("/test")
def test():
	from database import Feature, db_session;
	import json

	results = db_session.query(Feature).limit(1000).all() 
	
	# Ensembl ################################
	out = []
	for result in results:
		out.append({
			"Parent": result.transcript_id, # actually use gene ID instead
			"assembly_name": "GRCh38", # this should actually be strain name
			"biotype": "protein_coding", # result.type_id, 
			"description": None,
			"start": result.start,
			"end": result.end,
			"external_name": result.transcript_id, # probs doesn't matter
			"feature_type": "transcript", # without this, it won't draw the gene
			"id": result.transcript_id,
			"logic_name": "ensembl_havana",
			"seq_region_name": result.chromosome_id,
			"source": result.strain_id,
			"strand": 1, # whether it is + or -??
			"version": 1 # very unlikely to matter
		})
	buf = json.dumps(out)
	return buf
	# /Ensembl ###############################

	# convert results into format that the thingy understands

	# buf = ""
	# for result in results:
	# 	buf += result.chromosome_id+"\t"
	# 	buf += ".\t" # source
	# 	buf += result.type_id+"\t"
	# 	buf += str(result.start)+"\t"
	# 	buf += str(result.end)+"\t"
	# 	buf += "9\t" # score
	# 	buf += "+\t" # foward or reverse
	# 	buf += ".\t" # ???
	# 	buf += "Parent=Transcript:"+result.transcript_id+"\n" # attribs

	# null+"\t"
	# buf = json.dumps([
	# 	'foo', {'bar': ('baz', None, 1.0, 2)}])

	# buf = "["

	# for result in results:
	# 	buf += "{"
	# 	buf += "id: \""+result.transcript_id+"\", "
	# 	buf += "}"

	# buf += "]"

	# line = "{"
	# 	"id: "

	# feature.id     = fields.slice(0, 5).join('|');
 #    feature.start  = parseInt(fields[3], 10);
 #    feature.end    = parseInt(fields[4], 10);
 #    feature.source = fields[1];
 #    feature.type   = fields[2];
 #    feature.score  = fields[5];
 #    feature.strand = fields[6] + '1';


	return buf


if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`

	if (len(argv) > 1 and argv[1] == "resetdb"):
		# reset the database
		database.reset_db()
		
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)


