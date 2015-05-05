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
	from database import Feature, Transcript, db_session;
	from sqlalchemy import and_
	import json

	out = []

	chromosome_id = "Chr1"
	start = 0
	end = 100000


	sql = 	("SELECT *, MIN(start) min_start, MAX(end) max_end "
			 "FROM feature "
			 "WHERE chromosome_id = '"+chromosome_id+"' "
			 "AND start > '"+str(start)+"' "
			 "AND end < '"+str(end)+"' "
			 "GROUP BY transcript_id")

	# demonstrates how genes would be done ###########################
	results = database.engine.execute(sql)
	for result in results:
		out.append({
			"Parent": result.transcript_id, # actually use gene ID instead
			"feature_type": "transcript", # without this, it won't draw
			"logic_name": "ensembl_havana", # does not work without this label
			"start": result.min_start,
			"end": result.max_end,
			"id": result.transcript_id,
			"strand": 1, # whether it is + or -??
		})

	

	# /demonstrates how genes would be done ###########################

	results =  db_session.query(Feature).filter(and_(Feature.start >= start, Feature.end <= end)).all() 
	for feature in results:

		print(feature)

		out.append({
			"Parent": feature.transcript_id, # actually use gene ID instead
			"feature_type": feature.type_id, # without this, it won't draw the gene
			"logic_name": "ensembl_havana", # does not work without this label
			"start": feature.start,
			"end": feature.end,
			"id": feature.transcript_id+"-"+str(feature.id),
			"strand": 1, # whether it is + or -??
		})


	buf = json.dumps(out)
	return buf

if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`

	if (len(argv) > 1 and argv[1] == "resetdb"):
		# reset the database
		database.reset_db()
		
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)



