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

	results = db_session.query(Feature).limit(100).all() 
	
	# convert results into JSON that the thingy understands

	buf = ""
	for result in results:
		buf += result.chromosome_id+"\t"
		buf += ".\t" # source
		buf += result.type_id+"\t"
		buf += str(result.start)+"\t"
		buf += str(result.end)+"\t"
		buf += ".\t" # score
		buf += "+\t" # foward or reverse
		buf += ".\t" # ???
		buf += "Parent=Transcript:"+result.transcript_id+"\n" # attribs

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


