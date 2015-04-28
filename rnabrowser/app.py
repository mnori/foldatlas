from flask import Flask
from sys import argv
app = Flask(__name__)

@app.route("/")
def hello():
	return "Hello World test"

if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`

	if (argv[1] == "resetdb"):
		# reset the database
		from database import reset_db
		reset_db()
		
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)


