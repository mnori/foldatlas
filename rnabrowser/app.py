from flask import Flask
from flask import render_template
from sys import argv

import settings

app = Flask(__name__)

@app.route("/")
def hello():
	print(settings.base_path+"/templates/hello.html")
	return render_template("hello.html", message="Hello world!")

if __name__ == "__main__": 
	# if we're in here, we're using `python3 app.py [blah...]`

	if (len(argv) > 1 and argv[1] == "resetdb"):
		# reset the database
		from database import reset_db
		reset_db()
		
	else:
		# dev server: get the party started
		app.run(host='0.0.0.0', debug=True)


