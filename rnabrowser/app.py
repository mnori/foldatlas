from flask import Flask
app = Flask(__name__)

from flask.ext.sqlalchemy import SQLAlchemy
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:@127.0.0.1/rnabrowser?charset=utf8&use_unicode=0'
db = SQLAlchemy(app)

@app.route("/")
def hello():
	return "Hello World test"

if __name__ == "__main__": 
	# if we're in here, it means we're running the built in development web server.
    app.run(host='0.0.0.0', debug=True)
