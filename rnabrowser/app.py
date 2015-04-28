from flask import Flask
app = Flask(__name__)

from flask.ext.sqlalchemy import SQLAlchemy

app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:@127.0.0.1/rnabrowser?charset=utf8&use_unicode=0'
db = SQLAlchemy(app)

# create a new user class for testing purposes
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(80), unique=True)

    # Constructor
    def __init__(self, name):
        self.name = name

    # Return string representing the object
    def __repr__(self):
        return '<Test %r>' % self.username

@app.route("/")
def hello():
	return "Hello World test"

if __name__ == "__main__": 
	# if we're in here, it means we're running the built in development web server.
    app.run(host='0.0.0.0', debug=True)
