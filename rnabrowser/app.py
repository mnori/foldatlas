from flask import Flask
from flaskext.mysql import MySQL

app = Flask(__name__)

mysql = MySQL()
app.config['MYSQL_DATABASE_USER'] = 'root'
app.config['MYSQL_DATABASE_PASSWORD'] = ''
app.config['MYSQL_DATABASE_DB'] = 'rnabrowser'
app.config['MYSQL_DATABASE_HOST'] = 'localhost'
mysql.init_app(app)

@app.route("/")
def hello():
	cursor = mysql.connect().cursor()
	cursor.execute("SHOW TABLES")
	data = cursor.fetchone()

	return str(data)

if __name__ == "__main__": # if we called the script directly... (run server from command line)
    app.run(host='0.0.0.0')
