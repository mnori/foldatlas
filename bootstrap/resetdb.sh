# @author Matthew Norris
# Resets the database, without you having to fuck about with django's 
# shitty migration commands.

# Bit silly that django doesn't have it's own command for this

projectroot="/vagrant/rnabrowser"
dbname="rnabrowser"
appname="rnabrowserapp"

cd $projectroot

# wipe the migrations directory clean
rm -rf $appname/migrations 
mkdir $appname/migrations

# wipe migrations from DB
echo "USE $dbname; DELETE FROM django_migrations WHERE app = '$appname'" | mysql -u root

# delete all old app-specific tables
mysql -u root -D $dbname -e "show tables" -s | 
  egrep "^$appname" | 
  xargs -I "@@" echo "SET foreign_key_checks = 0; DROP TABLE @@; " | 
  mysql -u root -D $dbname

# reset django shiznits
python3 manage.py flush --noinput
python3 manage.py makemigrations $appname
python3 manage.py migrate $appname
