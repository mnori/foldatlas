dbname="rnabrowser"
appname="rnabrowserapp"
projectroot="/vagrant/rnabrowser"

cd $projectroot

# wipe the migrations directory clean
rm -rf $appname/migrations 
mkdir $appname/migrations

# also wipe migrations from DB
echo "USE $dbname; DELETE FROM django_migrations WHERE app = '$appname'" | mysql -u root

# delete all old tables
mysql -u root -D $dbname -e "show tables" -s | 
  egrep "^$appname" | 
  xargs -I "@@" echo "SET foreign_key_checks = 0; DROP TABLE @@; " | 
  mysql -u root -D $dbname

# reset django shiznits
python3 manage.py flush --noinput
python3 manage.py makemigrations $appname
python3 manage.py migrate $appname
