function pretty_print() {
    printf "\n-=-=-=-=[ $1 ]=-=-=-=-\n"
    len=${#1} ch='-'
    padding=$(printf '%*s' "$len" | tr ' ' "$ch")
    printf "          $padding\n\n"
}

pretty_print "PROVISIONING"
pretty_print "Installing pip"
apt-get update
apt-get install -y python3-pip


pretty_print "Installing Django"
pip3 install Django

pretty_print "Installing BioPython"
pip3 install biopython

pretty_print "Installing MySQL"
export DEBIAN_FRONTEND=noninteractive
apt-get -q -y install mysql-server

echo "create database rnabrowser" | mysql -u root
apt-get install -y git
apt-get install -y libmysqlclient-dev
pip3 install mysqlclient

pretty_print "Migrating Database"
cd /vagrant/rnabrowser
rm rnabrowserapp/migrations/* # otherwise old migrations might bork it
python3 manage.py migrate
python3 manage.py makemigrations rnabrowserapp
python3 manage.py migrate rnabrowserapp

pretty_print "Provisioning Complete"
