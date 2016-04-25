#!/bin/bash
# Bootstrap for vagrant browser server
# @author Matthew Norris

DBPASSWD="vagrant"

. /vagrant/bootstrap/functions.sh

# Installation stuff goes here.
function install() {

	pretty_print "PROVISIONING"

	apt-get update

	# Copy handy bash aliases to home folder. Must use explicit home folder path, otherwise 
	# it'll copy to super user's path instead of vagrant's
	cp /vagrant/bootstrap/.bash_aliases /home/vagrant/.bash_aliases

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing apache2" 
	apt-get install -y apache2
	a2enmod proxy
	a2enmod proxy_http 
	a2enmod headers 

	# copy apache config file
	cp /vagrant/bootstrap/000-default_vagrant.conf /etc/apache2/sites-available/000-default.conf

	# setup static folder
	ln -s /vagrant/static /var/www/static

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing MySQL specific packages and settings"
	export DEBIAN_FRONTEND=noninteractive

	# this bypasses the crappy GUI-based install
	echo "mysql-server mysql-server/root_password password $DBPASSWD" | debconf-set-selections
	echo "mysql-server mysql-server/root_password_again password $DBPASSWD" | debconf-set-selections
	echo "phpmyadmin phpmyadmin/dbconfig-install boolean true" | debconf-set-selections
	echo "phpmyadmin phpmyadmin/app-password-confirm password $DBPASSWD" | debconf-set-selections
	echo "phpmyadmin phpmyadmin/mysql/admin-pass password $DBPASSWD" | debconf-set-selections
	echo "phpmyadmin phpmyadmin/mysql/app-pass password $DBPASSWD" | debconf-set-selections
	echo "phpmyadmin phpmyadmin/reconfigure-webserver multiselect none" | debconf-set-selections
	apt-get -y install mysql-server phpmyadmin > /dev/null
	a2disconf phpmyadmin # switch off the PMA conf - not needed
	service apache2 restart
	echo "create database foldatlas" | mysql -u root -p$DBPASSWD

	# custom config - needed for DB export to work properly
	cp /vagrant/bootstrap/my.cnf /etc/mysql/my.cnf 

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing git"
	apt-get install -y git
	
	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing pip"
	apt-get install -y python3-pip
	
	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing sklearn"
	
	# lots of scikit dependencies. use apt-get - much faster to grab the binaries
	# than compile from scratch using pip3
	apt-get install -y build-essential python3-dev python3-setuptools 
	apt-get install -y libatlas-dev libatlas3gf-base
	apt-get install -y python3-numpy python3-scipy

	# must use pip3 for scikit itself
	pip3 install scikit-learn
	
	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing misc python dependencies"
	pip3 install biopython
	pip3 install Flask
	pip3 install Flask-Migrate

	# echo https://cdn.mysql.com/Downloads/Connector-Python/mysql-connector-python-2.1.3.tar.gz >> requirements.txt
	# sudo -H pip3 install -r ./requirements.txt

	# mysql-connector-python is a pain for some reason
	# see https://stackoverflow.com/questions/31748278/how-do-you-install-mysql-connector-python-development-version-through-pip
	git clone https://github.com/mysql/mysql-connector-python.git
	cd mysql-connector-python
	python3 ./setup.py build
	python3 ./setup.py install
	pip3 install Flask-SQLAlchemy
	# pip3 install forgi # autotranslated, for now

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing Clustalw2"
	cd /tmp/
	wget http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz >&/dev/null
	tar xvzf clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
	cd clustalw-2.1-linux-x86_64-libcppstatic/
	cp clustalw2 /usr/local/bin/clustalw2

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Installing ViennaRNA"
	apt-add-repository -y ppa:j-4/vienna-rna
	apt-get -y update
	apt-get -y install vienna-rna

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	if [ ! -f "/vagrant/static/downloads/foldatlas.sql.tar.gz" ]
		then
		dl_sauce
		hydrate_db # create the DB from the raw data
		export_db # export the database so it can be quickly imported again

	else
		import_db # DB was previous exported. Import it since this is way quicker.

	fi

	# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	pretty_print "Provisioning complete"
}

# get the party started
install
