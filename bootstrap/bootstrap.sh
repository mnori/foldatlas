# Bootstrap for vagrant browser server
# @author Matthew Norris

# Prints aesthetically pleasing messages into the terminal.
function pretty_print() {
	printf "\n-=-=-=-=[ $1 ]=-=-=-=-\n"
    len=${#1} ch='-'
    padding=$(printf '%*s' "$len" | tr ' ' "$ch")
    printf "          $padding\n\n"
}

# Fetch a file and save it.
function fetch_raw() {
	base_sauce_url=$1
	base_filename=$2
	sauce_url="$base_sauce_url$base_filename"

	echo "Processing [$base_filename]..."

	if [ ! -f "$base_filename" ] # if file does not exist...
		then
		echo "Fetching..."
		curl "$sauce_url" -o "$base_filename" >&/dev/null
		echo "...Done."
	else
		echo "...Already exists!"
	fi
	
}

# Fetch a .bz2 file and extract it
function fetch_extract_gff3() {
	base_sauce_url=$1
	base_filename=$2
	sauce_url="$base_sauce_url$base_filename"

	echo "Processing [$base_filename]..."
	if [ ! -f "$base_filename.gff3" ]
		then
		echo "Fetching..."
		curl "$sauce_url.gff3.bz2" -o "$base_filename.gff3.bz2" >&/dev/null # fetch compressed file
		bunzip2 "$base_filename.gff3.bz2" # decompress to get the .gff3 file
		echo "...Done."
	else
		echo "...Already exists!"
	fi
}

# Installation stuff goes here.
function install() {

	pretty_print "PROVISIONING"

	cp /vagrant/bootstrap/.bash_aliases ~/.bash_aliases

	pretty_print "Installing pip"
	apt-get update
	apt-get install -y python3-pip

	# Django old install
	# pretty_print "Installing Django"
	# pip3 install Django

	pretty_print "Installing BioPython"
	pip3 install biopython

	pretty_print "Installing MySQL"
	export DEBIAN_FRONTEND=noninteractive
	apt-get -q -y install mysql-server

	echo "create database rnabrowser" | mysql -u root
	apt-get install -y git
	apt-get install -y libmysqlclient-dev
	pip3 install mysqlclient

	# Install project
	cd /vagrant/rnabrowser 
	python3.3 setup.py develop # install all dependencies

	# easy_install-3.4 pyramid waitress

	#  python3.4 setup.py develop ??

	# Install

	# Django migrate shiz
	# pretty_print "Migrating Database"
	# cd /vagrant/rnabrowser
	# rm rnabrowserapp/migrations/* # otherwise old migrations might bork it
	# python3 manage.py migrate # django setup migration
	# python3 manage.py makemigrations rnabrowserapp
	# python3 manage.py migrate rnabrowserapp
}

# Grabs lots of genome data files. These will be parsed and used to seed the SNP database.
function dl_sauce() {
	# TODO
	# grab raw genome data from network if available. 
	# otherwise, download the data from its origin

	pretty_print "Grabbing Sauce Data"

	cd /vagrant/sauce_data

	# Grab sequence annotation files, relative to both strain-of-interest and reference, for each strain
	urlbase="http://mus.well.ox.ac.uk/19genomes/annotations/consolidated_annotation_9.4.2011/gene_models/"

	# This is the reference sequence annotation
	fetch_extract_gff3 $urlbase "consolidated_annotation.Col_0.Col_0"
	fetch_extract_gff3 $urlbase "consolidated_annotation.Col_0"

	# fetch_extract_gff3 $urlbase "consolidated_annotation.Bur_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Bur_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Can_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Can_0"
	
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ct_1.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ct_1"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Edi_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Edi_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Hi_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Hi_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Kn_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Kn_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ler_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ler_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Mt_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Mt_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.No_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.No_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Oy_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Oy_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Po_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Po_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Rsch_4.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Rsch_4"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Sf_2.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Sf_2"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Tsu_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Tsu_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wil_2.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wil_2"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ws_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ws_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wu_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wu_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Zu_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Zu_0"


	# Fetch the chromosomal TAIR10 reference sequence. We'll get the RNAs out using the .gff3 annotation
	urlbase="ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/"
	fetch_raw $urlbase "TAIR10_chr1.fas"
	fetch_raw $urlbase "TAIR10_chr2.fas"
	fetch_raw $urlbase "TAIR10_chr3.fas"
	fetch_raw $urlbase "TAIR10_chr4.fas"
	fetch_raw $urlbase "TAIR10_chr5.fas"
	fetch_raw $urlbase "TAIR10_chrC.fas"
	fetch_raw $urlbase "TAIR10_chrM.fas"

	# Grab chromosomal sequences for each of the 18 other strains
	# commented out to make it quicker for testing

	# urlbase="http://mus.well.ox.ac.uk/19genomes/fasta/MASKED/"
	# fetch_raw $urlbase "bur_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "can_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "ct_1.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "edi_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "hi_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "kn_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "ler_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "mt_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "no_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "oy_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "po_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "rsch_4.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "sf_2.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "tsu_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "wil_2.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "ws_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "wu_0.v7.PR_in_lowercase.fas"
	# fetch_raw $urlbase "zu_0.v7.PR_in_lowercase.fas"
}

# get the party started
install
dl_sauce

pretty_print "Provisioning Complete"

