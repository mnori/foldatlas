#!/bin/bash

function hydrate_db() {
	pretty_print "Hydrating database"
	cd /vagrant/foldatlas
	python3 app.py hydratedb
}

function export_db() {
	pretty_print "Dumping database"
	cd /vagrant/static/downloads/
	mysqldump -uroot -pvagrant --add-drop-database foldatlas > foldatlas.sql
	tar -czpf foldatlas.sql.tar.gz foldatlas.sql 
	rm foldatlas.sql
}

function import_db() {
	pretty_print "Importing database"
	cd /vagrant/static/downloads/
	tar xvzf foldatlas.sql.tar.gz

	# stops mysql crashing out on import
	echo "SET GLOBAL max_allowed_packet=1073741824;" | mysql -u root -pvagrant

	# does the actual import
	mysql -uroot -pvagrant foldatlas < foldatlas.sql
	rm foldatlas.sql
}

# Grabs lots of genome data files. These will be parsed and used to seed the SNP database.
function dl_sauce() {
	pretty_print "Grabbing sauce data"
	# TODO
	# grab raw genome data from network if available. 
	# otherwise, download the data from its origin
	cd /vagrant/sauce_data

	# Grab sequence annotation files, relative to both strain-of-interest and reference, for each strain
	urlbase="http://mus.well.ox.ac.uk/19genomes/annotations/consolidated_annotation_9.4.2011/gene_models/"

	fetch_extract_gff3 $urlbase "consolidated_annotation.Col_0"
	fetch_extract_gff3 $urlbase "consolidated_annotation.Bur_0"
	fetch_extract_gff3 $urlbase "consolidated_annotation.Can_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ct_1"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Edi_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Hi_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Kn_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ler_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Mt_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.No_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Oy_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Po_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Rsch_4"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Sf_2"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Tsu_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wil_2"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ws_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wu_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Zu_0"

	# Annotations relative to the reference
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Col_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Bur_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Can_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ct_1.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Edi_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Hi_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Kn_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ler_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Mt_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.No_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Oy_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Po_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Rsch_4.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Sf_2.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Tsu_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wil_2.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Ws_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Wu_0.Col_0"
	# fetch_extract_gff3 $urlbase "consolidated_annotation.Zu_0.Col_0"

	# Fetch the chromosomal TAIR10 reference sequence. We'll get the RNAs out using the .gff3 annotation
	urlbase="ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/"
	fetch_raw $urlbase "TAIR10_chr1.fas"
	fetch_raw $urlbase "TAIR10_chr2.fas"
	fetch_raw $urlbase "TAIR10_chr3.fas"
	fetch_raw $urlbase "TAIR10_chr4.fas"
	fetch_raw $urlbase "TAIR10_chr5.fas"
	# fetch_raw $urlbase "TAIR10_chrC.fas" # DLing chloroplast and mitochondrial genomes is pointless.
	# fetch_raw $urlbase "TAIR10_chrM.fas"

	# Combine all the chromosome files together - this will make the parsing easier
	cat TAIR10_chr* > TAIR10_combined.fas

	# Grab chromosomal sequences for each of the 18 other strains
	urlbase="http://mus.well.ox.ac.uk/19genomes/fasta/MASKED/"

	fetch_raw $urlbase "bur_0.v7.PR_in_lowercase.fas"
	fetch_raw $urlbase "can_0.v7.PR_in_lowercase.fas"
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

# Prints aesthetically pleasing messages into the terminal.
function pretty_print() {
	printf "\n-=-=-=-=[ $1 ]=-=-=-=-\n"
    len=${#1} ch='-'
    padding=$(printf '%*s' "$len" | tr ' ' "$ch")
    printf "          $padding\n\n"
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

# this is only really useful in the production environment
# since we prefer Flask's own web server for development
function install_apache_wsgi() {
	pretty_print "Installing Apache WSGI"
	apt-get install -y apache2
	apt-get install -y libapache2-mod-wsgi
	cp /vagrant/bootstrap/000-default.conf /etc/apache2/sites-available/000-default.conf
	service apache2 restart
}