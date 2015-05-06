
# This defines hostname, database name, username and password for connecting to the DB.
database_uri = "mysql+mysqlconnector://root:vagrant@127.0.0.1/rnabrowser?charset=utf8&use_unicode=0"

genomes_sauce_folder = "/vagrant/sauce_data"

base_path = "/vagrant/rnabrowser"

static_base = "http://static.rnabrowser.dev"
genoverse_base = static_base+"/genoverse"

# this is the one that will be displayed by the genome browser.
reference_strain_id = "Col_0"

# Strain metadata. This is used to parse from sauce files when hydrating the DB.
strains = [
	{
		"name": "Col_0",
		"description": "TAIR 10 Columbia reference ecotype",
		"sequence_filename": "TAIR10_combined.fas",
		"annotation_filename": "consolidated_annotation.Col_0.gff3"
	}
]
