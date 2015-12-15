live = False

if live:
	app_base_url = "http://www.foldatlas.com/" # if this is wrong, some ajax will fail
	static_base_url = "/static"
	static_path = "/var/www/foldatlas/static"
	dbuser = "root"
	dbpassword = "s7Alvwh801mcZ" # don't put the real live password here. change it on server instead.
else:
	app_base_url = "http://foldatlas.dev"
	static_base_url = "http://static.foldatlas.dev"
	static_path = "/vagrant/static"
	dbuser = "root"
	dbpassword = "vagrant"

# This defines hostname, database name, username and password for connecting to the DB.
database_uri = "mysql+mysqlconnector://"+dbuser+":"+dbpassword+"@127.0.0.1/foldatlas?charset=utf8&use_unicode=0"

# Points to the general data folder
data_folder = "/vagrant/sauce_data"

db_name="foldatlas"

# Points to structure data folder, which contains a *lot* of files
structure_data_folder = "/vagrant/structure_data"

dms_reactivities_experiment = {
	"nucleotide_measurement_run_id": 1,
	"strain_id": "Col_0",
	"nucleotides_filepath": data_folder+"/a_thaliana_compiled_counts.txt",
	"description": "DMS reactivities"
}

ribosome_profile_experiment = {
	"nucleotide_measurement_run_id": 2,
	"strain_id": "Col_0",
	"nucleotides_filepath": data_folder+"/p_site_counts_all.txt",
	# "coverage_filepath": data_folder+"/riboseq_coverage.txt", # just get coverage from summed counts
	"description": "Ribosome occupancies",
}

structures_in_silico = {
	"structure_prediction_run_id": 1,
	"strain_id": "Col_0",
	"description": "In silico structure prediction",
	"sauce_filepath": structure_data_folder+"/in_silico_structures",
	"sauce_ext": ".ct",
}

structures_in_vivo = {
	"structure_prediction_run_id": 2,
	"strain_id": "Col_0",
	"description": "In vivo experimental structure prediction",
	"sauce_filepath": structure_data_folder+"/in_vivo_structures",
	"sauce_ext": ".ct",
}

transcripts_fasta_filepath = data_folder+"/transcripts.fasta"

base_path = "/vagrant/foldatlas"

genoverse_base = static_base_url+"/genoverse"

# this is the one that will be displayed by the genome browser.
reference_strain_id = "Col_0"

# path of folder for temporary files
temp_folder = "/tmp"

# Strain metadata. This is used to parse from sauce files when hydrating the DB.
strains = [
	{
		"name": "Col_0",
		"description": "TAIR 10 Columbia reference ecotype",
		"sequence_filename": "TAIR10_combined.fas",
		"annotation_filename": "consolidated_annotation.Col_0.gff3"
	} # , 

	# {
	# 	"name": "Bur_0",
	# 	"description": "Bur_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "bur_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Bur_0.gff3"
	# }, {
	# 	"name": "Can_0",
	# 	"description": "Can_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "can_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Can_0.gff3"
	# }

	# comment these out for the real thing
	# , {
	# 	"name": "Ct_1",
	# 	"description": "Ct_1 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "ct_1.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Ct_1.gff3"
	# }, {
	# 	"name": "Edi_0",
	# 	"description": "Edi_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "edi_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Edi_0.gff3"
	# }, {
	# 	"name": "Hi_0",
	# 	"description": "Hi_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "hi_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Hi_0.gff3"
	# }, {
	# 	"name": "Kn_0",
	# 	"description": "Kn_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "kn_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Kn_0.gff3"
	# }, {
	# 	"name": "Ler_0",
	# 	"description": "Ler_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "ler_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Ler_0.gff3"
	# }, {
	# 	"name": "Mt_0",
	# 	"description": "Mt_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "mt_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Mt_0.gff3"
	# }, {
	# 	"name": "No_0",
	# 	"description": "No_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "no_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.No_0.gff3"
	# }, {
	# 	"name": "Oy_0",
	# 	"description": "Oy_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "oy_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Oy_0.gff3"
	# }, {
	# 	"name": "Po_0",
	# 	"description": "Po_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "po_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Po_0.gff3"
	# }, {
	# 	"name": "Rsch_4",
	# 	"description": "Rsch_4 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "rsch_4.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Rsch_4.gff3"
	# }, {
	# 	"name": "Sf_2",
	# 	"description": "Sf_2 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "sf_2.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Sf_2.gff3"
	# }, {
	# 	"name": "Tsu_0",
	# 	"description": "Tsu_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "tsu_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Tsu_0.gff3"
	# }, {
	# 	"name": "Wil_2",
	# 	"description": "Wil_2 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "wil_2.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Wil_2.gff3"
	# }, {
	# 	"name": "Ws_0",
	# 	"description": "Ws_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "ws_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Ws_0.gff3"
	# }, {
	# 	"name": "Wu_0",
	# 	"description": "Wu_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "wu_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Wu_0.gff3"
	# }, {
	# 	"name": "Zu_0",
	# 	"description": "Zu_0 strain, sequenced by the 19 Genomes project",
	# 	"sequence_filename": "zu_0.v7.PR_in_lowercase.fas",
	# 	"annotation_filename": "consolidated_annotation.Zu_0.gff3"
	# }
	
]

ignored_chromosomes = set(["chloroplast", "mitochondria"])
