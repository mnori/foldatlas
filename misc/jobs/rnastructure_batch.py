# Python script for doing batches of RNAstructure predictions.
# For running on HPC cluster.
# @author Matthew Norris <matthew.norris@jic.ac.uk>

import os
import sys

chunk_path = os.path.expanduser(sys.argv[1])
rnastructure_path = os.path.expanduser("~/RNAstructure")

def run_rnastructure(transcript_path):
	print(transcript_path)
	os.environ['DATAPATH'] = rnastructure_path+"/data_tables"
	os.system(
		rnastructure_path+"/exe/Fold"+ # the program
		" "+transcript_path+"/seq.fasta"+ # the sequence of interest
		" "+transcript_path+"/results.txt" # the output filename
	)

# Go through each transcript in turn, run rnastructure on each one.
for transcript_id in os.listdir(chunk_path):
	if os.path.isdir(chunk_path+"/"+transcript_id):
		run_rnastructure(chunk_path+"/"+transcript_id)
