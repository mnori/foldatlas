import time, os

# Useful for tracking down slow code.
# @author Matthew Norris
class Timeline():

	def __init__(self, name="Timeline"):
		self.name = name
		self.entries = []

	def log(self, name):
		self.entries.append(TimelineEntry(name))

	def dump(self):
		print("Timeline.dump() invoked")
		for i in range(0, len(self.entries) - 1):
			entry_a = self.entries[i]
			entry_b = self.entries[i + 1]
			t_diff = entry_b.time - entry_a.time;
			print("["+self.name+"] ["+entry_a.name+"] => ["+entry_b.name+"]: "+str(t_diff));

		t_tot = self.entries[len(self.entries) - 1].time - self.entries[0].time
		print("[TOTAL]: "+str(t_tot)+"")

# Helper class for Timeline
class TimelineEntry():

	def __init__(self, name):
		self.name = name
		self.time = time.time()

def ensure_dir(f):
	if not os.path.exists(f):
		os.makedirs(f)

# Create a big fasta containing all transcript sequences stored in the DB.
class FastaExporter():
	def export(self):
		print("Exporting...")

		from Bio import SeqIO
		from database import db_session
		from models import Transcript
		import settings 

		# Get all transcript IDs
		results = db_session \
			.query(Transcript.id) \
			.all()

		n = 0

		target_filepath = settings.data_folder+"/transcripts.fasta"

		output_handle = open(target_filepath, "w")
		for result in results:
			transcript_id = result[0]
			seq_record = Transcript(transcript_id).get_sequence()
			
			if seq_record == None:
				print ("Missing sequence for ["+transcript_id+"]")
				continue

			seq_record.id = transcript_id
			SeqIO.write(seq_record, output_handle, "fasta")

			n += 1
			if n % 100 == 0:
				print("["+str(n)+"] sequences written")

		output_handle.close()

		print("...Saved transcripts to ["+target_filepath+"]")

# Split big fasta file into smaller ones. This is for doing HPC structure predictions
class FastaSplitter():

	def __init__(self):
		import settings

		self.n_chunks = 256
		self.sauce_filepath = settings.data_folder+"/transcripts.fasta"
		self.target_dirpath = settings.data_folder+"/rnastructure_seqs"

	def split(self):
		from Bio import SeqIO
		from utils import ensure_dir

		# Count the sequences

		print("Counting sequences...")
		n_seqs = 0
		handle = open(self.sauce_filepath, "r")
		for record in SeqIO.parse(handle, "fasta"):
			n_seqs += 1
		handle.close()
		chunk_size = int(n_seqs / self.n_chunks)
		print("...There are ["+str(n_seqs)+"] sequences")

		# Iterate through chunks of sequences.
		# For each sequence, write a single fasta file in the right location.
		n_seqs = 0
		handle = open(self.sauce_filepath, "r")
		for record in SeqIO.parse(handle, "fasta"):
			chunk_n = int(n_seqs / chunk_size)
			chunk_dir = self.target_dirpath+"/chunk_"+str(chunk_n)+"/"+record.id
			ensure_dir(chunk_dir+"/")

			f = open(chunk_dir+"/seq.fasta", "w")
			SeqIO.write(record, f, "fasta")
			f.close()

			n_seqs += 1

			if n_seqs % 100 == 0:
				print("["+str(n_seqs)+"] fasta files written")

		handle.close()		

# Grab structures from storage so they can be easily imported
def grab_structures():
	import os, shutil
	from database import db_session
	from models import Gene

	print("grab_structures() invoked");

	# this should be a symlink to /media/shares/Research-Groups/Yiliang-Ding/data_analysis_Ding_2013/MAC/Yin/Mapping_F/raw_data/
	sauce_folder = "/vagrant/raw_structures_folder"
	dest_folder = os.path.expanduser("~/foldatlas/sauce_data/structures")

	# grab ALL genes in the DB.

	genes_to_grab = []

	data = db_session.query(Gene).all()
	for gene in data:
		genes_to_grab.append(gene.id)

	print(genes_to_grab)

def insert_newlines(string, every=80):
	lines = []
	for i in range(0, len(string), every):
		lines.append(string[i:i+every])
	return '\n'.join(lines)
	
def build_dot_bracket(positions):
	# build dot bracket string
	n_reverse = n_forward = 0
	dot_bracket_str = ""

	for curr_position in range(1, len(positions) + 1):
		paired_to_position = positions[curr_position - 1]
		if paired_to_position == 0:
			dot_bracket_str += "."
		elif paired_to_position < curr_position:
			n_reverse += 1
			dot_bracket_str += ")"
		elif paired_to_position > curr_position:
			n_forward += 1
			dot_bracket_str += "("
		else:
			# should never happen
			print("Error: cannot do self pairing!");
			dot_bracket_str += "."

	if n_reverse != n_forward:
		return "ERROR: n_reverse != n_forward!"
	return dot_bracket_str


	# genes_to_grab = {"AT1G01010", "AT1G01020", "AT1G01030", "AT1G01040", "AT1G01046", "AT1G01050", "AT1G01060", "AT1G01070", "AT1G01073", "AT1G01080", "AT1G01090", "AT1G01100", "AT1G01110", "AT1G01115", "AT1G01120", "AT1G01130", "AT1G01140", "AT1G01150", "AT1G01160", "AT1G01170", "AT1G01180", "AT1G01183", "AT1G01190", "AT1G01200", "AT1G01210", "AT1G01220", "AT1G01225", "AT1G01230", "AT1G01240", "AT1G01250", "AT1G01260", "AT1G01270", "AT1G01280", "AT1G01290", "AT1G01300", "AT1G01305", "AT1G01310", "AT1G01320", "AT1G01340", "AT1G01350", "AT1G01355", "AT1G01360", "AT1G01370", "AT1G01380", "AT1G01390", "AT1G01400", "AT1G01410", "AT1G01420", "AT1G01430", "AT1G01440", "AT1G01448", "AT2G01010", "AT2G01020", "AT2G01021", "AT2G01023", "AT2G01050", "AT2G01060", "AT2G01070", "AT2G01080", "AT2G01090", "AT2G01100", "AT2G01110", "AT2G01120", "AT2G01130", "AT2G01140", "AT2G01150", "AT2G01160", "AT2G01170", "AT2G01175", "AT2G01180", "AT2G01190", "AT2G01200", "AT2G01210", "AT2G01220", "AT2G01240", "AT2G01250", "AT2G01260", "AT2G01270", "AT2G01275", "AT2G01280", "AT2G01290", "AT2G01300", "AT2G01310", "AT2G01320", "AT2G01330", "AT2G01340", "AT2G01350", "AT2G01360", "AT2G01370", "AT2G01390", "AT2G01400", "AT2G01410", "AT2G01420", "AT2G01422", "AT2G01430", "AT2G01440", "AT2G01450", "AT2G01460", "AT2G01470", "AT2G01480", "AT2G01490", "AT3G01015", "AT3G01020", "AT3G01030", "AT3G01040", "AT3G01050", "AT3G01060", "AT3G01070", "AT3G01080", "AT3G01085", "AT3G01090", "AT3G01100", "AT3G01120", "AT3G01130", "AT3G01140", "AT3G01142", "AT3G01150", "AT3G01160", "AT3G01170", "AT3G01175", "AT3G01180", "AT3G01185", "AT3G01190", "AT3G01200", "AT3G01202", "AT3G01210", "AT3G01220", "AT3G01230", "AT3G01240", "AT3G01250", "AT3G01260", "AT3G01270", "AT3G01280", "AT3G01290", "AT3G01300", "AT3G01310", "AT3G01311", "AT3G01313", "AT3G01316", "AT3G01319", "AT3G01320", "AT3G01322", "AT3G01323", "AT3G01324", "AT3G01325", "AT3G01326", "AT3G01327", "AT3G01328", "AT3G01329", "AT3G01330", "AT3G01331", "AT4G00020", "AT4G00026", "AT4G00030", "AT4G00040", "AT4G00050", "AT4G00060", "AT4G00070", "AT4G00080", "AT4G00085", "AT4G00090", "AT4G00100", "AT4G00110", "AT4G00120", "AT4G00124", "AT4G00130", "AT4G00140", "AT4G00150", "AT4G00160", "AT4G00165", "AT4G00170", "AT4G00180", "AT4G00190", "AT4G00200", "AT4G00210", "AT4G00220", "AT4G00230", "AT4G00231", "AT4G00232", "AT4G00234", "AT4G00238", "AT4G00240", "AT4G00250", "AT4G00260", "AT4G00270", "AT4G00280", "AT4G00290", "AT4G00300", "AT4G00305", "AT4G00310", "AT4G00315", "AT4G00320", "AT4G00330", "AT4G00335", "AT4G00340", "AT4G00342", "AT4G00350", "AT4G00355", "AT4G00360", "AT4G00370", "AT4G00380", "AT5G01015", "AT5G01020", "AT5G01030", "AT5G01040", "AT5G01050", "AT5G01060", "AT5G01070", "AT5G01075", "AT5G01080", "AT5G01090", "AT5G01100", "AT5G01110", "AT5G01120", "AT5G01130", "AT5G01140", "AT5G01150", "AT5G01160", "AT5G01170", "AT5G01175", "AT5G01180", "AT5G01190", "AT5G01200", "AT5G01210", "AT5G01215", "AT5G01220", "AT5G01225", "AT5G01230", "AT5G01240", "AT5G01250", "AT5G01260", "AT5G01270", "AT5G01280", "AT5G01290", "AT5G01300", "AT5G01310", "AT5G01320", "AT5G01330", "AT5G01340", "AT5G01350", "AT5G01360", "AT5G01365", "AT5G01370", "AT5G01380", "AT5G01390", "AT5G01400", "AT5G01410", "AT5G01420", "AT5G01430", "AT5G01440", "AT5G01445"}

	# def process_folder(sauce_folder, dest_folder):
	# 	files = os.listdir(sauce_folder)
	# 	for filename in files:
	# 		for gene in genes_to_grab:
	# 			if gene in filename:
	# 				shutil.copyfile(sauce_folder+"/"+filename, dest_folder+"/"+filename)

	# process_folder(sauce_folder+"/in_silico_structures", dest_folder+"/in_silico")
	# process_folder(sauce_folder+"/in_vivo_structures", dest_folder+"/in_vivo")
