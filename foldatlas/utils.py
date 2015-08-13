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
		output_handle = open(settings.genomes_sauce_folder+"/transcripts.fasta", "w")
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
		print("...Finished exporting")

# Split big fasta file into smaller ones. This is for doing HPC structure predictions
class FastaSplitter():

	def __init__(self):
		import settings

		self.n_chunks = 256
		self.sauce_filepath = settings.genomes_sauce_folder+"/transcripts.fasta"
		self.target_dirpath = settings.genomes_sauce_folder+"/rnastructure_seqs"

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



