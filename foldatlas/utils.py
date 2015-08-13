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
			seq_record.id = transcript_id
			SeqIO.write(seq_record, output_handle, "fasta")

			n += 1
			if n % 100 == 0:
				print("["+str(n)+"] sequences written")

		output_handle.close()
		print("...Finished exporting")
