import time

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

class TimelineEntry():

	def __init__(self, name):
		self.name = name
		self.time = time.time()

time.time()
