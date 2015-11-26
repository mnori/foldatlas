import os
import shutil
from os.path import expanduser

# Renames .ct files to their transcript IDs.

def rename_folder(folder_filepath):
	folder_filepath = expanduser(folder_filepath)
	filenames = os.listdir(folder_filepath)
	for filename in filenames:
		if filename[-3:] != ".ct":
			print("Skipped "+filename)
			continue

		sauce_filepath = folder_filepath+"/"+filename
		with open(sauce_filepath, "r") as f:
			line = f.readline()
			transcript_id = line.strip().split()[-1]
			dest_filepath = folder_filepath+"/"+transcript_id+".ct"
			shutil.move(sauce_filepath, dest_filepath)
			print("Moved to ["+dest_filepath+"]")

rename_folder("~/data_input/foldatlas_structures/in_vivo_structures")
rename_folder("~/data_input/foldatlas_structures/in_silico_structures")


