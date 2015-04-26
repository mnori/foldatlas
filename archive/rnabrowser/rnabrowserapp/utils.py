import os
import errno
from Bio import SeqIO
from urllib.request import urlopen

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Save sequences to a Fasta file
def save_seqs_fasta(sequences, path):
    print("Writing sequences...")
    output_handle = open(path, "w")
    count = SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()
    print("Saved " + str(count) + " sequences to " + path)

# Download whole genomes
def dl_fastas(data_dir, prefix, sauce_url, sauce_filenames) :
    for sauce_filename in sauce_filenames:
        target_file = data_dir + prefix + "_" + sauce_filename

        # skip the file if it's alreadt downloaded
        if os.path.isfile(target_file):
            print("Skipped " + sauce_filename+" - already exists")
            continue

        # open the remote file
        remote_file = urlopen(sauce_url+sauce_filename)

        # open the local file
        output = open(target_file, 'wb')

        # write data to local file & close it
        output.write(remote_file.read())
        output.close()
        print("Downloaded " + sauce_filename)

def load_sequences_fasta(filepath):
    input_handle = open(filepath, "rU")
    seqs = list(SeqIO.parse(input_handle, "fasta"))
    input_handle.close()
    return seqs