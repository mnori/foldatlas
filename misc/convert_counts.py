# Script to convert  2-line counts format into a nicer 1-line format

input_file = open("log_react_f_N1A_cdna_n_before.txt", "r")
output_file = open("log_react_f_N1A_cdna_n_before.formatted.txt", "w")

line = input_file.readline()
while(True):
    transcript_id = line.strip()
    line = input_file.readline().strip()

    if line[:2] == "AT": # another transcript id
        continue;

    if line == "": # if empty, it means end of file was reached
        break

    output_file.write(transcript_id+"\t"+line+"\n")

input_file.close()
output_file.close()
