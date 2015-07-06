from database import db_session;
from sqlalchemy import and_

import json, database, settings, uuid, os, subprocess
from models import Feature, Transcript, AlignmentEntry, NucleotideMeasurement, Experiment, \
    TranscriptCoverage, StructurePosition

from utils import ensure_dir

# Fetches sequence annotation data from the DB and sends it to the genome
# browser front end as JSON.

class GenomeBrowser():

    def get_transcripts(self, request):

        chromosome_id = "Chr"+str(int(request.args.get('chr'))) # SQL-injection safe
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        # Retrieve features using the gene location cache table
        sql = ( "SELECT feature.* "
                "FROM gene_location, transcript, feature "
                "WHERE gene_location.strain_id = '"+settings.reference_strain_id+"' "
                "AND gene_location.chromosome_id = '"+chromosome_id+"' "
                "AND gene_location.end > '"+str(start)+"' "
                "AND gene_location.start < '"+str(end)+"' "
                "AND gene_location.gene_id = transcript.gene_id "
                "AND transcript.id = feature.transcript_id "
                "AND feature.strain_id = '"+settings.reference_strain_id+"'")

        results = database.engine.execute(sql)

        # collect transcript data
        transcripts = {}
        feature_rows = []
        for result in results:
            if result.transcript_id not in transcripts:
                transcripts[result.transcript_id] = {
                    "Parent": result.transcript_id,
                    "feature_type": "transcript", # without this, it won't draw
                    "direction": result.direction,
                    "start": None,
                    "end": None,
                    "id": result.transcript_id
                }

            transcript = transcripts[result.transcript_id]

            # keep track of total start and end
            if transcript["start"] == None or result.start < transcript["start"]:
                transcript["start"] = result.start
            if transcript["end"] == None or result.end > transcript["end"]:
                transcript["end"] = result.end

            feature_rows.append(result)

        out = []

        # add the transcript metadata to the output. make sure the transcripts are added
        # in alphabetical order
        transcript_ids = []
        for transcript_id in transcripts:
            transcript_ids.append(transcript_id)

        transcript_ids = sorted(transcript_ids)
        for transcript_id in transcript_ids:
            out.append(transcripts[transcript_id])

        # also add all the feature metadata to the output
        for feature_row in feature_rows:
            out.append({
                "Parent": feature_row.transcript_id,
                "feature_type": feature_row.type_id,
                "direction": result.direction,
                "start": feature_row.start,
                "end": feature_row.end,
                "id": feature_row.transcript_id+"-"+str(feature_row.id)
            })

        return json.dumps(out)

    def get_genes(self, request):
        
        from utils import Timeline

        chromosome_id = "Chr"+str(int(request.args.get('chr'))) # SQL-injection safe
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        # fetch gene data from the location cache table.
        sql = ( "SELECT * FROM gene_location "
                "WHERE strain_id = '"+settings.reference_strain_id+"' "
                "AND chromosome_id = '"+chromosome_id+"' "
                "AND end > '"+str(start)+"' "
                "AND start < '"+str(end)+"'")

        results = database.engine.execute(sql)

        out = []
        for result in results:
            out.append({ 
                "feature_type": "gene", # without this, it won't draw
                "direction": result.direction,
                "id": result.gene_id,
                "start": result.start,
                "end": result.end,
            })        

        buf = json.dumps(out)
        return buf

    # Fetch chromosome IDs and their lengths. Used for chromosome menu and also initialising the genome browser.
    def get_chromosomes(self):

        sql = ( "SELECT chromosome_id, CHAR_LENGTH(sequence) length FROM chromosome "
                "WHERE strain_id = '"+settings.reference_strain_id+"' "
                "ORDER BY chromosome_id ASC")

        results = database.engine.execute(sql)

        out = []
        for result in results:
            out.append({
                "chromosome_id": result.chromosome_id,
                "length": result.length,
                "int_id": int(result.chromosome_id[3])
            })

        return out

class TranscriptView():

    def __init__(self, transcript_id):

        # TODO try to fetch transcript ID first, then if it works do the rest of the stuff
        self.transcript_id = transcript_id
        self.structure_view = StructureView(self.transcript_id, settings.reference_strain_id)
        self.nucleotide_measurement_view = NucleotideMeasurementView(self.transcript_id, settings.reference_strain_id)
        self.alignment_view = AlignmentView(self.transcript_id)

class NucleotideMeasurementView():

    def __init__(self, transcript_id, strain_id):
        self.transcript_id = transcript_id
        self.strain_id = strain_id
        self.build_entries([1, 2])

    def build_entries(self, experiment_ids):

        # Load experiments
        experiments = db_session \
            .query(Experiment) \
            .filter(Experiment.id.in_(experiment_ids)) \
            .all()

        # Load measurements
        seq_str = str(Transcript(self.transcript_id).get_sequence(self.strain_id).seq)
        measurements = db_session \
            .query(NucleotideMeasurement) \
            .filter(
                NucleotideMeasurement.experiment_id.in_(experiment_ids),
                NucleotideMeasurement.strain_id==self.strain_id,
                NucleotideMeasurement.transcript_id==self.transcript_id
            ) \
            .all()

        data = {}
        # populate experiment rows
        for experiment in experiments:
            experiment_data = {
                "id": experiment.id,
                "type": experiment.type,
                "description": experiment.description,
                "data": []
            }
            
            for n in range(len(seq_str)): # initialise the array
                experiment_data["data"].append({
                    "position": n,
                    "nuc": seq_str[n],
                    "measurement": None
                })
            data[experiment.id] = experiment_data

        # add measurements to each experiment
        for measurement_row in measurements: # add values where present
            experiment_id = measurement_row.experiment_id
            pos = measurement_row.position - 1
            data[experiment_id]["data"][pos]["measurement"] = measurement_row.measurement

        self.data_json = json.dumps(data)
        
class AlignmentView():

    alignment_line_length = 80

    def __init__(self, transcript_id):
        self.transcript_id = transcript_id
        self.build_alignment_entries()

    def build_alignment_entries(self):
        self.alignment_rows = []

        # fetch the alignment rows from the DB, using the ORM
        alignment_entries = db_session \
            .query(AlignmentEntry) \
            .filter(AlignmentEntry.transcript_id==self.transcript_id) \
            .all()

        if (len(alignment_entries) == 0):
            return # not enough transcripts to align

        aln_len = len(alignment_entries[0].sequence) # length of alignment, including gaps
        row_n = 0
        reached_end = False
        seq_len_processed = 0

        # initialise tot_nucs counters. these are for showing nuc counts at the ends of each alignment row.
        nuc_counts = {}
        for alignment_entry in alignment_entries:
            nuc_counts[alignment_entry.strain_id] = 0

        while(True): # Each iteration builds 1 row of alignment data

            start = row_n * self.alignment_line_length
            end = start + self.alignment_line_length

            if aln_len < end:
                reached_end = True
                end = aln_len

            self.alignment_rows.append({
                "strain_data": {},
                "diff": list("*" * (end - start))
            })

            # create diff - as "*" - then change to "." when a difference is encountered
            # create alignment entries data structure, for showing the sequences        
            for alignment_entry in alignment_entries:
                self.alignment_rows[row_n]["strain_data"][alignment_entry.strain_id] = {
                    "nuc_count": 0, # TODO fill this shiz out
                    "sequence": list(alignment_entry.sequence[start : end])
                }

            # Loop through each nucleotide in the sequence. Determine any differences between the 
            # strains at the position of interest. Store in "diff" variable
            for n in range(start, end):
                different = False
                old_nuc = None
                for alignment_entry in alignment_entries:
                    new_nuc = alignment_entry.sequence[n]

                    if new_nuc != "-": # keep track of nucleotide counts, for showing on the end
                        nuc_counts[alignment_entry.strain_id] += 1

                    if old_nuc != None and new_nuc != old_nuc:
                        self.alignment_rows[row_n]["diff"][n - start] = "."
                    old_nuc = new_nuc

            # add nucleotide counts to the ends of the sequence alignment.
            for alignment_entry in alignment_entries:
                self.alignment_rows[row_n]["strain_data"][alignment_entry.strain_id]["nuc_count"] = nuc_counts[alignment_entry.strain_id]

            if reached_end:
                break

            row_n += 1

class TranscriptSearcher():
    def search(self, search_string):
        transcripts = db_session \
            .query(Transcript) \
            .filter(Transcript.id.like("%"+search_string+"%")) \
            .all()

        out = []
        for transcript in transcripts:
            out.append(transcript.id)

        return json.dumps(out)

class CoverageSearcher():
    def __init__(self):

        # size of pages
        self.page_size = 25

        # The experiment ID to sort by. Ideally this should have a value for each
        # transcript, otherwise there will be some missing transcripts...
        self.experiment_id = 2

    def fetch_page_count(self):
        # better to do the imports closer to where they are needed
        from sqlalchemy import func
        from math import ceil

        transcript_count = db_session \
            .query(func.count('*')) \
            .select_from(TranscriptCoverage) \
            .filter(TranscriptCoverage.experiment_id==self.experiment_id) \
            .scalar()

        page_count = ceil(transcript_count / self.page_size)
        return page_count

    def fetch_transcript_data(self, page_num):

        from models import Structure

        results = db_session \
            .query(
                TranscriptCoverage
            ) \
            .filter(
                TranscriptCoverage.experiment_id==self.experiment_id
            ) \
            .outerjoin((
                Structure, 
                Structure.transcript_id==TranscriptCoverage.transcript_id,
            )) \
            .add_entity(Structure) \
            .group_by(TranscriptCoverage.transcript_id) \
            .order_by(TranscriptCoverage.measurement.desc()) \
            .offset((int(page_num) - 1) * self.page_size) \
            .limit(str(self.page_size)) \
            .all()

        out = []
        for result in results:
            out.append({
                "coverage": result[0],
                "structure": result[1]
            })

        return out

class StructureView():
    def __init__(self, transcript_id, strain_id):
        self.transcript_id = transcript_id
        self.strain_id = strain_id
        self.build_entries([3, 4])

    def build_entries(self, experiment_ids):

        from models import Structure

        # Load experiments
        experiments = db_session \
            .query(Experiment) \
            .filter(Experiment.id.in_(experiment_ids)) \
            .all()

        data = {}

        for experiment in experiments:

            experiment_data = {
                "id": experiment.id,
                "type": experiment.type,
                "description": experiment.description,
                "data": []
            }

            # fetch all Structure objects that match the IDs
            results = db_session \
                .query(Structure) \
                .filter(
                    Structure.experiment_id==experiment.id,
                    Structure.strain_id==self.strain_id,
                    Structure.transcript_id==self.transcript_id
                ) \
                .all()

            # add the structures to output json
            for structure in results:
                experiment_data["data"].append({
                    "id": structure.id,
                    "energy": structure.energy,
                    "pc1": structure.pc1,
                    "pc2": structure.pc2
                })

            data[experiment.id] = experiment_data

        self.data_json = json.dumps(data)

# Plots an RNA structure using the RNAplot program from the ViennaRNA package.
class StructurePlotView():
    def __init__(self, structure_id):
        self.structure_id = structure_id
        self.build_plot()

    def build_plot(self):
        unique_folder = "/tmp/"+str(uuid.uuid4())
        ensure_dir(unique_folder)

        # generate dot bracket file for RNAplot
        dot_bracket_filepath = unique_folder+"/dotbracket.txt"
        dot_bracket_data = self.build_dot_bracket()
        
        self.dot_bracket_json = json.dumps(dot_bracket_data)

    # Deprecated method to make Vienna RNAplot image. These plots are a bit shit
    # compared to those of forna.
    def get_vienna_svg(self, dot_bracket_str):
        dot_bracket_file = open(dot_bracket_filepath, "w")
        dot_bracket_file.write(dot_bracket_str+"\n")
        dot_bracket_file.close()

        # change to tmp folder
        os.chdir(unique_folder)

        # use RNAplot to generate the postscript image
        os.system("RNAplot < "+dot_bracket_filepath)

        # convert the plot to SVG, capturing the SVG string
        result = str(subprocess.check_output("pstoedit -f plot-svg rna.ps", shell=True))

        # Chop out some unwanted text
        result = result[result.find("<svg") : ]
        result = result[0 : result.find("</svg>")]

        # return the SVG so we can show it on the front end
        return result

    def build_dot_bracket(self):

        # get all the positions
        results = db_session \
            .query(StructurePosition) \
            .filter(
                StructurePosition.structure_id==self.structure_id
            ) \
            .order_by(StructurePosition.position) \
            .all()

        # build dot bracket string
        n_reverse = n_forward = 0
        dot_bracket_str = seq_str = ""
        for result in results:
            seq_str += result.letter
            if result.paired_to_position == 0:
                dot_bracket_str += "."
            elif result.paired_to_position < result.position:
                n_reverse += 1
                dot_bracket_str += ")"
            elif result.paired_to_position > result.position:
                n_forward += 1
                dot_bracket_str += "("
            else:
                # should never happen
                print("Error: cannot do self pairing!");
                dot_bracket_str += "."

        # # these must be indentical!
        # print("n_reverse: "+str(n_reverse))
        # print("n_forward: "+str(n_forward))

        return {
            "seq": seq_str,
            "dot_bracket": dot_bracket_str
        }

            
