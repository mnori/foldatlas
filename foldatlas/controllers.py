from database import db_session;
from sqlalchemy import and_

# # try to use python2 module in python3
# from past import autotranslate
# autotranslate(['forgi'])

# this import is slow. almost certainly due to autotranslate.
# might have to build a pre-translated version of this module.
# import forgi.graph.bulge_graph as fgb

import json, database, settings, uuid, os, subprocess
from models import Feature, Transcript, NucleotideMeasurement, \
    NucleotideMeasurementSet, Structure, StructurePosition, GeneLocation, \
    NucleotideExperiment, StructurePredictionRun

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

        self.transcript_id = transcript_id

        # Get the coords of the associated gene
        data = db_session \
            .query(Transcript, GeneLocation) \
            .filter(
                Transcript.id==transcript_id,
                Transcript.gene_id==GeneLocation.gene_id,
                GeneLocation.strain_id==settings.reference_strain_id
            ) \
            .all()

        self.gene_id = data[0][1].gene_id

        self.transcript_data = json.dumps({
            "gene_id": self.gene_id,
            "transcript_id": transcript_id,
            "chromosome_id": data[0][1].chromosome_id,
            "start": data[0][1].start, 
            "end": data[0][1].end
        })

        self.structure_view = StructureView(self.transcript_id, settings.reference_strain_id)
        self.nucleotide_measurement_view = NucleotideMeasurementView(self.transcript_id, settings.reference_strain_id)
        
        self.empty = self.structure_view.empty and self.nucleotide_measurement_view.empty

        # disable alignment view... revisit later with SNPstructure
        # self.alignment_view = AlignmentView(self.transcript_id)

class NucleotideMeasurementView():

    def __init__(self, transcript_id, strain_id):
        self.transcript_id = transcript_id
        self.strain_id = strain_id
        self.build_entries([1, 2])

    def build_entries(self, experiment_ids):

        from models import NucleotideExperiment

        # Load experiments
        experiments = db_session \
            .query(NucleotideExperiment) \
            .filter(NucleotideExperiment.id.in_(experiment_ids)) \
            .all()

        # Load measurements
        seq_str = str(Transcript(self.transcript_id).get_sequence(self.strain_id).seq)
        measurements_data = db_session \
            .query(NucleotideMeasurementSet, NucleotideMeasurement) \
            .filter(
                NucleotideMeasurementSet.nucleotide_experiment_id.in_(experiment_ids),
                NucleotideMeasurementSet.transcript_id==self.transcript_id,
                NucleotideMeasurementSet.id==NucleotideMeasurement.nucleotide_measurement_set_id
            ) \
            .all()

        data = {}

        # Populate experiment rows
        for experiment in experiments:
            experiment_data = {
                "id": experiment.id,
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

        # Add measurements to each experiment
        for row in measurements_data: # add values where present

            measurement_set = row[0]
            measurement = row[1]

            experiment_id = measurement_set.nucleotide_experiment_id
            pos = measurement.position - 1
            data[experiment_id]["data"][pos]["measurement"] = measurement.measurement

        # For each experiment, check whether there is no data and set empty flags accordingly.
        self.empty = True # all empty flag
        for experiment_id in data:
            entry = data[experiment_id]

            empty = True
            for pos in entry["data"]:
                if pos["measurement"] != 0 and pos["measurement"] != None:
                    empty = False
                    self.empty = False

            if empty:
                del entry["data"]
                entry["empty"] = True
            else:
                entry["empty"] = False

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
        from flask import abort

        transcripts = db_session \
            .query(Transcript) \
            .filter(Transcript.id.like("%"+search_string+"%")) \
            .all()

        if len(transcripts) == 0: # no transcripts found
            abort(404)

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
        self.nucleotide_experiment_id = 1

    def fetch_page_count(self):
        # better to do the imports closer to where they are needed
        from sqlalchemy import func
        from math import ceil

        transcript_count = db_session \
            .query(func.count('*')) \
            .select_from(NucleotideMeasurementSet) \
            .filter(NucleotideMeasurementSet.nucleotide_experiment_id==self.nucleotide_experiment_id) \
            .scalar()

        page_count = ceil(transcript_count / self.page_size)
        return page_count

    def fetch_transcript_data(self, page_num):

        from sqlalchemy import func, and_
        from models import Structure, GeneLocation

        # optional in vivo query - will include some non-in vivo results
        # results = db_session \
        #     .query(NucleotideMeasurementSet, Transcript, GeneLocation,) \
        #     .filter(
        #         NucleotideMeasurementSet.nucleotide_experiment_id==self.nucleotide_experiment_id,
        #         Transcript.id==NucleotideMeasurementSet.transcript_id,
        #         Transcript.gene_id==GeneLocation.gene_id,
        #         GeneLocation.strain_id==settings.reference_strain_id # get this for gene len
        #     ) \
        #     .outerjoin(( # Left join to find in-vivo structures for structure indicator
        #         Structure, 
        #         and_(
        #             Structure.transcript_id==NucleotideMeasurementSet.transcript_id,

        #             # this filters so it's only in vivo considered
        #             Structure.structure_prediction_run_id==2 
        #         ) 
        #     )) \
        #     .add_entity(Structure) \
        #     .group_by(NucleotideMeasurementSet.transcript_id) \
        #     .order_by(NucleotideMeasurementSet.coverage.desc()) \
        #     .offset((int(page_num) - 1) * self.page_size) \
        #     .limit(str(self.page_size)) \
        #     .all()

        # mandatory in vivo query
        results = db_session \
            .query(NucleotideMeasurementSet, Transcript, GeneLocation, Structure, ) \
            .filter(
                NucleotideMeasurementSet.nucleotide_experiment_id==self.nucleotide_experiment_id,
                Transcript.id==NucleotideMeasurementSet.transcript_id,
                Transcript.gene_id==GeneLocation.gene_id,
                GeneLocation.strain_id==settings.reference_strain_id, # get this for gene len
                Structure.transcript_id==NucleotideMeasurementSet.transcript_id,

                # this filters so it's only in vivo considered
                Structure.structure_prediction_run_id==2    
            ) \
            .add_entity(Structure) \
            .group_by(NucleotideMeasurementSet.transcript_id) \
            .order_by(NucleotideMeasurementSet.coverage.desc()) \
            .offset((int(page_num) - 1) * self.page_size) \
            .limit(str(self.page_size)) \
            .all()

        out = []
        for result in results:
            out.append({
                "measurement_set": result[0],
                "structure": result[3],
                "gene_length": (result[2].end - result[2].start) + 1
            })

        return out

class StructureView():
    def __init__(self, transcript_id, strain_id):
        self.transcript_id = transcript_id
        self.strain_id = strain_id
        self.build_entries([1, 2])

    def build_entries(self, structure_prediction_run_ids):

        from models import Structure, StructurePredictionRun

        # Load experiments
        runs = db_session \
            .query(StructurePredictionRun) \
            .filter(StructurePredictionRun.id.in_(structure_prediction_run_ids)) \
            .all()

        data = {}

        for run in runs:

            run_data = {
                "id": run.id,
                "description": run.description,
                "data": []
            }

            # fetch all Structure objects that match the experiment ID and the transcript ID
            results = db_session \
                .query(Structure) \
                .filter(
                    Structure.structure_prediction_run_id==run.id,
                    Structure.transcript_id==self.transcript_id
                ) \
                .all()

            # add the structures to output json
            for structure in results:
                run_data["data"].append({
                    "id": structure.id,
                    "energy": structure.energy,
                    "pc1": structure.pc1,
                    "pc2": structure.pc2
                })

            data[run.id] = run_data

        self.empty = True
        for experiment_id in data:
            entry = data[experiment_id]
            if len(entry["data"]) > 0:
                self.empty = False

        if not self.empty:
            self.data_json = json.dumps(data)

# Plots an RNA structure using the RNAplot program from the ViennaRNA package.
class StructureDiagramView():
    def __init__(self, structure_id):
        self.structure_id = structure_id
        self.build_plot()

    def build_plot(self):
        # convert entities to dot bracket string
        data = self.build_dot_bracket()

        # use ViennaRNA to get 2d plot coords    
        data["coords"] = self.get_vienna_layout(data)

        # return the results as a json string
        self.data_json = json.dumps(data)

    def build_dot_bracket(self):

        # get all the positions
        results = db_session \
            .query(StructurePosition) \
            .filter(StructurePosition.structure_id==self.structure_id) \
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

        return {
            "sequence": seq_str.replace("T", "U"),
            "structure": dot_bracket_str
        }

        # RNA.cvar.rna_plot_type = 1
        # coords = RNA.get_xy_coordinates(bp_string)
        # xs = np.array([coords.get(i).X for i in range(len(bp_string))])
        # ys = np.array([coords.get(i).Y for i in range(len(bp_string))])

        # return zip(xs,ys)

    # Grab 2d coords from viennaRNA
    # There is a python2 wrapper for vienna RNA but not python 3 compatible
    def get_vienna_layout(self, data):

        temp_folder = "/tmp/"+str(uuid.uuid4())
        ensure_dir(temp_folder)
        dot_bracket_filepath = temp_folder+"/dotbracket.txt"

        f = open(dot_bracket_filepath, "w")
        f.write(data["sequence"]+"\n"+data["structure"]+"\n")
        f.close()

        # change to tmp folder
        os.chdir(temp_folder)

        # use RNAplot CLI to generate the xrna tab delimited file
        os.system("RNAplot -o xrna < "+dot_bracket_filepath)

        # get the coords out by parsing the file
        coords = []
        with open(temp_folder+"/rna.ss") as f:
            for line in f:
                line = line.strip()
                if line == "" or line[0] == "#":
                    continue

                bits = line.split()
                x = float(bits[2])
                y = float(bits[3])
                coords.append([x, y])

        os.system("rm -rf "+temp_folder)

        return coords
        
        # return result

class StructureCirclePlotView():
    def __init__(self, structure_id):
        self.structure_id = structure_id
        self.get_values()

    def get_values(self):
        # convert entities to dot bracket string
        out = [];

        # get all the positions
        results = db_session \
            .query(StructurePosition) \
            .filter(StructurePosition.structure_id==self.structure_id) \
            .order_by(StructurePosition.position) \
            .all()

        # build the output. backward facing links are left blank
        # results must be shifted back to array indexes, since they start at 1 in the DB.
        for result in results:
            if      result.paired_to_position == None or \
                    result.paired_to_position < result.position:

                link = None
            else:
                link = result.paired_to_position - 1

            out.append({
                "name": str(result.position - 1),
                "label": result.letter,
                "link": link
            })

        self.data_json = json.dumps(out)

# Generates plaintext structure text files for download
class StructureDownloader():
    def __init__(self, structure_prediction_run_ids, transcript_id):
        self.structure_prediction_run_ids = structure_prediction_run_ids
        self.transcript_id = transcript_id

    def generateTxt(self):

        # NOTE - more than 20 structures is likely to crash this shizzle
        # probably need to look at a chunk based system

        results = db_session \
            .query(Structure, StructurePosition, StructurePredictionRun) \
            .filter(
                Structure.structure_prediction_run_id.in_(self.structure_prediction_run_ids),
                StructurePredictionRun.id==Structure.structure_prediction_run_id,
                Structure.transcript_id==self.transcript_id,
                Structure.id==StructurePosition.structure_id
            ) \
            .order_by(Structure.structure_prediction_run_id, Structure.id, StructurePosition.position) \
            .all()

        buf = ""
        for result in results:
            structure = result[0]
            structure_position = result[1]
            run = result[2]

            buf +=  str(structure.id)+"\t"+ \
                    str(run.description)+"\t"+ \
                    str(structure.transcript_id)+"\t"+ \
                    str(structure.energy)+"\t"+ \
                    str(structure.pc1)+"\t"+ \
                    str(structure.pc2)+"\t"+ \
                    str(structure_position.letter)+"\t"+ \
                    str(structure_position.position)+"\t"+ \
                    str(structure_position.paired_to_position)+"\n"

        return buf

# Generates plain text nucleotide measurements for user download
class NucleotideMeasurementDownloader():
    def __init__(self, experiment_id, transcript_id):
        self.experiment_id = experiment_id
        self.transcript_id = transcript_id

    def generateTxt(self):
        strain_id = settings.reference_strain_id

        # Grab sequence string
        seq_str = str(Transcript(self.transcript_id) \
            .get_sequence(settings.reference_strain_id).seq)

        # Use the ORM to grab all the measurements
        results = db_session \
            .query(NucleotideMeasurement) \
            .filter(
                NucleotideMeasurement.experiment_id==self.experiment_id,
                NucleotideMeasurement.transcript_id==self.transcript_id
            ) \
            .all()
        
        # index measurements by pos
        measurements = {}
        for result in results:
            measurements[result.position] = result.measurement

        # build the output string
        buf = ""
        n = 0
        for n in range(0, len(seq_str)): 
            pos = n + 1
            measurement = "NA" if pos not in measurements else measurements[pos]
            buf +=  str(pos)+"\t"+ \
                    seq_str[n]+"\t"+ \
                    str(measurement)+"\n"
            n += 1

        return buf

