from database import db_session;
from sqlalchemy import and_

import json, database, settings
from models import Feature, Transcript, AlignmentEntry, ReactivityMeasurement

# Fetches sequence annotation data from the DB and sends it to the genome
# browser front end as JSON.

class GenomeBrowser():

    def get_transcripts(self, request):

        chromosome_id = "Chr"+str(int(request.args.get('chr'))) # SQL-injection safe
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        # need to rework this. get gene IDs first, then expand using the gene ID list.
        # using feature table and start/end, grab the gene IDs that fall within the range.
        sql = ( "SELECT DISTINCT transcript.gene_id "
                "FROM feature, transcript "
                "WHERE strain_id = '"+settings.reference_strain_id+"' "
                "AND chromosome_id = '"+chromosome_id+"' "
                "AND end > '"+str(start)+"' "
                "AND start < '"+str(end)+"' "
                "AND feature.transcript_id = transcript.id ");

        results = database.engine.execute(sql)
        gene_ids = [];
        for result in results:
            gene_ids.append(result["gene_id"])

        # use the gene IDs to get feature data
        sql = ( "SELECT * "
                "FROM transcript, feature "
                "WHERE transcript.gene_id IN ('"+"','".join(gene_ids)+"') "
                "AND strain_id = '"+settings.reference_strain_id+"' "
                "AND chromosome_id = '"+chromosome_id+"' "
                "AND transcript.id = feature.transcript_id ")

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
        
        out = []

        chromosome_id = "Chr"+str(int(request.args.get('chr'))) # SQL-injection safe
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        sql = ( "SELECT *, MIN(start) min_start, MAX(end) max_end "
                "FROM feature, transcript "
                "WHERE feature.transcript_id = transcript.id "
                "AND feature.strain_id = '"+settings.reference_strain_id+"'"
                "AND chromosome_id = '"+chromosome_id+"' "
                "AND end > '"+str(start)+"' "
                "AND start < '"+str(end)+"' "
                "GROUP BY transcript.gene_id")

        # Add gene rows to the output
        results = database.engine.execute(sql)
        for result in results:
            out.append({ 
                "feature_type": "gene", # without this, it won't draw
                "direction": result.direction,
                "id": result.gene_id,
                "start": result.min_start,
                "end": result.max_end,
                "strand": 1, # whether it is + or -?? this isn't actually used by the looks of things
            })        

        return json.dumps(out)

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
        self.reactivities_view = ReactivitiesView(self.transcript_id, settings.reference_strain_id)
        self.alignment_view = AlignmentView(self.transcript_id)

class ReactivitiesView():

    def __init__(self, transcript_id, strain_id):
        self.transcript_id = transcript_id
        self.strain_id = strain_id
        self.build_reactivity_entries()

    def build_reactivity_entries(self):
        seq_str = str(Transcript(self.transcript_id).get_sequence(self.strain_id).seq)
        reactivities = db_session \
            .query(ReactivityMeasurement) \
            .filter(
                ReactivityMeasurement.strain_id==self.strain_id,
                ReactivityMeasurement.transcript_id==self.transcript_id
            ) \
            .all()

        reactivity_data = []
        for n in range(len(seq_str)): # initialise the array
            reactivity_data.append({
                "position": n,
                "nuc": seq_str[n],
                "reactivity": None
            })

        for reactivity in reactivities: # add values where present
            reactivity_data[reactivity.position - 1]["reactivity"] = reactivity.reactivity

        self.reactivity_data_json = json.dumps(reactivity_data)
        
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
