from database import Feature, Transcript, db_session;
from sqlalchemy import and_
import json, database, settings

# Fetches sequence annotation data from the DB and sends it to the genome
# browser front end as JSON.

class GenomeBrowser():

    def get_transcripts(self, request):

        out = []

        chromosome_id = "Chr"+str(int(request.args.get('chr'))) # SQL-injection safe
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        sql =   ("SELECT *, MIN(start) min_start, MAX(end) max_end "
                 "FROM feature "
                 "WHERE chromosome_id = '"+chromosome_id+"' "
                 "AND start > '"+str(start)+"' "
                 "AND end < '"+str(end)+"' "
                 "GROUP BY transcript_id")

        # Add the transcript rows to the output
        results = database.engine.execute(sql)
        for result in results:
            out.append({
                "Parent": result.transcript_id,
                "feature_type": "transcript", # without this, it won't draw
                "direction": result.direction,
                "start": result.min_start,
                "end": result.max_end,
                "id": result.transcript_id
            })

        # Use the ORM to get feature details
        results = db_session \
            .query(Feature) \
            .filter(and_( \
                Feature.start >= start, Feature.end <= end, Feature.chromosome_id == chromosome_id)) \
            .all() 

        # Add transcript feature rows to the output
        for feature in results:
            out.append({
                "Parent": feature.transcript_id,
                "feature_type": feature.type_id,
                "direction": result.direction,
                "start": feature.start,
                "end": feature.end,
                "id": feature.transcript_id+"-"+str(feature.id)
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
                "AND chromosome_id = '"+chromosome_id+"' "
                "AND start > '"+str(start)+"' "
                "AND end < '"+str(end)+"' "
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

class AlignmentViewer():

    def get_alignment_entries(self, transcript_id):

        # fetch the alignment rows from the DB, using the ORM
        alignment_entries = db_session \
            .query(AlignmentEntry) \
            .filter(transcript_id=transcript_id) \
            .all() 

        alignments = {}
        for alignment_entry in alignment_entries:
            alignments[alignment_entry] = alignment_entry

        return alignments
