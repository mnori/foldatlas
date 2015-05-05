from database import Feature, Transcript, db_session;
from sqlalchemy import and_
import json
import database

# Fetches sequence annotation data from the DB and sends it to the genome
# browser front end as JSON.

class GenomeBrowser():

    def get_transcripts(self, request):

        print("request")
        print(request)

        out = []
        chromosome_id = "Chr"+str(int(request.args.get('chr'))) # SQL-injection safe
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        print("chr: "+chromosome_id)
        print("start: "+str(start))
        print("end: "+str(end))

        sql =   ("SELECT *, MIN(start) min_start, MAX(end) max_end "
                 "FROM feature "
                 "WHERE chromosome_id = '"+chromosome_id+"' "
                 "AND start > '"+str(start)+"' "
                 "AND end < '"+str(end)+"' "
                 "GROUP BY transcript_id")

        # demonstrates how genes would be done ###########################
        results = database.engine.execute(sql)
        for result in results:
            out.append({
                "Parent": result.transcript_id, # actually use gene ID instead
                "feature_type": "transcript", # without this, it won't draw
                "logic_name": "ensembl_havana", # does not work without this label
                "start": result.min_start,
                "end": result.max_end,
                "id": result.transcript_id,
                "strand": 1, # whether it is + or -??
            })

        results =  db_session.query(Feature).filter(and_(Feature.start >= start, Feature.end <= end, Feature.chromosome_id == chromosome_id)).all() 
        for feature in results:
            
            out.append({
                "Parent": feature.transcript_id, # actually use gene ID instead
                "feature_type": feature.type_id, # without this, it won't draw the gene
                "logic_name": "ensembl_havana", # does not work without this label
                "start": feature.start,
                "end": feature.end,
                "id": feature.transcript_id+"-"+str(feature.id),
                "strand": 1, # whether it is + or -??
            })


        buf = json.dumps(out)
        return buf

    