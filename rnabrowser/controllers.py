
# Fetches sequence annotation data from the DB and sends it to the genome
# browser front end as JSON.

class GenomeBrowser(Base):
    __tablename__ = 'transcript'

    