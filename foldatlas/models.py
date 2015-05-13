# Schema definitions for RNA browser database.
# @author Matthew Norris <matthew.norris@jic.ac.uk

from sqlalchemy import Column, Integer, String, Text, Enum, Float, ForeignKey, ForeignKeyConstraint
from database import Base
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import database

# A Gene describes a locus identifier for a gene, plus any metadata associated with the locus.
# Genes are generic - they can be associated with multiple strains.
class Gene(Base):
    __tablename__ = "gene"

    id = Column(String(256), primary_key=True) # TAIR locus ID (e.g. AT1G01225)

    def __init__(self, id=None):
        self.id = id

    def __repr__(self):
        return "<Gene %r>" % (self.id)

# A Transcript is effectively an RNA sequence identifier, which can be shared amongst multiple strains.
# Sequences are mapped to Transcripts via the Feature entity.
class Transcript(Base):
    __tablename__ = "transcript"

    id = Column(String(256), primary_key=True) # TAIR transcript ID (e.g. AT1G01225.1)
    gene_id = Column(String(256), ForeignKey("gene.id"), nullable=False)

    def __init__(self, id=None, gene_id=None):
        self.id = id
        self.gene_id = gene_id

    def __repr__(self):
        return "<Transcript %r>" % (self.id)

    def get_sequences(self, strain_id=None):

        if strain_id != None:
            strain_sql = "AND feature.strain_id = :strain_id"
        else:
            strain_sql = ""


        # given the transcript ID, fetch the feature sequences in the correct order.
        # TODO wrap this up - have a Transcript.get_sequences() method
        sql = """
            SELECT 
                feature.strain_id, 
                feature.direction,
                SUBSTR(
                    chromosome.sequence, 
                    feature.start,
                    feature.end - feature.start + 1
                ) seq
            FROM chromosome, feature
            WHERE 
                feature.strain_id = chromosome.strain_id AND
                feature.chromosome_id = chromosome.chromosome_id AND
                feature.type_id = 'exon' 
                AND feature.transcript_id = :transcript_id
                {0}
            ORDER BY feature.strain_id, start
        """
        sql = sql.format(strain_sql)

        sql_params = {"transcript_id": str(self.id)}
        if strain_id != None:
            sql_params["strain_id"] = strain_id

        results = database.db_session.execute(sql, sql_params)

        # collect data about the sequences
        transcript_seqs = {}
        for row in results:
            strain_id = row["strain_id"]
            # print("Found ["+strain_id+"]")
            if strain_id  not in transcript_seqs:
                transcript_seqs[strain_id] = {} 
                transcript_seqs[strain_id]["seq"] = Seq("")

            transcript_seqs[strain_id]["seq"] += row["seq"]
            transcript_seqs[strain_id]["direction"] = row["direction"]

        # make collection of SeqRecord objects.
        seqs_out = {}
        for strain_id in transcript_seqs:
            seq = transcript_seqs[strain_id]["seq"]

            # if direction is reverse, do reverse complement
            if transcript_seqs[strain_id]["direction"] == "reverse":
                seq.reverse_complement()

            seqs_out[strain_id] = SeqRecord(seq, id=strain_id, description="")

        return seqs_out

    # convenience method to fetch a single SeqRecord.
    def get_sequence(self, strain_id):
        return list(self.get_sequences(strain_id).values())[0]

# Describes a strain.
class Strain(Base):
    __tablename__ = "strain"
    id = Column(String(256), nullable=False, primary_key=True)
    description = Column(Text, nullable=False)

    def __init__(self, id=None, description=None):
        self.id = id
        self.description = description

    def __repr__(self):
        return "<Strain %r>" % (self.id)

# Describes the sequence of a chromosome for a particular strain. This is the only place
# where nucleotide sequence data is stored.
class Chromosome(Base):
    __tablename__ = "chromosome"

    strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=True)

    chromosome_id = Column(String(256), primary_key=True)

    sequence = Column(Text(4294967295), nullable=False)

    def __init__(self, strain_id=None, chromosome_id=None, sequence=None):
        self.strain_id = strain_id
        self.chromosome_id = chromosome_id
        self.sequence = sequence

    def __repr__(self):
        return "<Chromosome "+self.strain_id+", "+self.chromosome_id+" >"

# TranscriptSequenceFeatures annotations of the ChromosomeSequence. This is the main destination of 
# all the *.gff3 data.
class Feature(Base):
    __tablename__ = "feature"

    # This constraint maps the Feature to a unique Chromosome entry.
    __table_args__ = (
        ForeignKeyConstraint(
            ["strain_id", "chromosome_id"],
            ["chromosome.strain_id", "chromosome.chromosome_id"]
        ),
    )

    # A unique identifier for this feature.
    id = Column(Integer, primary_key=True)

    # Transcript identifier - this is a string
    transcript_id = Column(String(256), ForeignKey("transcript.id"), nullable=False)

    # What kind of annotation is this? Maybe have a foreign key pointing to a special meta table?
    # Or use an enum
    type_id = Column(String(256), nullable=False)

    # These properties describe where we can find the associated sequence.
    strain_id = Column(String(256), ForeignKey("strain.id"), nullable=False)
    chromosome_id = Column(String(256), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    direction = Column(Enum("forward", "reverse"), nullable=False)

    def __init__(self, transcript_id=None, type_id=None, strain_id=None, chromosome_id=None, start=None, end=None, direction=None):
        self.transcript_id = transcript_id
        self.type_id = type_id
        self.strain_id = strain_id
        self.chromosome_id = chromosome_id
        self.start = start
        self.end = end
        self.direction = direction

    def __repr__(self):
        return "<Feature %r>" % (self.id)

# Represents a single row of ClustalW alignment data. 
# Each row is a combination of strain and transcript.
# TODO add type field(s), so we can store splices/unspliced/different methods etc.
class AlignmentEntry(Base):
    __tablename__ = "alignment_entry"

    transcript_id = Column(String(256), ForeignKey("transcript.id"), primary_key=True)
    strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=True)
    sequence = Column(Text, nullable=False)

    def __init__(self, transcript_id=None, strain_id=None, sequence=None):
        self.transcript_id = transcript_id
        self.strain_id = strain_id
        self.sequence = sequence

    def __repr__(self):
        return "<AlignmentEntry %r-%r>" % (self.transcript_id, self.strain_id)

# Represents one reactivity measurement, at a particular nucleotide position.
class ReactivityMeasurement(Base):
    __tablename__ = "reactivity_measurement"

    strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=True)
    transcript_id = Column(String(256), ForeignKey("transcript.id"), primary_key=True)
    position = Column(Integer, autoincrement=False, primary_key=True) # if there's no reactivity at a position, there is no corresponding row.
    reactivity = Column(Float, nullable=False) 

    def __init__(self, strain_id=None, transcript_id=None, position=None, reactivity=None):

        self.strain_id = strain_id
        self.transcript_id = transcript_id
        self.position = position
        self.reactivity = reactivity

    def __repr__(self):
        return "<ReactivityMeasurement %r-%r-%r>" % (self.strain_id, self.transcript_id, self.position)

    # 
    # strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=True)
    # sequence = Column(Text, nullable=False)

    # def __init__(self, transcript_id=None, strain_id=None, sequence=None):
    #     self.transcript_id = transcript_id
    #     self.strain_id = strain_id
    #     self.sequence = sequence

    # def __repr__(self):
    #     return "<AlignmentEntry %r-%r>" % (self.transcript_id, self.strain_id)

