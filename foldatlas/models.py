# Schema definitions for RNA browser database.
# @author Matthew Norris <matthew.norris@jic.ac.uk

from sqlalchemy import Column, Integer, String, Text, Enum, Float, ForeignKey, ForeignKeyConstraint
from database import Base
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import database
import settings

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

    # Retrieve sequences for a transcript, keyed by strain ID.
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

            # APPEND the feature sequence
            transcript_seqs[strain_id]["seq"] += row["seq"]
            transcript_seqs[strain_id]["direction"] = row["direction"]

        # make collection of SeqRecord objects.
        seqs_out = {}
        for strain_id in transcript_seqs:
            seq = transcript_seqs[strain_id]["seq"]

            # if direction is reverse, do reverse complement
            if transcript_seqs[strain_id]["direction"] == "reverse":
                seq = seq.reverse_complement()

            seq = Seq(str(seq).replace("T", "U"))
            seqs_out[strain_id] = SeqRecord(seq, id=strain_id, description="")

        return seqs_out

    # convenience method to fetch a single SeqRecord sequence.
    # Sequence is always reverse complemented if it's a backwards gene
    def get_sequence(self, strain_id=None):
        if strain_id == None:
            strain_id = settings.reference_strain_id

        vals = list(self.get_sequences(strain_id).values())
        if len(vals) > 0:
            return vals[0]
        else:
            return None

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

# GeneLocation describes the location of a gene for a particular strain. This table is redundant
# since everything needed is already in the Feature table. But it is cached here for speed.
class GeneLocation(Base):
    __tablename__ = "gene_location"

    # This constraint maps the GeneLocation to a unique Chromosome entry.
    __table_args__ = (
        ForeignKeyConstraint(
            ["strain_id", "chromosome_id"],
            ["chromosome.strain_id", "chromosome.chromosome_id"]
        ),
    )

    gene_id = Column(String(256), ForeignKey("gene.id"), nullable=False, primary_key=True)
    strain_id = Column(String(256), ForeignKey("strain.id"), nullable=False, primary_key=True)
    chromosome_id = chromosome_id = Column(String(256), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    direction = Column(Enum("forward", "reverse"), nullable=False)

    def __init__(self, gene_id=None, strain_id=None, chromosome_id=None, start=None, end=None, direction=None):
        self.gene_id = gene_id
        self.strain_id = strain_id
        self.chromosome_id = chromosome_id
        self.start = start
        self.end = end
        self.direction = direction

    def __repr__(self):
        return "<GeneLocation "+self.gene_id+", "+self.strain_id+">";

class NucleotideExperiment(Base):
    __tablename__ = "nucleotide_experiment"

    id = Column(Integer, primary_key=True, autoincrement=False)
    strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=False)
    description = Column(Text, nullable=False)

    def __init__(self, id, strain_id, description):
        self.id = id
        self.strain_id = strain_id
        self.description = description

    def __repr__(self):
        return "<NucleotideExperiment %r>" % (self.id)

class StructurePredictionRun(Base):
    __tablename__ = "structure_prediction_run"

    id = Column(Integer, primary_key=True, autoincrement=False)
    strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=False)
    description = Column(Text, nullable=False)

    def __init__(self, id, strain_id, description):
        self.id = id
        self.strain_id = strain_id
        self.description = description

    def __repr__(self):
        return "<StructurePredictionRun %r>" % (self.id)

# Represents a coverage measurement for a single transcript
class NucleotideMeasurementSet(Base):
    __tablename__ = "nucleotide_measurement_set"

    id = Column(Integer, primary_key=True, autoincrement=True)
    nucleotide_experiment_id = Column(Integer, ForeignKey("nucleotide_experiment.id"))
    transcript_id = Column(String(256), ForeignKey("transcript.id"))
    coverage = Column(Float, nullable=False) 

    def __init__(self, nucleotide_experiment_id=None, transcript_id=None, coverage=None):

        self.nucleotide_experiment_id = nucleotide_experiment_id
        self.transcript_id = transcript_id
        self.coverage = coverage

    def __repr__(self):
        return "<NucleotideMeasurementSet %r-%r-%r>" % (
            self.nucleotide_experiment_id, self.transcript_id
        )

# Represents one measurement, at a particular nucleotide position.
class NucleotideMeasurement(Base):
    __tablename__ = "nucleotide_measurement"

    nucleotide_measurement_set_id = Column(Integer, 
        ForeignKey("nucleotide_measurement_set.id"), primary_key=True, autoincrement=True)

    # if there's no measurement at a position, there is no corresponding row.
    position = Column(Integer, autoincrement=False, primary_key=True) 
    measurement = Column(Float, nullable=False) 

    def __init__(self, nucleotide_measurement_set_id=None, position=None, measurement=None):
        self.nucleotide_measurement_set_id = nucleotide_measurement_set_id
        self.position = position
        self.measurement = measurement

    def __repr__(self):
        return "<NucleotideMeasurement %r-%r-%r-%r>" % (
            self.nucleotide_measurement_set_id, self.position, self.measurement
        )

# Represents a structure prediction for a single RNA sequence
class Structure(Base):

    __tablename__ = "structure"

    id = Column(Integer, primary_key=True, autoincrement=True)
    structure_prediction_run_id = Column(Integer, ForeignKey("structure_prediction_run.id"), nullable=False)
    transcript_id = Column(String(256), ForeignKey("transcript.id"), nullable=False)
    energy = Column(Float, nullable=False)
    pc1 = Column(Float, nullable=False, default=0)
    pc2 = Column(Float, nullable=False, default=0)

    def __init__(self, structure_prediction_run_id, transcript_id, energy, pc1=0, pc2=0):
        self.structure_prediction_run_id = structure_prediction_run_id
        self.transcript_id = transcript_id
        self.energy = energy
        self.pc1 = pc1
        self.pc2 = pc2

    def __repr__(self):
        return "<Structure %r>" % (self.id)

# Represents a single position within an RNA sequence structure prediction
class StructurePosition(Base):

    __tablename__ = "structure_position"

    structure_id = Column(Integer, ForeignKey("structure.id"), primary_key=True)
    position = Column(Integer, primary_key=True, autoincrement=False)
    paired_to_position = Column(Integer, nullable=True) # null means unpaired
    letter = Column(Enum("T", "A", "C", "G"), nullable=False)

    def __init__(self, structure_id, position, paired_to_position, letter):
        self.structure_id = structure_id
        self.position = position
        self.paired_to_position = paired_to_position
        self.letter = letter

    def __repr__(self):
        return "<StructurePosition %r-%r>" % (self.structure_id, self.position)

# Represents a single row of ClustalW alignment data. 
# Each row is a combination of strain and transcript.
# TODO add type field(s), so we can store splices/unspliced/different methods etc.
# class AlignmentEntry(Base):
#     __tablename__ = "alignment_entry"

#     transcript_id = Column(String(256), ForeignKey("transcript.id"), primary_key=True)
#     strain_id = Column(String(256), ForeignKey("strain.id"), primary_key=True)
#     sequence = Column(Text, nullable=False)

#     def __init__(self, transcript_id=None, strain_id=None, sequence=None):
#         self.transcript_id = transcript_id
#         self.strain_id = strain_id
#         self.sequence = sequence

#     def __repr__(self):
#         return "<AlignmentEntry %r-%r>" % (self.transcript_id, self.strain_id)