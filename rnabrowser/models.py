# Schema definitions for RNA browser database.
# @author Matthew Norris <matthew.norris@jic.ac.uk

from sqlalchemy import Column, Integer, String, Text, Enum, ForeignKey
from database import Base


# Describes a strain - can have different sequences for each gene/strain combination
class Strain(Base):
    __tablename__ = 'strain'
    id = Column(String(256), nullable=False, primary_key=True)
    description = Column(Text, nullable=False)

    def __init__(self, id=None, description=None):
        self.id = id
        self.description = description

    def __repr__(self):
        return '<Strain %r>' % (self.id)

# A Gene describes a locus identifier for a gene, plus any metadata associated with the locus.
class Gene(Base):
    __tablename__ = 'gene'

    id = Column(String(256), primary_key=True) # TAIR locus ID (e.g. AT1G01225)

    def __init__(self, id=None):
        self.id = id

    def __repr__(self):
        return '<Gene %r>' % (self.id)


class Transcript(Base):
    __tablename__ = 'transcript'

    id = Column(String(256), primary_key=True) # TAIR transcript ID (e.g. AT1G01225.1)
    gene_id = Column(String(256), ForeignKey('gene.id'), nullable=False)

    def __init__(self, id=None, gene_id=None):
        self.id = id
        self.gene_id = gene_id

    def __repr__(self):
        return '<Transcript %r>' % (self.id)

# Contains the sequence of a Gene or intergenic region. Each Sequence is associated with 1 strain.
class Sequence(Base):
    __tablename__ = 'sequence'

    # Each Sequence gets a unique Sequence ID
    id = Column(Integer, nullable=False, primary_key=True)

    # which Strain the sequence belongs to
    strain_id = Column(String(256), ForeignKey('strain.id'), nullable=False)

    # gene_id is NULL when it's an intergenic sequence.
    gene_id = Column(String(256), ForeignKey('gene.id'), nullable=True)

    # The actual sequence string - lower and uppercase letters - lowercase means unmapped and imputed from the reference strain.
    sequence = Column(Text, nullable = False)

    # Chromosomal coords, relative to strain of interest's sequence
    chromosome = Column(String(256), nullable = False)
    start = Column(Integer, nullable = False)
    end = Column(Integer, nullable = False)
    
    def __init__(self, id=None, strain_id=None, gene_id=None, sequence=None, chromosome=None, start=None, end=None):
        self.id = id
        self.strain_id = strain_id
        self.gene_id = gene_id
        self.sequence = sequence
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __repr__(self):
        return '<Sequence %r>' % (self.id)

# need TranscriptSequence to store the direction.
class TranscriptSequence(Base):
    __tablename__ = 'transcript_sequence'

    sequence_id = Column(Integer, ForeignKey('sequence.id'), primary_key=True)

    transcript_id = Column(String(256), ForeignKey('transcript.id'), primary_key=True)

    # these coords are relative to the associated Sequence entry
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)

    direction = Column(Enum("forward", "reverse"), nullable = False)

    def __init__(self, transcript_id=None, sequence_id=None, start=None, end=None, direction=None):
        self.transcript_id = transcript_id
        self.sequence_id = sequence_id
        self.start = start
        self.end = end
        self.direction = direction

    def __repr__(self):
        return '<TranscriptSequence %r>' % (self.id)


# TranscriptSequenceFeatures annotations of the TranscriptSequence. These can be either UTRs, CDSs or miRNAs.
# Uncovered regions denote introns. To get the RNA sequence of a Transcript, simply concatenate the 
# sequences of the Features together.
class TranscriptSequenceFeature(Base):
    __tablename__ = "transcript_sequence_feature"

    id = Column(Integer, primary_key=True)
    
    sequence_id = Column(Integer, ForeignKey('sequence.id'), nullable=False)

    # Transcript identifier - this is a string
    transcript_id = Column(String(256), ForeignKey('transcript.id'), nullable=False)

    # These coords are relative to the parent TranscriptSequence entity's coords.
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)

    type_id = Column(String(256), nullable=False)

    def __init__(self, sequence_id=None, transcript_id=None, start=None, end=None, type_id=None):
        self.sequence_id = sequence_id
        self.transcript_id = transcript_id
        self.start = start
        self.end = end
        self.type_id = type_id

    def __repr__(self):
        return '<TranscriptSequenceFeature %r>' % (self.id)

# note - if a sequence straddles a gene boundary, it will be troublesome (does this even happen?)
# but it's best to break it down into Sequence objects, for the sake of database design / querying etc.
# class SequenceFeature(Base):
#     __tablename__ = 'sequence_feature'
#     id = Column(Integer, primary_key=True)

#     # which Sequence entity this Feature belongs to
#     sequence_id = Column(Integer, ForeignKey('sequence.id'))

#     # these coords are relative to the parent Sequence entity
#     start = Column(Integer, nullable = False)
#     end = Column(Integer, nullable = False)

    # TODO add an enum describing what kind of feature it is.


# A TranscriptSequence describes the coordinates of a particular transcript, relative to Sequence coords.
# This is basically everything between the transcription start and end sites, including introns.
# Not actually needed - can simply use TranscriptSequenceFeature to do queries
# class TranscriptSequence(Base):
#     transcript_id = Column(String(256), ForeignKey('transcript.id'), nullable=False, primary_key=True)
#     sequence_id = Column(Integer, ForeignKey('sequence.id'), primary_key=True)
#     start = Column(Integer, nullable = False)
#     end = Column(Integer, nullable = False)