# Schema definitions for RNA browser database.
# @author Matthew Norris <matthew.norris@jic.ac.uk

from sqlalchemy import Column, Integer, String, Text, ForeignKey
from database import Base

# Describes the strain - can have different sequences for each gene/strain combination
class Strain(Base):
    __tablename__ = 'strain'
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(256), nullable=False)
    description = Column(Text, nullable=False)

class Transcript(Base):
    __tablename__ = 'transcript'
    id = Column(String(256), primary_key=True)
    # we probably want more stuff here - e.g. description and other gene metadata.

# Describes the sequence of an RNA gene for a particular strain
class Sequence(Base):
    __tablename__ = 'sequence'

    # about 1/2 of sequences will be intergenic. so we need a unique ID 
    # can't rely on a composite strain/transcript ID
    id = Column(Integer, primary_key=True)

    strain_id = Column(Integer, ForeignKey('strain.id'), nullable=False)
    transcript_id = Column(String(256), ForeignKey('transcript.id'), nullable=False)

    # The actual sequence string - lower and uppercase letters
    sequence = Column(Text, nullable = False)

    # Chromosomal coords, relative to strain of interest's sequence
    chromosome = Column(Integer, nullable = False)
    start = Column(Integer, nullable = False)
    end = Column(Integer, nullable = False)


    def __init__(self, sequence=None, chromosome=None, start=None, end=None):
        self.sequence = sequence
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __repr__(self):
        return '<Sequence %r>' % (self.id)

# note - if a sequence straddles a gene boundary, it will be troublesome (does this even happen?)
# but it's best to break it down into Sequence objects, for the sake of database design / querying etc.
class SequenceFeature(Base):
    __tablename__ = 'sequence_feature'

    id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, ForeignKey('sequence.id'))

    # these coords are relative to the parent Sequence entity
    start = Column(Integer, nullable = False)
    end = Column(Integer, nullable = False)
