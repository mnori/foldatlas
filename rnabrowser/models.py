# Schema definitions for RNA browser database.
# @author Matthew Norris <matthew.norris@jic.ac.uk

from sqlalchemy import Column, Integer, String, Text, Enum, ForeignKey
from database import Base

# A Gene describes a locus identifier for a gene, plus any metadata associated with the locus.
class Gene(Base):
    __tablename__ = 'gene'
    id = Column(String(256), primary_key=True)

# A Transcript describes an RNA molecule generated from a Gene sequence.
class Transcript(Base):

    __tablename__ = 'transcript'
    id = Column(String(256), primary_key=True)
    gene_id = Column(String(256), ForeignKey('gene.id'), nullable=False)

# Describes the strain - can have different sequences for each gene/strain combination
class Strain(Base):
    __tablename__ = 'strain'
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(256), nullable=False)
    description = Column(Text, nullable=False)

    def __init__(self, name=None, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return '<Strain %r>' % (self.id)


    # we probably want more stuff here - e.g. description and other gene metadata.

# Contains the sequence of a Gene or intergenic region. Each Sequence is associated with 1 strain.
class Sequence(Base):
    __tablename__ = 'sequence'

    # about 1/2 of sequences will be intergenic. so we need a unique ID  can't rely on a composite strain/transcript ID
    id = Column(Integer, primary_key=True)

    # which Strain the sequence belongs to
    strain_id = Column(Integer, ForeignKey('strain.id'), nullable=False)

    # If gene_id is missing, it means it's an intergenic sequence.
    gene_id = Column(String(256), ForeignKey('gene.id'), nullable=True)

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

class 

# note - if a sequence straddles a gene boundary, it will be troublesome (does this even happen?)
# but it's best to break it down into Sequence objects, for the sake of database design / querying etc.
class SequenceFeature(Base):
    __tablename__ = 'sequence_feature'

    id = Column(Integer, primary_key=True)

    # which Sequence entity this Feature belongs to
    sequence_id = Column(Integer, ForeignKey('sequence.id'))

    # these coords are relative to the parent Sequence entity
    start = Column(Integer, nullable = False)
    end = Column(Integer, nullable = False)

    # TODO add an enum describing what kind of feature it is.
