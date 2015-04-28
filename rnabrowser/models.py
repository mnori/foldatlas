from sqlalchemy import Column, Integer, String, Text
from database import Base

class Sequence(Base):
    # Primary key ID for the sequence

    __tablename__ = 'sequence'

    id = Column(Integer, primary_key=True, autoincrement=True)

    # The actual sequence string - lower and uppercase letters
    sequence = Column(Text, nullable = False)

    # Strain identifier
    # strain_id = db.Column(db.Integer, nullable = False, )

    # Chromosomal coordinates, relative to strain of interest's sequence
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
