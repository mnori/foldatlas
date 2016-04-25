
# Database package
# Includes all the code needed to import plaintext files into the DB.
# @author Matthew Norris <matthew.norris@jic.ac.uk>

from sqlalchemy import create_engine, and_
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql import and_
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

import settings, os
import sys
import re
import math

# from app import sqla
engine = create_engine(settings.database_uri, convert_unicode=True)

# Autoflush = true is important to prevent errors on EC2
db_session = scoped_session(sessionmaker(autocommit=False, autoflush=True, bind=engine))
# Base = sqla.Model # declarative_base()
# Base.query = db_session.query_property()
