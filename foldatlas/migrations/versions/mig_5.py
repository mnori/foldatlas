"""empty message

Revision ID: bd9b4f637408
Revises: 09ddd9e7a518
Create Date: 2016-04-21 17:38:51.183670

"""

# revision identifiers, used by Alembic.
revision = 'bd9b4f637408'
down_revision = '09ddd9e7a518'

from alembic import op
import sqlalchemy as sa
from importers import BppmImporter

# This script inserts the BPPM data. Takes a long time!

def upgrade():
	print("upgrade() invoked")
	BppmImporter().run()
	# exit() # don't forget to remove this!

def downgrade():
	print("Nothing to do for BPPM downgrade")
    
