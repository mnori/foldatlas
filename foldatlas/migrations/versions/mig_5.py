"""empty message

Revision ID: b890906ee17d
Revises: 09ddd9e7a518
Create Date: 2016-04-24 13:41:19.526609

"""

# revision identifiers, used by Alembic.
revision = 'b890906ee17d'
down_revision = '09ddd9e7a518'

from alembic import op
import sqlalchemy as sa
from importers import BppmImporter

def upgrade():
	BppmImporter().run()
	exit()

def downgrade():
	print("Nothign to do for BPPM import reversal")
	pass
	
