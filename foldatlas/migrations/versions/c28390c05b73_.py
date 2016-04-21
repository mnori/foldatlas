"""empty message

Revision ID: c28390c05b73
Revises: 6fe4187803d7
Create Date: 2016-04-21 14:05:22.679010

"""

# revision identifiers, used by Alembic.
revision = 'c28390c05b73'
down_revision = '6fe4187803d7'

from alembic import op
import sqlalchemy as sa


def upgrade():
	from importers import MinusPlusCompiler
	MinusPlusCompiler().run()

def downgrade():
	# not implemented, since we won't ever need to revert to the old counts
	print("WARNING: Undoing minus plus counts is not implemented.")
	