# the import will fail unless we tell python3 where to find app.py
import sys
sys.path.insert(0, '/vagrant/rnabrowser')

# now do the import
from app import app as application
