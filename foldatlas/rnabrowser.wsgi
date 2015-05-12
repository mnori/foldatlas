# the import will fail unless we tell python3 where to find app.py
# note - this is not used in the development environment, since it's better to use 
# Flask's built in web server for that.

import sys
sys.path.insert(0, '/vagrant/foldatlas')

# now do the import
from app import app as application
