# the import will fail unless we tell python3 where to find app.py
# note - this is not used in the development environment, since it's better to use 
# Flask's built in web server for that.

import logging, sys
logging.basicConfig(stream=sys.stderr)

sys.path.insert(0, '/var/www/foldatlas/foldatlas')

print(sys.version)

# now do the import
from app import app as application