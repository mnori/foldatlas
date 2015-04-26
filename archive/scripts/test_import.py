# Test FASTA import class
import sys, os
sys.path.append('/vagrant/rnabrowser')
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
from django.conf import settings
import rnabrowser.settings as settings

from rnabrowserapp.lib import importers
# print("sys.version"+sys.version);

importer = importers.Importer(settings)
importer.do_import()



# from rna_browpwdser_app.lib import importers

# import sys, os
# sys.path.append('/path/to/your/django/app')
# os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
# from django.conf import settings

# BASE_PATH = expanduser("~")+"/ding/browser/temp/ATH_EST_sequences_20101108"
# SEQS_PATH = BASE_PATH+"/ding/browser/temp"
# SEQUENCE_SOURCES = {
#     "Col-0": SEQS_PATH+"/ding/browser/temp/ATH_EST_sequences_20101108.fas"
# }

# importer = Importer()
# importer.do_import()

# wat happen?

