from django.core.management.base import BaseCommand, CommandError
from rnabrowserapp.lib import importers
import rnabrowser.settings as settings

class Command(BaseCommand):
    help = 'Runs some test code'

    def handle(self, *args, **options):
        importer = importers.Importer(settings)
        importer.do_import()


