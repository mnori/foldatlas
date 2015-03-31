from django.db import models

# Describes the strain - can have different sequences for each 
# gene/strain combination
class Strain(models.Model):
    # id is default auto increment
    name = models.CharField(max_length=255)
    description = models.TextField()

# A Transcript describes an RNA gene (can be shared across multiple Strains)
class Transcript(models.Model):
    # note: this ID is the accession number (accession_number)
    id = models.CharField("accession number", max_length=255, primary_key=True) 
    gi_number = models.CharField(max_length=255)
    def get_absolute_url(self):
    	return "/transcript/"+self.id+"/"

# Describes the sequence of an RNA gene for a particular strain
# note - this uses the default "id" key. would be better to use a composite
# of [rna_gene_id, strain_id], but django doesn't support that :(
class TranscriptSequence(models.Model):
    strain = models.ForeignKey(Strain)
    transcript = models.ForeignKey(Transcript)
    sequence = models.TextField()
