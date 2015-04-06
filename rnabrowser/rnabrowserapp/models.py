# @author Matthew Norris

from django.db import models

# Describes the strain - can have different sequences for each 
# gene/strain combination
class Strain(models.Model):
    # id is default auto increment
    name = models.CharField(max_length=255)
    description = models.TextField()

# A gene is a generalised identifier for a gene sequence - can have multiple sequences
# associated, one for each strain.
class Transcript(models.Model):
    # note: this ID is the accession number (accession_number)
    id = models.CharField("TAIR gene ID", max_length=255, primary_key=True) 
    def get_absolute_url(self):
    	return "/transcript/"+self.id+"/"

# Describes the sequence of an RNA gene for a particular strain
# note - this uses the default "id" key. would be better to use a composite
# of [rna_gene_id, strain_id], but django doesn't support that :(
class Sequence(models.Model):
    strain = models.ForeignKey(Strain)
    transcript = models.ForeignKey(Transcript, null=True)
    sequence = models.TextField()
    start_ref = models.IntegerField("Start position relative to reference sequence")
    end_ref = models.IntegerField("End position relative to reference sequence")
    start_strain = models.IntegerField("Start position relative to specific strain")
    end_strain = models.IntegerField("End position relative to specific strain")

# note - if a sequence straddles a gene boundary, it will be troublesome
# but it's best to break it down into Sequence objects, for the sake of database
# design / querying etc.
class SequenceFeature(models.Model):
    sequence = models.ForeignKey(Sequence)
    type = models.CharField(max_length=256, choices=(
        ("UTR", "Untranslated region"),
        ("CDS", "Coding sequence"),
        ("Intron", "Intron")
    ))

     ## enter enum values here
    start_ref = models.IntegerField("Start position relative to reference sequence")
    end_ref = models.IntegerField("End position relative to reference sequence")
    start_strain = models.IntegerField("Start position relative to specific strain")
    end_strain = models.IntegerField("End position relative to specific strain")