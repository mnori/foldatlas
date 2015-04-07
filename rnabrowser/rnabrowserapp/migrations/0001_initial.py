# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('sequence', models.TextField()),
                ('start_ref', models.IntegerField(verbose_name='Start position relative to reference sequence')),
                ('end_ref', models.IntegerField(verbose_name='End position relative to reference sequence')),
                ('start_strain', models.IntegerField(verbose_name='Start position relative to specific strain')),
                ('end_strain', models.IntegerField(verbose_name='End position relative to specific strain')),
            ],
        ),
        migrations.CreateModel(
            name='SequenceFeature',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.CharField(max_length=256, choices=[('UTR', 'Untranslated region'), ('CDS', 'Coding sequence'), ('Intron', 'Intron')])),
                ('start_ref', models.IntegerField(verbose_name='Start position relative to reference sequence')),
                ('end_ref', models.IntegerField(verbose_name='End position relative to reference sequence')),
                ('start_strain', models.IntegerField(verbose_name='Start position relative to specific strain')),
                ('end_strain', models.IntegerField(verbose_name='End position relative to specific strain')),
                ('sequence', models.ForeignKey(to='rnabrowserapp.Sequence')),
            ],
        ),
        migrations.CreateModel(
            name='Strain',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='Transcript',
            fields=[
                ('id', models.CharField(verbose_name='TAIR gene ID', primary_key=True, serialize=False, max_length=255)),
            ],
        ),
        migrations.AddField(
            model_name='sequence',
            name='strain',
            field=models.ForeignKey(to='rnabrowserapp.Strain'),
        ),
        migrations.AddField(
            model_name='sequence',
            name='transcript',
            field=models.ForeignKey(null=True, to='rnabrowserapp.Transcript'),
        ),
    ]
