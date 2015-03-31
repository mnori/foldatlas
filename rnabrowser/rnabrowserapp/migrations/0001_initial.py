# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Strain',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('description', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Transcript',
            fields=[
                ('id', models.CharField(serialize=False, max_length=255, verbose_name='accession number', primary_key=True)),
                ('gi_number', models.CharField(max_length=255)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TranscriptSequence',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('sequence', models.TextField()),
                ('strain', models.ForeignKey(to='rnabrowserapp.Strain')),
                ('transcript', models.ForeignKey(to='rnabrowserapp.Transcript')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
