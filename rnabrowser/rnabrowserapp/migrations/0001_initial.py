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
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
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
                ('id', models.CharField(max_length=255, serialize=False, primary_key=True, verbose_name='accession number')),
                ('gi_number', models.CharField(max_length=255)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TranscriptSequence',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('sequence', models.TextField()),
                ('strain', models.ForeignKey(to='rnabrowserapp.Strain')),
                ('transcript', models.ForeignKey(to='rnabrowserapp.Transcript')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
