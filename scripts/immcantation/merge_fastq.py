#!/usr/bin/env python3
"""
Concatentate multiple fastq files
"""
from sys import argv
from Bio import SeqIO

out_file = argv[1]
in_files = argv[2:]

with open(out_file, 'w') as out_handle:
    for f in in_files:
        records = SeqIO.parse(f, 'fastq')
        SeqIO.write(records, out_handle, 'fastq')

#print(out_file)