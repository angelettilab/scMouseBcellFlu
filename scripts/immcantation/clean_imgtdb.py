#!/usr/bin/env python3
"""
Clean IMGT germline fasta files for IgBLAST database build
"""
from sys import argv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Get input and output file names
in_file = argv[1]
out_file = argv[2]

# Load sequences into memory and process them
name_set = set()
seq_list = list()
for rec in SeqIO.parse(in_file, 'fasta'):
    name = rec.description.split('|')[1]
    if name not in name_set:
        name_set.add(name)
        seq = SeqRecord(rec.seq.ungap('.').upper(), id=name, name=name, description=name)
        seq_list.append(seq)

# Overwrite file
with open(out_file, 'w') as out_handle:
    writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
    writer.write_file(seq_list)
