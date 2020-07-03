#!/usr/bin/env python3
"""
Join the 10X annotation table and Change-O table
"""

# Imports
import pandas as pd
import sys

# Parse arguments
changeo_file = sys.argv[1]
annotation_file = sys.argv[2]
out_file = sys.argv[3]

# Merge annotations
annfile = pd.read_csv(annotation_file).set_index('contig_id')
annfile.columns = [column.upper() + '_10X' for column in annfile.columns]

pd.read_csv(changeo_file, dtype = 'object', sep = '\t').set_index('SEQUENCE_ID')\
.join(annfile, lsuffix='', rsuffix='_10X')\
.to_csv(out_file, sep = '\t')