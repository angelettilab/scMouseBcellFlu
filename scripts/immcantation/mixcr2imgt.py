#!/usr/bin/env python3
"""
Add IMGT gaps to the V-region of MiXCR exported clones
"""

# Imports
import csv
import sys
import pandas as pd
from changeo.IO import readGermlines
from changeo.Gene import getGene
from changeo.Alignment import gapV

# Parse arguments
clone_file = sys.argv[1]
repo_file = sys.argv[2]
out_file = sys.argv[3]

# Load data
repo_dict = readGermlines(repo_file)
data = pd.read_table(clone_file, low_memory=False)

# Extract MiXCR alignment reference points
# From https://mixcr.readthedocs.io/en/latest/export.html
anchor_regex = '^(?P<V5UTRBegin>-?[0-9]*):' \
               '(?P<L1Begin>-?[0-9]*):' \
               '(?P<VIntronBegin>-?[0-9]*):' \
               '(?P<L2Begin>-?[0-9]*):' \
               '(?P<FR1Begin>-?[0-9]*):' \
               '(?P<CDR1Begin>-?[0-9]*):' \
               '(?P<FR2Begin>-?[0-9]*):' \
               '(?P<CDR2Begin>-?[0-9]*):' \
               '(?P<FR3Begin>-?[0-9]*):' \
               '(?P<CDR3Begin>-?[0-9]*):' \
               '(?P<V3Deletion>-?[0-9]*):' \
               '(?P<VEnd>-?[0-9]*):' \
               '(?P<DBegin>-?[0-9]*):' \
               '(?P<D5Deletion>-?[0-9]*):' \
               '(?P<D3Deletion>-?[0-9]*):' \
               '(?P<DEnd>-?[0-9]*):' \
               '(?P<JBegin>-?[0-9]*):' \
               '(?P<J5Deletion>-?[0-9]*):' \
               '(?P<CDR3End>-?[0-9]*):' \
               '(?P<FR4End>-?[0-9]*):' \
               '(?P<CBegin>-?[0-9]*):' \
               '(?P<CExon1End>-?[0-9]*)$'
data = pd.concat([data, data.refPoints.str.extract(anchor_regex, expand=True).apply(pd.to_numeric)], axis=1)

# Extract MiXCR V alignment positions
valign_regex = '^(?P<targetFrom>[0-9]*)\|' \
               '(?P<targetTo>[0-9]*)\|' \
               '(?P<targetLength>[0-9]*)\|' \
               '(?P<queryFrom>[0-9]*)\|' \
               '(?P<queryTo>[0-9]*)\|'
data = pd.concat([data, data.allVAlignments.str.extract(valign_regex, expand=True).apply(pd.to_numeric)], axis=1)

# Open output file
out_handle = open(out_file, 'w')
writer = csv.DictWriter(out_handle,
                        fieldnames=['SEQUENCE_IMGT', 'V_GERM_START_IMGT', 'V_GERM_LENGTH_IMGT', 'GERMLINE_IMGT_V_REGION'],
                        delimiter='\t')
writer.writeheader()

for i, row in data.iterrows():
    # Build dictionary of required fields
    v_call = '%s*01' % getGene(row['allVHitsWithScore'])
    db = {'SEQUENCE_ID': row['cloneId'],
          'SEQUENCE_VDJ': row['clonalSequence'],
          'V_GERM_START_VDJ': row['targetFrom'] + 1,
          'V_GERM_LENGTH_VDJ': row['targetTo'] - row['targetFrom'],
          'V_CALL': v_call}

    # Create and write IMGT gapped sequence
    seq = gapV(db, repo_dict)
    if seq['V_GERM_LENGTH_IMGT']:
        seq['GERMLINE_IMGT_V_REGION'] = repo_dict[v_call][:seq['V_GERM_LENGTH_IMGT']]
    writer.writerow(seq)

# Close open handles
out_handle.close()