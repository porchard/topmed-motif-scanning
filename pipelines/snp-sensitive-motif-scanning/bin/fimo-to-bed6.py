#!/usr/bin/env python

import sys
import csv

FIMO = sys.argv[1]

line_count = 0

with open(FIMO, 'r') as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for line in reader:
        motif_id = line['#pattern name'] if '#pattern name' in line else line['motif_id']
        chrom = line['sequence name'] if 'sequence name' in line else line['sequence_name']
        start = str(int(line['start']) - 1)
        end = str(int(line['stop']))
        out = [chrom, start, end, motif_id, line['score'], line['strand']]
        print('\t'.join(out))
