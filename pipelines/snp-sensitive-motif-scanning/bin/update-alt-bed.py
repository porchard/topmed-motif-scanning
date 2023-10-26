#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd
import re

# ALT_BED = '/lab/work/porchard/2022-muscle-sn/work/snp-sensitive-motif-scanning/work/62/6220e8d00db2cf6b908f322ecd73cc/CTCF_known2.bed'
ALT_BED = sys.argv[1]

# variant positions are given 1-indexed, vcf style
# absolute start (encoded in the 'chrom' is 1-indexed, vcf style)
# motif start and end have already been adjusted to be 0-indexed, bed style

def parse_sequence_name(s):
    RE = '^(.*)_(\d+)A@(\d+)(.*)\/(.*):(\d+)-(\d+)$'
    chrom, flanking, pos, ref, alt, start, end = re.match(RE, s).groups()
    return {'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt, 'absolute_start': int(start), 'absolute_end': int(end)}



def update_coordinate(absolute_start, variant_pos, variant_ref, variant_alt, motif_coordinate):
    # if there's an insertion, then len(ref) - len(alt) = negative value, and need
    # assume that len(ref) == len(alt), then adjust as need be
    absolute_coordinate = absolute_start + motif_coordinate
    if len(variant_ref) > len(variant_alt):
        # deletion
        # so if the coordinate happens beyond the variant position, add the difference in lengths
        if absolute_coordinate > variant_pos:
            absolute_coordinate += (len(variant_ref) - len(variant_alt))
    elif len(variant_alt) > len(variant_ref):
        # insertion
        # if the coordinate happens within the variant, set to variant pos
        if absolute_coordinate > variant_pos and absolute_coordinate <= (variant_pos + len(variant_alt) - len(variant_ref)):
            absolute_coordinate = variant_pos
        # if the coordinate happens beyond the variant, sutract the difference in lengths
        if absolute_coordinate > variant_pos and absolute_coordinate > (variant_pos + len(variant_alt) - len(variant_ref)):
            absolute_coordinate -= (len(variant_alt) - len(variant_ref))
    return (absolute_coordinate - 1)


bed = pd.read_csv(ALT_BED, header=None, sep='\t', names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

motif_chrom = []
motif_start = []
motif_end = []
for chrom, start, end in zip(bed.chrom, bed.start, bed.end):
    parsed = parse_sequence_name(chrom)
    motif_start.append(update_coordinate(parsed['absolute_start'], variant_pos=parsed['pos'], variant_ref=parsed['ref'], variant_alt=parsed['alt'], motif_coordinate=int(start)))
    motif_end.append(update_coordinate(parsed['absolute_start'], variant_pos=parsed['pos'], variant_ref=parsed['ref'], variant_alt=parsed['alt'], motif_coordinate=int(end)))
    motif_chrom.append(parsed['chrom'])
bed['motif_chrom'] = motif_chrom
bed['motif_start'] = motif_start
bed['motif_end'] = motif_end


# should correspond to deletions
# bed[bed.motif_end-bed.motif_start>17].head()

bed[['motif_chrom', 'motif_start', 'motif_end', 'name', 'score', 'strand']].to_csv(sys.stdout, sep='\t', header=False, index=False)
