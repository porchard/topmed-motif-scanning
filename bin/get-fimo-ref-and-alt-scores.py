#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pybedtools as bt
import pandas as pd
import numpy as np
import glob
import re
import logging

# for now, use only SNPs, not indels

#logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
#bt.set_tempdir('/localscratch/porchard/tmp')

FIMO_REF, FIMO_ALT, VCF = sys.argv[1:]

# VCF = '/net/topmed11/working/porchard/variant-sensitive-motif-scanning/work/vcf/results/bcfs-merged/merged.vcf.gz'
# FIMO_REF = '/net/topmed11/working/porchard/variant-sensitive-motif-scanning/work/scan/results/scan/ref/IKZF1_2.fimo.txt'
# FIMO_ALT = '/net/topmed11/working/porchard/variant-sensitive-motif-scanning/work/scan/results/scan/alt/IKZF1_2.fimo.txt'

variants = pd.read_csv(VCF, header=None, names=['chrom', 'pos', 'id', 'ref', 'alt'], sep='\t', usecols=[0, 1, 2, 3, 4], comment='#')
fimo_ref = pd.read_csv(FIMO_REF, sep='\t')
fimo_alt = pd.read_csv(FIMO_ALT, sep='\t')


def parse_sequence_name(s):
    RE = '^(.*)_(\d+)A@(\d+)(.*)/(.*):(\d+)-(\d+)$'
    chrom, flanking, pos, ref, alt, start, end = re.match(RE, s).groups()
    return {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'absolute_start': start, 'absolute_end': end}


for i in ['chrom', 'pos', 'ref', 'alt', 'absolute_start', 'absolute_end']:
    fimo_alt[i] = fimo_alt.sequence_name.map(lambda x: parse_sequence_name(x)[i])

fimo_alt.start = fimo_alt.start - 1 + fimo_alt.absolute_start.astype(int)
fimo_alt['end'] = fimo_alt.stop - 1 + fimo_alt.absolute_start.astype(int)
# bed format (end not inclusive)
fimo_ref.stop = fimo_ref.stop + 1
fimo_alt.end = fimo_alt.end + 1

fimo_ref = fimo_ref[fimo_ref.sequence_name.isin([f'chr{i}' for i in range(1, 23)] + ['chrX'])]
fimo_alt = fimo_alt[fimo_alt.chrom.isin([f'chr{i}' for i in range(1, 23)] + ['chrX'])]

fimo_ref['motif_instance'] = fimo_ref.motif_id + ':' + fimo_ref.sequence_name + ':' + fimo_ref.start.astype(str) + ':' + fimo_ref.stop.astype(str) + ':' + fimo_ref.strand
fimo_alt['motif_instance'] = fimo_alt.motif_id + ':' + fimo_alt.chrom + ':' + fimo_alt.start.astype(str) + ':' + fimo_alt.end.astype(str) + ':' + fimo_alt.strand

# get all motif hits that overlap SNPs
all_motif_hits = pd.concat([fimo_ref[['sequence_name', 'start', 'stop', 'motif_instance']].rename(columns={'sequence_name': 'chrom', 'stop': 'end'}), fimo_alt[['chrom', 'start', 'end', 'motif_instance']]]).drop_duplicates()
motif_hits_bt = bt.BedTool().from_dataframe(all_motif_hits).sort()
variants_bt = bt.BedTool().from_dataframe(variants.assign(start = lambda df: df.pos - 1, end = lambda df: df.pos).loc[:,['chrom', 'start', 'end']]).sort()

intersect = motif_hits_bt.intersect(variants_bt, wa=True, wb=True).to_dataframe()
# for each overlapping motif, get the ref score and the alt score
intersect = intersect.drop(columns=['chrom', 'start', 'end']).rename(columns={'name': 'motif_instance', 'score': 'chrom', 'thickStart': 'pos'}).loc[:,['motif_instance', 'chrom', 'pos']]
intersect = intersect.merge(variants)
# add scores
ref_scores = dict(zip(fimo_ref.motif_instance, fimo_ref.score))
alt_scores = dict(zip(fimo_alt[['motif_instance', 'chrom', 'pos', 'ref', 'alt']].astype(str).apply(lambda x: ':'.join(x), axis=1),
                     fimo_alt['score']))
intersect['ref_score'] = intersect.motif_instance.map(lambda x: ref_scores[x] if x in ref_scores else np.nan)
intersect['alt_score'] = intersect[['motif_instance', 'chrom', 'pos', 'ref', 'alt']].astype(str).apply(lambda x: ':'.join(x), axis=1).map(lambda x: alt_scores[x] if x in alt_scores else np.nan)

intersect = intersect[(intersect.ref.str.len()==1) & (intersect.alt.str.len()==1)]
intersect.to_csv(sys.stdout, index=False, sep='\t')