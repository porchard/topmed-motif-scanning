#!/usr/bin/env python
# coding: utf-8

import gzip
import logging
from Bio import SeqIO
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', required=True, help='Must not be zipped')
parser.add_argument('--vcf', required=True, help='Must be gzipped')
parser.add_argument('--flank', required=True, type=int, default=40, help='Default = 40')
args = parser.parse_args()

#fasta = '/lab/work/porchard/snp-aware-pwm-scanning/data/fasta/hg19.fa' # must not be gzipped
#VCF = '/lab/work/porchard/snp-aware-pwm-scanning/data/fasta/1000G.vcf.gz' # must be gzipped
#FLANK = 29
FLANK = args.flank

seqs = dict()

with open(args.fasta, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        logging.info('Loaded sequence {}'.format(record.id))
        seqs[record.id] = record

logging.info('Finished loading sequences\n')


# now make the new fasta
variant_count = 0
with gzip.open(args.vcf, 'rt') as vcf:
    for line in vcf:
        if line.startswith('#'):
            continue
        variant_count += 1
        if variant_count % 1000 == 0:
            logging.info(f'Processed {variant_count} variants')
        CHROM, POS, ID, REF, ALT = line.rstrip().split('\t')[:5]
        CHROM = f'chr{CHROM}' if CHROM not in seqs and f'chr{CHROM}' in seqs else CHROM
        POS = int(POS)
        if not CHROM in seqs:
            logging.warning(f'Skipping variant {CHROM} {POS} {ID} {REF} {ALT}; {CHROM} not in fasta')
            continue
        
        PYTHON_POS = POS - 1 # pos is indexed from 1; python indexes from 0

        # check that the ref matches as it should
        assert(REF == str(seqs[CHROM][PYTHON_POS:(PYTHON_POS+len(REF))].seq).upper())

        # get the flanking sequence
        LEFT_FLANK_START = PYTHON_POS - FLANK
        LEFT_FLANK_END = PYTHON_POS
        RIGHT_FLANK_START = PYTHON_POS + len(REF)
        RIGHT_FLANK_END = RIGHT_FLANK_START + FLANK
        LEFT_FLANK = str(seqs[CHROM][LEFT_FLANK_START:LEFT_FLANK_END].seq)
        RIGHT_FLANK = str(seqs[CHROM][RIGHT_FLANK_START:RIGHT_FLANK_END].seq)

        NAME = f'>{CHROM}_{FLANK}A@{POS}{REF}/{ALT}:{LEFT_FLANK_START+1}-{RIGHT_FLANK_END}'
        ALT_SEQUENCE = LEFT_FLANK + ALT + RIGHT_FLANK
        print(NAME)
        print(ALT_SEQUENCE)

