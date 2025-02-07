# TOPMed variant-sensitive motif scanning


## Dependencies
* Singularity (v. >=3)
* NextFlow (v. >= 21.04.0)

## Setup
 
1. Update nextflow.config as necessary for your compute environment.
2. Update the path at the top of the Makefile to reflect your current working directory.
3. Unpack the motifs: `make motifs`
4. Make a directory, `data/vcf`, and place the (gzipped) VCF into that directory (should be one VCF containing all variants; genotypes need not be included).
5. Make a directory, `data/fasta`, and place the TOPMed hg38 fasta file (`Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta`) into that directory.

After the above, an example `data` directory will look something like:

```
data
├── fasta
    └── Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
├── motifs
│   ├── AFP_1.meme
│   ├── AHR_1.meme
│   ├── AHR_2.meme
│   ├── ...
│   ├── ZSCAN26_1.meme
│   ├── ZSCAN4_2.meme
│   └── ZSCAN4_3.meme
├── motifs.tar.gz
└── vcf
    └── merged.vcf.gz
```

## Running

1. Make VCF `make vcf`
2. Scan hg38: `make scan`
3. Wrangle motif ref/alt variant scores: `make fimo-ref-and-alt-scores`.
