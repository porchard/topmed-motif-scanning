# TOPMed variant-sensitive motif scanning


## Dependencies
* Singularity (v. >=3)
* NextFlow (v. >= 21.04.0)

## Setup
 
<!-- 1. Update nextflow.config as necessary for your compute environment.
2. Prep data: `make data` -->

## Running

1. Make VCF `make vcf`
2. Scan hg38: `make scan`
3. Wrangle motif ref/alt variant scores: `make fimo-ref-and-alt-scores`