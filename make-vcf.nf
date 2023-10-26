#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VCF_GLOB = params.vcf_glob
CHROMS = (1..22)
chroms_ordered = CHROMS.collect({it -> "chr" + it}).sort() + ['chrX']

// https://www.kingrelatedness.com/manual.shtml
// Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING.

process remove_genotypes {

    cache 'lenient'
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(chrom), path(vcf)

    output:
    path("${chrom}.bcf")

    """
    bcftools view --drop-genotypes -o ${chrom}.bcf -Ou $vcf
    """

}


process merge_bcfs {

    memory '40 GB'
    time '168h'
    publishDir "${params.results}/bcfs-merged"
    cache 'lenient'
    container 'library://porchard/default/general:20220107'

    input:
    path(vcf)

    output:
    path('merged.vcf.gz')

    script:
    tmp = chroms_ordered.collect({x -> x + ".bcf"}).join(' ')

    """
    bcftools concat -Oz -o merged.vcf.gz $tmp
    """

}


workflow {
    vcf_in = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[1], it]})
    remove_genotypes(vcf_in).toSortedList() | merge_bcfs
}
