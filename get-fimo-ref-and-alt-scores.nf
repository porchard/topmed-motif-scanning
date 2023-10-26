#!/usr/bin/env nextflow

nextflow.enable.dsl=2

FIMO_REF_GLOB = params.fimo_ref_glob
FIMO_ALT_GLOB = params.fimo_alt_glob
VCF = params.vcf


process get_scores {

    publishDir "${params.results}"
    cache 'lenient'
    container 'library://porchard/default/general:20220107'
    memory '30 GB'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(motif), path('fimo_ref.txt'), path('fimo_alt.txt'), path(vcf)

    output:
    path("${motif}.txt")

    """
    get-fimo-ref-and-alt-scores.py fimo_ref.txt fimo_alt.txt $vcf > ${motif}.txt
    """

}


workflow {
    fimo_ref = Channel.fromPath(FIMO_REF_GLOB).map({it -> [it.getName().replace('.fimo.txt', ''), it]})
    fimo_alt = Channel.fromPath(FIMO_ALT_GLOB).map({it -> [it.getName().replace('.fimo.txt', ''), it]})
    vcf = Channel.fromPath(VCF)
    fimo_ref.combine(fimo_alt, by: 0).combine(vcf) | get_scores
}
