#!/usr/bin/env nextflow

nextflow.enable.dsl=2

REF_FASTA = params.fasta
MOTIF_GLOB = params.meme_glob // must be MEME format, with basename {motif_name}.meme
VCF = params.vcf

// generate a bed file for reference
// and then generate a file for 

process make_alt_fasta {

    publishDir "${params.results}/alt-fasta"
    container 'library://porchard/default/general:20220107'
    memory '20 GB'

    input:
    path('ref.fasta')
    path(vcf)

    output:
    path('alt.fasta')

    """
    make-alt-fasta.py --fasta ref.fasta --vcf $vcf --flank 40 > alt.fasta
    """

}


process make_background {

    publishDir "${params.results}/background"
    container 'docker.io/memesuite/memesuite:5.5.3'
    memory '20 GB'

    input:
    path(fasta)

    output:
    path("background.txt")

    """
    fasta-get-markov $fasta > background.txt
    """

}


process fimo {

    publishDir "${params.results}/scan/${ref_or_alt}"
    maxForks 500
    container 'docker.io/memesuite/memesuite:5.5.3'
    memory '10 GB'
    time '6h'

    input:
    tuple val(motif), path(meme), val(ref_or_alt), path(fasta), path(bgfile)

    output:
    tuple val(motif), val(ref_or_alt), path("${motif}.fimo.txt")

    """
    fimo --text --bgfile $bgfile $meme $fasta > ${motif}.fimo.txt
    """

}


process make_ref_bed6 {

    maxForks 100
    publishDir "${params.results}/bed/ref"
    errorStrategy 'ignore'
    container 'library://porchard/default/general:20220107'
    memory '10 GB'

    input:
    tuple val(motif), val(ref_or_alt), path(fimo)

    output:
    tuple val(motif), val(ref_or_alt), path("${motif}.bed")

    when:
    ref_or_alt == 'ref'

    """
    fimo-to-bed6.py $fimo > ${motif}.bed
    """

}


process make_alt_bed6 {

    maxForks 100
    publishDir "${params.results}/bed/alt"
    errorStrategy 'ignore'
    container 'library://porchard/default/general:20220107'
    memory '10 GB'

    input:
    tuple val(motif), val(ref_or_alt), path(fimo)

    output:
    tuple val(motif), val(ref_or_alt), path("${motif}.bed")

    when:
    ref_or_alt == 'alt'

    """
    fimo-to-bed6.py $fimo > tmp.relative.bed
    update-alt-bed.py tmp.relative.bed > ${motif}.bed
    """

}


workflow {
    ref_fasta = Channel.fromPath(REF_FASTA)
    vcf = Channel.fromPath(VCF)
    motifs = Channel.fromPath(MOTIF_GLOB).map({it -> [it.getName().replaceAll('.meme', ''), it]})

    alt_fasta = make_alt_fasta(ref_fasta, vcf)
    fastas = ref_fasta.map({it -> ['ref', it]}).mix(alt_fasta.map({it -> ['alt', it]}))

    background = make_background(ref_fasta)
    fimo(motifs.combine(fastas).combine(background)) | (make_ref_bed6 & make_alt_bed6)

}