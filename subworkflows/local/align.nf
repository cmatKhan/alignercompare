
include { HISAT2_EXTRACTSPLICESITES } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD              } from '../../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN              } from '../../modules/nf-core/hisat2/align/main'
include { STAR_GENOMEGENERATE       } from '../../modules/local/starlong/genomegenerate/main'
include { MINIMAP2_ALIGN            } from '../../modules/nf-core/minimap2/align/main'
include { BAM_SORT_STATS_SAMTOOLS   } from '../nf-core/bam_sort_stats_samtools/main'
include { FASTQ_ALIGN_STAR          } from '../nf-core/fastq_align_star/main'

workflow ALIGN {
    take:
    fastq // [ meta, file(fastq) ]
    genome // file: /path/to/genome.fasta
    gtf

    main:

    ch_versions = Channel.empty()
    ch_unsorted_bam = Channel.empty()

    MINIMAP2_ALIGN(
        fastq,
        genome,
        true,
        '',
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_unsorted_bam = ch_unsorted_bam.mix(
        MINIMAP2_ALIGN.out.bam.map{
            meta, bam -> [augment_meta(meta, 'minimap2'), bam]})


    if(!params.star_index){

        STAR_GENOMEGENERATE (
            genome,
            gtf
        )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        ch_star_index = STAR_GENOMEGENERATE.out.index

    } else{

        ch_star_index = Channel.fromPath(params.star_index).collect()

    }

    FASTQ_ALIGN_STAR(
        fastq,                         // channel: [ val(meta), [ reads ] ]
        ch_star_index,                 // channel: /path/to/star/index/
        gtf,                           // channel: /path/to/genome.gtf
        params.star_ignore_sjdbgtf,    // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
        params.seq_platform,           // string : sequencing platform
        params.seq_center,             // string : sequencing center
        genome                         // channel: /path/to/reference.fasta

    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)
    ch_unsorted_bam = ch_unsorted_bam.mix(
        FASTQ_ALIGN_STAR.out.bam.map{
            meta, bam -> [augment_meta(meta,'star'),bam]})

    if(!params.hisat2_splicesites){

        HISAT2_EXTRACTSPLICESITES (
            gtf
        )
        ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        ch_hisat2_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt

    } else{

        ch_hisat2_splicesites = Channel.fromPath(params.hisat2_splicesites).collect()

    }

    if(!params.hisat2_index){

        HISAT2_BUILD (
            genome,
            gtf,
            ch_hisat2_splicesites
        )
        ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        ch_hisat2_index = HISAT2_BUILD.out.index

    } else{

        ch_hisat2_index = Channel.fromPath(params.hisat2_index).collect()

    }

    HISAT2_ALIGN(
        fastq,
        ch_hisat2_index,
        ch_hisat2_splicesites
    )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
    ch_unsorted_bam = ch_unsorted_bam.mix(
        HISAT2_ALIGN.out.bam.map{
            meta, bam -> [augment_meta(meta, 'hisat2'), bam]})

    BAM_SORT_STATS_SAMTOOLS (
        ch_unsorted_bam,
        genome
    )

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]

    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions // channel: [ versions.yml ]
}

def augment_meta(Map meta, aligner) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.id = new_meta.id + "_" + aligner

    new_meta["aligner"] = aligner

    return new_meta
}
