/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAlignercompare.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input,
                           params.multiqc_config,
                           params.fasta,
                           params.gtf,
                           params.star_index,
                           params.bwamem2_index ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK            } from '../subworkflows/local/input_check'
include { GTF2BED                } from '../modules/local/gtf2bed'
include { ALIGN                  } from '../subworkflows/local/align'
include { MULTIQC_CUSTOM_BIOTYPE } from '../modules/local/multiqc_custom_biotype'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { BAM_RSEQC                   } from '../subworkflows/nf-core/bam_rseqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ALIGNERCOMPARE {

    ch_versions = Channel.empty()

    ch_fasta = Channel.fromPath(params.fasta).collect()

    ch_gtf = Channel.fromPath(params.gtf).collect()

    ch_biotypes_header_multiqc   = file("$projectDir/assets/biotypes_header.txt", checkIfExists: true)

    rseqc_modules = params.rseqc_modules ?
        params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } :
        []

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    GTF2BED (
        ch_gtf
    )
    ch_version = ch_versions.mix(GTF2BED.out.versions)

    //
    // SUBWORKFLOW: Align with a bunch of different aligners
    //
    ALIGN(
        INPUT_CHECK.out.reads,
        ch_fasta,
        ch_gtf
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    // set up downstream channels
    ALIGN.out.bam
        .join(ALIGN.out.bai)
        .set{ ch_bam_bai }

    GTF2BED.out.bed.view()

    ALIGN.out.bam
        .combine(ch_gtf)
        .set{ ch_featurecounts }

    //
    // SUBWORKFLOW: Run selected RSeQC modules
    //
    BAM_RSEQC(
        ch_bam_bai,
        GTF2BED.out.bed,
        rseqc_modules
    )
    ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)
    ch_bamstat_multiqc            = BAM_RSEQC.out.bamstat_txt
    ch_inferexperiment_multiqc    = BAM_RSEQC.out.inferexperiment_txt
    ch_innerdistance_multiqc      = BAM_RSEQC.out.innerdistance_freq
    ch_junctionannotation_multiqc = BAM_RSEQC.out.junctionannotation_log
    ch_junctionsaturation_multiqc = BAM_RSEQC.out.junctionsaturation_rscript
    ch_readdistribution_multiqc   = BAM_RSEQC.out.readdistribution_txt
    ch_readduplication_multiqc    = BAM_RSEQC.out.readduplication_pos_xls
    ch_tin_multiqc                = BAM_RSEQC.out.tin_txt

    //
    // SUBWORKFLOW: Run subread featurecounts to quantify coverage of various
    //              annotated features
    //
    SUBREAD_FEATURECOUNTS(
        ch_featurecounts
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    MULTIQC_CUSTOM_BIOTYPE (
        SUBREAD_FEATURECOUNTS.out.counts,
        ch_biotypes_header_multiqc
    )
    ch_featurecounts_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv
    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAlignercompare.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowAlignercompare.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    // SAMTOOLS metrics
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.idxstats.collect{it[1]}.ifEmpty([]))

    // featurecounts
    ch_multiqc_files = ch_multiqc_files.mix(ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]))

    // RSEQC metrics
    ch_multiqc_files = ch_multiqc_files.mix(ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(ch_tin_multiqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
