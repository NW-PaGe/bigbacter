/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { TIMESTAMP              } from '../modules/local/timestamp/main'
include { PREPARE                } from '../subworkflows/local/prepare/main'
include { CLASSIFY               } from '../subworkflows/local/classify/main'
include { VARIANTS               } from '../subworkflows/local/variants/main'
include { REPORT                 } from '../subworkflows/local/report/main'
include { PUSH                   } from '../subworkflows/local/push/main'

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bigbacter_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BIGBACTER {

    take:
    ch_samplesheet // channel: [ val(meta), path(reads)?, path(assembly)?, val(sra)?, val(genbank)? ]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /* 
    =============================================================================================================================
        TIMESTAMP GENERATION
        - Ensures consistent timestamps when resuming workflows
    =============================================================================================================================
    */

    ch_timestamp = TIMESTAMP()

    ch_samplesheet = ch_samplesheet
        .combine(ch_timestamp)
        .map { rec, timestamp ->
            rec.meta.timestamp = timestamp
            return rec
        }

    /* 
    =============================================================================================================================
        DATA PREPARATION
        - Fetch data from NCBI (SRA/GenBank)
        - Quality control and trimming
        - De novo assembly
    =============================================================================================================================
    */

    PREPARE(ch_samplesheet)
    ch_versions      = ch_versions.mix(PREPARE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE.out.fastqc_zip.collect { it[1] })

    /* 
    =============================================================================================================================
        TAXONOMIC CLASSIFICATION AND CLUSTERING
        - Species prediction (GAMBIT)
        - Genomic clustering (FLOC)
    =============================================================================================================================
    */

    CLASSIFY(PREPARE.out.manifest)
    ch_versions = ch_versions.mix(CLASSIFY.out.versions)

    /* 
    =============================================================================================================================
        VARIANT CALLING AND PHYLOGENETICS
        - Reference selection
        - Per-sample SNP calling
        - Core genome alignment
        - Recombination masking (optional)
        - Phylogenetic reconstruction
    =============================================================================================================================
    */

    VARIANTS(
        CLASSIFY.out.manifest,
        PREPARE.out.fastp_json
    )
    ch_versions = ch_versions.mix(VARIANTS.out.versions)

    /* 
    =============================================================================================================================
        REPORT GENERATION
        - Cluster-level summaries
        - Microreact visualizations
    =============================================================================================================================
    */

    REPORT(
        CLASSIFY.out.db_dist,
        VARIANTS.out.snp_summary,
        VARIANTS.out.snp_dist,
        VARIANTS.out.snp_tree,
        VARIANTS.out.ref_data
    )
    ch_versions = ch_versions.mix(REPORT.out.versions)

    /* 
    =============================================================================================================================
        DATABASE UPDATE (OPTIONAL)
        - Push results to database directory structure
    =============================================================================================================================
    */

    if (params.push) {
        PUSH(
            CLASSIFY.out.manifest,
            CLASSIFY.out.db,
            VARIANTS.out.ref_data,
            VARIANTS.out.snp_files, 
            REPORT.out.report
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'bigbacter_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
        ch_timestamp
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: [ path(multiqc_report.html) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/