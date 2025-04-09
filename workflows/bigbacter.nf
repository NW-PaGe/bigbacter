/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBigbacter.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.db, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (! (params.input || params.ncbi)) { exit 1, 'Input not specified!' } else{ manifest = Channel.empty()}
if (params.db) { ch_db = file(params.db) } else { exit 1, 'BigBacter database not specified!' }
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
include { PREPARE_INPUT    } from '../subworkflows/local/prepare_inputs'
include { CLUSTER          } from '../subworkflows/local/cluster'
include { CORE             } from '../subworkflows/local/core'
include { ACCESSORY        } from '../subworkflows/local/accessory'
include { PUSH_FILES       } from '../subworkflows/local/push_files'

include { MRFIGS           } from '../modules/local/mrfigs'
include { SUMMARY_TABLE    } from '../modules/local/summary-tables'
include { COMBINED_SUMMARY } from '../modules/local/combined-summary'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TIMESTAMP                   } from '../modules/local/get-timestamp'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// File to save timestamp
timestamp_file = file(workflow.workDir).resolve("${workflow.sessionId}_bb-timestamp")
// Info required for completion email and summary
def multiqc_report = []

workflow BIGBACTER {

    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        GET EPOCH TIMESTAMP
    =============================================================================================================================
    */
    // MODULE: Get epoch timestamp
    TIMESTAMP()
        .set { ch_timestamp }

    // provide custom run ID that replaces the timestamp - set up this way to avoid being a value channel 
    params.run_id ? ch_timestamp.map{ timestamp -> params.run_id }.set{ ch_timestamp } : ch_timestamp
    // save timestamp to file for reference outside the workflow - is there a better way?
    ch_timestamp.subscribe{ it -> file(timestamp_file).text = it }

    /*
    =============================================================================================================================
        PREPARE INPUT
    =============================================================================================================================
    */
    PREPARE_INPUT(
        file(params.input),
        ch_timestamp
    )
    // Create manifest channels
    PREPARE_INPUT.out.manifest.set{ ch_man }
    PREPARE_INPUT.out.manifest_path.set{ ch_man_path }

    /*
    =============================================================================================================================
        ASSIGN CLUSTERS
        - clusters assigned by PopPUNK when not manually provided in the samplesheet
    =============================================================================================================================
    */
    // SUBWORKFLOW: Assign PopPUNK clusters
    CLUSTER(
        ch_man.filter{ it -> ! it.cluster }.map{ it -> [ it.sample, it.taxa, it.assembly ] },
        ch_timestamp
    )
    ch_versions = ch_versions.mix(CLUSTER.out.versions)
    // Add cluster info back and combine with manual clusters
    ch_man
        .filter{ it -> it.cluster }
        .concat( CLUSTER.out.clusters )
        .map{ it -> [ it.sample, it.taxa, it.cluster ] }
        .set{ ch_clusters }

    /*
    =============================================================================================================================
        CORE GENOME ANALYSIS
    =============================================================================================================================
    */
    //// SUBWORKFLOW: Core genome analysis
    CORE(
        ch_man.map{ it -> [ it.sample, it.taxa, it.assembly, it.fastq_1, it.fastq_2 ] }.join( ch_clusters, by: [ 0, 1 ] ),
        ch_man_path,
        ch_timestamp
    )
    ch_versions = ch_versions.mix(CORE.out.versions)

    /*
    =============================================================================================================================
        ACCESSORY GENOME ANALYSIS
    =============================================================================================================================
    */
    // SUBWORKFLOW: Accessory genome analysis
    ACCESSORY(
        CLUSTER.out.dist,
        CORE.out.tree,
        ch_man_path,
        ch_timestamp
    )
    ch_versions = ch_versions.mix(ACCESSORY.out.versions)

    /*
    =============================================================================================================================
        SUMMARIZE RESULTS
    =============================================================================================================================
    */
    // Consolidate results for summarylines
    CORE
        .out
        .dist
        .map{ taxa, cluster, source, dist -> [ taxa, cluster, dist ] }
        .groupTuple(by: [0,1])
        .join(CORE.out.stats, by: [0,1])
        .set{ch_core_summary}

    // MODULE: Make individual summary tables
    SUMMARY_TABLE(
        ch_core_summary.combine(ch_man_path),
        ch_timestamp
    )

    // MODULE: Combine summary tables
    COMBINED_SUMMARY(
       SUMMARY_TABLE.out.summary.map{taxa, cluster, summary -> [ summary ]}.collect(),
       ch_man_path,
       ch_timestamp
   )

    // Consolidate results for Microreact
    // Accessory genome
    CORE
        .out
        .meta
        .join( CORE.out.dist, by: [0,1,2] )
        .join( CORE.out.tree, by: [0,1,2] )
        .combine( ACCESSORY.out.dist.map{ taxa, cluster, source, dist -> [ taxa, cluster, dist ] }.concat( ch_man.filter{ it -> it.cluster }.map{ it -> [ it.taxa, it.cluster, [] ] } ), by: [0,1] )
        .combine( SUMMARY_TABLE.out.summary, by: [0,1] )
        .set{ ch_microreact }
    // MODULE: Create Microreact figures
    MRFIGS (
        ch_microreact,
        file("$projectDir/assets/microreact.json", checkIfExists: true),
        ch_timestamp
    )

    /*
    =============================================================================================================================
        PUSH FILES TO BIGBACTER DATABASE
    =============================================================================================================================
    */
    // Consolidate taxa-specific files
    CLUSTER
        .out
        .new_db
        .set{ ch_taxa_files }

    // Consolidate cluster-specific files
    ch_man
        .map{ it -> [ it.sample, it.taxa, it.assembly ] }
        .combine( ch_clusters, by: [0,1] )
        .map{ sample, taxa, assembly, cluster -> [ taxa, cluster, assembly ] }
        .groupTuple(by: [0,1])
        .set{ ch_cluster_assemblies }
    CORE
        .out
        .snp_files
        .map{ taxa, cluster, ref, new_snippy, old_snippy -> [ taxa, cluster, ref, new_snippy ] }
        .join(ch_cluster_assemblies, by: [0,1])
        .set{ ch_cluster_files }

    // if push is 'true'
    ch_wait = ch_cluster_files.concat(ch_taxa_files).collect().flatten().last() // force pipeline to wait till after all new database files have been created
    if(params.push){
        // SUBWORKFLOW: Push new BigBacter database
        PUSH_FILES(
            ch_cluster_files,
            ch_taxa_files
        )
        PUSH_FILES.out.push_files.last().set{ch_wait} // update the wait channel
    }

    /*
    =============================================================================================================================
        NEXTFLOW DEFAULTS
    =============================================================================================================================
    */
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        ch_timestamp
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBigbacter.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowBigbacter.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_timestamp
    )
    multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
    // Reminder to push files
    if ( ! params.push & !(workflow.commandLine =~ "PREPARE_DB") ) { 
        msg = """
              \033[1;33m------------------------------------------------------
              
              Next steps:
              1. Check your results
              2. Run the command below to push this run to your database
              
              \033[0;37m> "${workflow.commandLine.replaceAll(/ -resume/, '')} --push true -resume
              
              \033[1;33m------------------------------------------------------\033[0m
              """.stripIndent()
        println msg }
    // Get database info
    if ( params.db_info ) { dbInfo( params.db, params.outdir, timestamp_file.text ) }
    else { println "\033[1;36mTip: Run BigBacter with\033[0m `--db_info true` \033[1;36mto get a summary of your BigBacter database!\033[0m\n" }
    // nf-core email
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function for gathering BigBacter database info
def dbInfo ( db_path, outdir_path, timestamp ) {
    // Print message
    msg = """
          Gathering BigBacter database info. This can take a while if your database is large and/or in the cloud. \033[0;31mPress 'control + c' to skip.\033[0m
          You can skip this in the future by running BigBacter with '--db_info false'
          """.stripIndent()
    println msg
    // Define output summary file
    def db_info = file(outdir_path).resolve(timestamp.toString()).resolve(timestamp.toString()+"-db-info.csv")
    // Gather File Info
    //// Taxa
    def taxon    = file(db_path).exists() ? file(db_path).list().collect { [ taxon: it, taxon_path: file(db_path).resolve(it).exists() ? file(db_path).resolve(it) : null ] }.findAll{ it } : null
    //// Clusters
    def clusters = taxon ? taxon.collectMany { 
        it -> if( it.taxon_path.resolve('clusters').exists() ) { 
            it.taxon_path.resolve('clusters').list().collect { 
                cluster -> [ taxon: it.taxon, 
                             cluster: cluster, 
                             cluster_path: it.taxon_path.resolve('clusters').resolve(cluster) ] 
                }
            } else { [] }
        }.findAll{ it ? ( it.cluster_path ? true : false ) : false } : []
    //// PopPUNK Databases
    def pp_dbs  = taxon ? taxon.collectMany { 
        it -> if( it.taxon_path.resolve('pp_db').exists() ) {
            it.taxon_path.resolve('pp_db').list().collect { 
                pp_db ->  [ taxon: it.taxon, cluster: null, fileType: "PopPUNK_Database", fileName: pp_db, filePath: it.taxon_path.resolve('pp_db').resolve(pp_db).exists() ? it.taxon_path.resolve('pp_db').resolve(pp_db) : null ] 
                }
            } else { [] }
        } : []
    //// Snippy Files
    def snippy  = clusters ? clusters.collectMany { 
        it -> if( it.cluster_path.resolve('snippy').exists() ) {
            it.cluster_path.resolve('snippy').list().collect { 
                snippy_file -> [ taxon: it.taxon, cluster: it.cluster, fileType: "SNP_Files", fileName: snippy_file, filePath: it.cluster_path.resolve('snippy').resolve(snippy_file).exists() ? it.cluster_path.resolve('snippy').resolve(snippy_file) : null ] 
                }
            } else { [] } 
        } : []
    //// SNP References
    def refs  = clusters ? clusters.collectMany { 
        it -> if( it.cluster_path.resolve('ref').exists() ) {
            it.cluster_path.resolve('ref').list().collect { 
                ref -> [ taxon: it.taxon, cluster: it.cluster, fileType: "SNP_Reference", fileName: ref, filePath: it.cluster_path.resolve('ref').resolve(ref).exists() ? it.cluster_path.resolve('ref').resolve(ref) : null ] 
                }
            } else { [] } 
        } : []
    //// Assembly Files
    def assemblies  = clusters ? clusters.collectMany { 
        it -> if( it.cluster_path.resolve('assembly').exists() ) {
            it.cluster_path.resolve('assembly').list().collect { 
                assembly -> [ taxon: it.taxon, cluster: it.cluster, fileType: "Assembly_File", fileName: assembly, filePath: it.cluster_path.resolve('assembly').resolve(assembly).exists() ? it.cluster_path.resolve('assembly').resolve(assembly) : null ] 
                }
            } else { [] } 
        } : []
    //// Combine it all together
    def all_files = ( pp_dbs + snippy + refs + assemblies ).findAll { it }.collect { it + [ fileSize: it.filePath.size(), fileDate: new Date(it.filePath.lastModified()) ] }

    // Report Summary (if anything was collected)
    if( all_files ){
        // Save to file
        db_info.text = ( [ all_files[0].keySet().join(',') ] + all_files.collect{ it.values().join(',') } ).join('\n') + '\n'

        // Print to stdout
        def db_summary = all_files.groupBy{ it.fileType }.collectEntries{ group, list -> [ group, list.sum{ it.fileSize } ] }
        db_summary = ( db_summary + [ Total: db_summary.values().sum() ] ).entrySet().collect{ it -> "${it.key}: ${ it.value >= 1e9 ? MemoryUnit.of(it.value).toGiga() : MemoryUnit.of(it.value).toMega() } ${ it.value >= 1e9 ? 'GB' : 'MB' }" }
        ( [ '=' * 27, 'BigBacter Database Summary:', '=' * 27 ] + db_summary ).each{ println it }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
