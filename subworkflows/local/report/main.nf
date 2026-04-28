//
// Generate cluster-level reports and visualizations
//

include { CLUSTER } from '../../../modules/local/report/cluster/main'
include { ALL     } from '../../../modules/local/report/all/main'

workflow REPORT {
    
    take:
    ch_db_dist     // channel: [ val(meta), path(db_distance_matrix) ]
    ch_snp_summary // channel: [ val(meta), path(snp_summary_csv) ]
    ch_snp_dist    // channel: [ val(meta), path(snp_distance_matrix) ]
    ch_snp_tree    // channel: [ val(meta), path(tree) ]
    ch_ref_data    // channel: [ val(meta), path(ref_assembly), path(ref_meta) ]
    
    main:
    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        COMBINE CLUSTER-LEVEL DATA FOR REPORTING
    =============================================================================================================================
    */

    ch_report_data = ch_snp_summary
        .join(ch_snp_dist)
        .join(ch_snp_tree, remainder: true)
        .map { meta, snp_summary, snp_dist, snp_tree -> 
            [
                Utils.cluster_meta(meta),
                snp_summary,
                snp_dist,
                snp_tree,
                meta.recomb_masked
            ]
        }
        .combine(
            ch_ref_data
                .map { meta, ref_path, ref_meta -> 
                    [Utils.cluster_meta(meta), ref_meta] 
                }
                .combine(ch_db_dist, by: 0),
            by: 0
        )
        .map { meta, snp_summary, snp_dist, snp_tree, recomb_masked, ref_meta, db_dist ->
            // Update metadata with recombination masking status
            def meta_updated = meta + [recomb_masked: recomb_masked]
            
            // Collect all available files (filter out null values)
            // Note: snp_dist is now a list [dist_wide, dist_long]
            def all_files = [
                snp_summary,
                snp_tree,
                ref_meta,
                db_dist
            ].findAll { it }
            
            // Add dist files (may be a list)
            if (snp_dist instanceof List) {
                all_files.addAll(snp_dist)
            } else if (snp_dist) {
                all_files.add(snp_dist)
            }

            [meta_updated, all_files]
        }

    /*
    =============================================================================================================================
        GENERATE CLUSTER REPORTS
    =============================================================================================================================
    */

    CLUSTER(
        ch_report_data,
        file(params.input),
        file(params.microreact_template)
    )
    ch_versions = ch_versions.mix(CLUSTER.out.versions.first())


    /*
    =============================================================================================================================
        GENERATE COMBINED REPORT
    =============================================================================================================================
    */

    ALL(
        CLUSTER.out.summary.map{ meta, summary ->
            def new_meta = [ timestamp: meta.timestamp ]
            [new_meta, summary]
        }
        .groupTuple(by: 0)
    )

    emit:
    report   = CLUSTER.out.report  // channel: [ val(meta), path(report_files) ]
    versions = ch_versions         // channel: [ path(versions.yml) ]
}