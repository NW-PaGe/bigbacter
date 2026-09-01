//
// Push analysis results to database directory structure
//

include { FILES as TAXA_FILES    } from '../../../modules/local/push/files/main'
include { FILES as CLUSTER_FILES } from '../../../modules/local/push/files/main'


workflow PUSH {
    
    take:
    ch_manifest  // channel: [ val(meta), path(assembly), path(reads) ]
    ch_db        // channel: [ val(meta), path(signature_db) ]
    ch_ref_data  // channel: [ val(meta), path(ref_assembly), path(ref_meta) ]
    ch_snp_files // channel: [ val(meta), path(snp_tar) ]
    ch_report    // channel: [ val(meta), path(report_files) ] - used as sentinel to prevent premature pushing
    
    main:
    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        PREPARE TAXA-LEVEL FILES
    =============================================================================================================================
    */

    ch_taxa = ch_db
        .map { meta, db -> 
            [[taxa: meta.taxa], db] 
        }
        .combine(
            ch_report
                .map { meta, files -> 
                    [[taxa: meta.taxa], files] 
                }
                .groupTuple(),
            by: 0
        )
        .map { meta, db, report_files ->
            def all_files = [db, report_files]
                .flatten()
                .findAll { it } // Remove null values
            
            [meta, all_files]
        }

    /*
    =============================================================================================================================
        PREPARE CLUSTER-LEVEL FILES
    =============================================================================================================================
    */

    ch_cluster = ch_manifest
        .map { [Utils.cluster_meta(it.meta), it.assembly] }
        .groupTuple()
        .join(
            ch_ref_data
                .map { meta, ref_path, ref_meta -> 
                    [Utils.cluster_meta(meta), ref_path, ref_meta] 
                },
            remainder: true
        )
        .join(
            ch_snp_files
                .map { meta, snp_files -> 
                    [Utils.cluster_meta(meta), snp_files] 
                }
                .groupTuple()
        )
        .combine(
            ch_report
                .map { meta, files -> 
                    [Utils.cluster_meta(meta), files] 
                }
                .groupTuple(),
            by: 0
        )
        .map { meta, assemblies, ref_path, ref_meta, snp_files, report_files ->
            def all_files = [
                assemblies,
                ref_path,
                ref_meta,
                snp_files,
                report_files // sentinel - not actually pushed
            ]
                .flatten()
                .findAll { it } // Remove null values
            
            [meta, all_files]
        }

    /*
    =============================================================================================================================
        PUSH FILES TO DATABASE
    =============================================================================================================================
    */

    TAXA_FILES( ch_taxa )
    CLUSTER_FILES( ch_cluster )
    
    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}