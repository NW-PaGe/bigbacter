//
// Prepare the database files of a run.
// They are always saved to the run directory (<outdir>/<timestamp>) in a layout that mirrors the
// database, so a run can be pushed later, and are also pushed to the database with --push.
// The publish paths are set in conf/modules.config (TAXA_FILES, CLUSTER_FILES).
//

include { FILES as TAXA_FILES    } from '../../../modules/local/push/files/main'
include { FILES as CLUSTER_FILES } from '../../../modules/local/push/files/main'


workflow PUSH {
    
    take:
    ch_manifest  // channel: [ val(meta), path(assembly), path(reads) ]
    ch_db        // channel: [ val(meta), path(signature_db + config.json) ]
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

    // Samples in this run, per taxa. FLOC returns every signature in the taxa database, but only
    // the signatures of this run's samples are saved to the run directory.
    ch_ids = ch_manifest
        .map { [[taxa: it.meta.taxa], it.meta.id.toString()] }
        .groupTuple()

    ch_taxa = ch_db
        .map { meta, db -> 
            [[taxa: meta.taxa], meta.timestamp, db] 
        }
        .combine(
            ch_report
                .map { meta, files -> 
                    [[taxa: meta.taxa], files] 
                }
                .groupTuple(),
            by: 0
        )
        .combine(ch_ids, by: 0)
        .map { meta, timestamp, db, report_files, ids ->
            // sorted so the task hash is stable when resuming
            def run_ids = ids.collect().sort()

            // without --push nothing goes to the database, so only this run's signatures are staged
            def db_files = params.push
                ? db
                : db.findAll { it.name == 'config.json' || run_ids.contains(it.name - '.sig.gz') }

            def all_files = [db_files, report_files]
                .flatten()
                .findAll { it } // Remove null values
            
            // timestamp and ids are needed to save config.json and this run's signatures to the outdir
            [meta + [timestamp: timestamp, ids: run_ids], all_files]
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
        SAVE FILES TO THE RUN DIRECTORY AND (WITH --push) THE DATABASE
    =============================================================================================================================
    */

    TAXA_FILES( ch_taxa )
    CLUSTER_FILES( ch_cluster )
    
    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}