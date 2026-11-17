//
// Accessory genome analysis (optional, --accessory true)
// Annotates every assembly in a cluster with Bakta, re-using annotations already in the database,
// then builds a gene presence-absence matrix for the cluster with Panaroo and derives pairwise
// accessory distances from it with pan-dists.
//

include { BAKTA_DOWNLOAD } from '../../../modules/local/bakta/download/main'
include { BAKTA          } from '../../../modules/local/bakta/annotate/main'
include { PANAROO        } from '../../../modules/local/panaroo/main'
include { PAN_DISTS      } from '../../../modules/local/pan_dists/main'

workflow ACCESSORY {

    take:
    ch_manifest // channel: [ val(meta), path(reads), path(assembly) ]

    main:
    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        SAMPLES IN EACH CLUSTER (THIS RUN + DATABASE) AND THEIR EXISTING ANNOTATIONS
    =============================================================================================================================
    */

    ch_cluster_samples = ch_manifest
        .map { [Utils.cluster_meta(it.meta), [it.meta.id.toString(), it.assembly]] }
        .groupTuple()
        .map { meta, run_samples ->
            def cluster_path = file(params.db) / meta.taxa / 'clusters' / meta.cluster
            def asm_dir      = cluster_path / 'asm'
            def ann_dir      = cluster_path / 'ann'

            // annotations already in the database: <cluster>/ann/<sample>.gff.gz
            def db_ann = ann_dir.exists()
                ? ann_dir.listFiles()
                    .findAll { it.name.endsWith('.gff.gz') }
                    .collectEntries { [it.name.replaceAll(/\.gff\.gz$/, ''), it] }
                : [:]

            // samples from this run. They replace database samples of the same name, as for variants.
            def samples = run_samples.collect { id, assembly ->
                [id: id, assembly: assembly, ann: params.ann_overwrite ? null : db_ann[id]]
            }
            def run_ids = samples.collect { it.id }

            // samples from earlier runs: <cluster>/asm/<sample>.fa.gz
            if (asm_dir.exists()) {
                samples += asm_dir.listFiles()
                    .findAll { it.name.endsWith('.fa.gz') && it.name != 'ref.fa.gz' }
                    .collect { [id: it.name.replaceAll(/\.fa\.gz$/, ''), assembly: it] }
                    .findAll { !run_ids.contains(it.id) }
                    .collect { it + [ann: db_ann[it.id]] }
            }

            // sorted so task hashes are stable when resuming
            [meta, samples.sort { it.id }]
        }

    // [ val(cluster_meta), [ path(gff) ] ] - may be an empty list
    ch_ann_db = ch_cluster_samples
        .map { meta, samples -> [meta, samples.findAll { it.ann }.collect { it.ann }] }

    // [ val(meta), path(assembly) ] - one per sample without an annotation
    ch_to_annotate = ch_cluster_samples
        .flatMap { meta, samples ->
            samples.findAll { !it.ann }.collect { [meta + [id: it.id], it.assembly] }
        }

    /*
    =============================================================================================================================
        BAKTA DATABASE
        - --bakta_db if given, else <db>/bakta/db-<type>
        - Downloaded and saved to the BigBacter database the first time it is needed
    =============================================================================================================================
    */

    def bakta_db = params.bakta_db
        ? file(params.bakta_db, checkIfExists: true)
        : file(params.db) / 'bakta' / "db-${params.bakta_db_type}"

    if (params.bakta_db || (bakta_db / 'version.json').exists()) {
        ch_bakta_db = Channel.value(bakta_db)
    }
    else {
        // only runs when there is something to annotate
        BAKTA_DOWNLOAD(ch_to_annotate.first().map { params.bakta_db_type })
        ch_versions = ch_versions.mix(BAKTA_DOWNLOAD.out.versions)
        ch_bakta_db = BAKTA_DOWNLOAD.out.db
    }

    /*
    =============================================================================================================================
        ANNOTATION
    =============================================================================================================================
    */

    BAKTA(ch_to_annotate.combine(ch_bakta_db))
    ch_versions = ch_versions.mix(BAKTA.out.versions.first())

    // existing and new annotations of each cluster
    ch_ann = ch_ann_db
        .join(
            BAKTA.out.gff
                .map { meta, gff -> [Utils.cluster_meta(meta), gff] }
                .groupTuple(),
            remainder: true // clusters with nothing new to annotate
        )
        .map { meta, ann_db, ann_new -> [meta, ann_db, ann_new ?: []] }

    /*
    =============================================================================================================================
        GENE PRESENCE-ABSENCE
    =============================================================================================================================
    */

    PANAROO(
        ch_ann
            .map { meta, ann_db, ann_new ->
                def gffs = (ann_db + ann_new).sort { it.name }
                [meta + [n: gffs.size()], gffs]
            }
            .filter { meta, gffs -> meta.n > 1 }
    )
    ch_versions = ch_versions.mix(PANAROO.out.versions.first())

    /*
    =============================================================================================================================
        ACCESSORY DISTANCES
        - Accessory genes within N genes of one another are grouped into gain/loss events
        - Each pair is scored by how many events differ, rather than by how many genes differ
    =============================================================================================================================
    */

    PAN_DISTS(PANAROO.out.results)
    ch_versions = ch_versions.mix(PAN_DISTS.out.versions.first())

    emit:
    ann_new  = ch_ann.map { meta, ann_db, ann_new -> [meta, ann_new] } // channel: [ val(meta), [ path(gff) ] ] - saved by PUSH
    results  = PANAROO.out.results                                     // channel: [ val(meta), path(panaroo_files) ]
    rtab     = PANAROO.out.rtab                                        // channel: [ val(meta), path(gene_presence_absence.Rtab) ]
    csv      = PANAROO.out.csv                                         // channel: [ val(meta), path(gene_presence_absence.csv) ]
    acc_dist = PAN_DISTS.out.dist_events_all                           // channel: [ val(meta), path(dist_events_all.tsv) ]
    versions = ch_versions                                             // channel: [ path(versions.yml) ]
}
