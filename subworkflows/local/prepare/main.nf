//
// Prepare reads and assemblies from multiple sources
//

include { FASTERQDUMP   } from '../../../modules/local/fasterq_dump/main'
include { SEQTK_SAMPLE  } from '../../../modules/nf-core/seqtk/sample/main'
include { FASTQC        } from '../../../modules/nf-core/fastqc/main'
include { FASTP         } from '../../../modules/nf-core/fastp/main'

include { NCBI_DATASETS } from '../../../modules/local/ncbi_datasets/main'
include { SHOVILL       } from '../../../modules/local/shovill/main'
include { STD_ASSEMBLY  } from '../../../modules/local/standardize/assembly/main'

workflow PREPARE {

    take:
    ch_samplesheet // channel: [ val(meta), path(reads)?, path(assembly)?, val(sra)?, val(genbank)? ]

    main:
    ch_versions = Channel.empty()

    /*
    ==================================================================================================
        USER-PROVIDED READS & ASSEMBLIES
    ==================================================================================================
    */

    ch_reads = ch_samplesheet
        .filter { it.reads }
        .map { [it.meta, it.reads] }

    ch_assemblies = ch_samplesheet
        .filter { it.assembly }
        .map { [it.meta, it.assembly] }

    /*
    ==================================================================================================
        NCBI DATA RETRIEVAL
    ==================================================================================================
    */

    // Fetch reads from SRA
    FASTERQDUMP(
        ch_samplesheet
            .filter { it.sra }
            .map { [it.meta, it.sra] }
    )
    ch_versions = ch_versions.mix(FASTERQDUMP.out.versions.first())
    ch_reads    = ch_reads.concat(FASTERQDUMP.out.reads)

    // Fetch assemblies from GenBank
    NCBI_DATASETS(
        ch_samplesheet
            .filter { it.genbank }
            .map { [it.meta, it.genbank] }
    )
    ch_versions   = ch_versions.mix(NCBI_DATASETS.out.versions.first())
    ch_assemblies = ch_assemblies.concat(NCBI_DATASETS.out.assembly)

    /*
    ==================================================================================================
        READ QUALITY CONTROL
    ==================================================================================================
    */

    // Drop any samples that still do not have reads
    ch_reads = ch_reads.filter{ meta, reads -> reads[0] }

    // Downsample high-coverage samples
    if (params.max_reads) {
        ch_reads
            .map { meta, reads ->
                def read_count
                try {
                    read_count = reads[0].countFastq() * 2
                } catch (Exception e) {
                    log.warn "${meta.id}: Failed to count reads - will attempt downsampling anyway!"
                    read_count = params.max_reads
                }
                [meta: meta, reads: reads, n: read_count]
            }
            .branch {
                ok:   it.n <= params.max_reads
                high: it.n > params.max_reads
            }
            .set { ch_reads_branched }

        // Downsample R1 and R2 separately
        SEQTK_SAMPLE(
            ch_reads_branched.high
                .map { [it.meta, it.reads[0], params.max_reads] }
                .concat(
                    ch_reads_branched.high
                        .map { [it.meta, it.reads[1], params.max_reads] }
                )
        )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

        // Recombine downsampled R1/R2 pairs
        ch_reads = SEQTK_SAMPLE
            .out
            .reads
            .groupTuple(by: 0)
            .concat(ch_reads_branched.ok.map { [it.meta, it.reads] })
    }

    // Run FastQC for quality assessment
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Trim and filter reads with fastp
    FASTP(
        ch_reads.map { meta, reads -> [meta, reads, []] },
        false, // save_trimmed_fail
        false, // save_merged
        false  // save_split
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_reads    = FASTP.out.reads

    /*
    ==================================================================================================
        DE NOVO ASSEMBLY
    ==================================================================================================
    */

    // Assemble reads for samples without pre-existing assemblies
    SHOVILL(
        ch_samplesheet
            .filter { !it.assembly && !it.genbank }
            .map { [it.meta.id] }
            .join(
                ch_reads.map { meta, reads -> [meta.id, meta, reads] },
                by: 0
            )
            .map { id, meta, reads -> [meta, reads] }
    )
    ch_versions   = ch_versions.mix(SHOVILL.out.versions.first())
    ch_assemblies = ch_assemblies.concat(SHOVILL.out.contigs)

    /*
    ==================================================================================================
        STANDARDIZE ASSEMBLY FORMAT
    ==================================================================================================
    */

    ch_assemblies_std = STD_ASSEMBLY(
            ch_assemblies
                .filter { meta, assembly -> 
                    assembly.name != "${meta.id}.fa.gz" 
                }
        )
        .concat(
            ch_assemblies.filter { meta, assembly -> 
                assembly.name == "${meta.id}.fa.gz" 
            }
        )

    /*
    ==================================================================================================
        COMBINE READS AND ASSEMBLIES
    ==================================================================================================
    */

    ch_manifest = ch_assemblies_std
        .map { meta, assembly -> [meta.id, meta, assembly] }
        .join(
            ch_reads.map { meta, reads -> [meta.id, reads] },
            by: 0,
            remainder: true
        )
        .map { id, meta, assembly, reads ->
            [meta: meta, reads: reads ? reads : [[],[]], assembly: assembly]
        }

    emit:
    manifest   = ch_manifest         // channel: [ val(meta), path(reads), path(assembly) ]
    fastqc_zip = FASTQC.out.zip      // channel: [ val(meta), path(zip) ]
    fastp_json = FASTP.out.json      // channel: [ val(meta), path(json) ]
    versions   = ch_versions         // channel: [ path(versions.yml) ]
}