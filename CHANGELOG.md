# nf-core/bigbacter: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - 09/01/2026

Major release. Clustering, core-genome analysis, and assembly all changed, and the
surveillance database layout is not backward compatible with v1. 

### `Added`
- **Repo and Documentation Locations**
  - Main repository migrated from `DOH-JDJ0303/bigbacter-nf` to `nw-page/bigbacter`. Links to old repo should automatically redirect.
  - Documentation moved from GitHub wiki to [GitHub Pages](https://nw-page.github.io/bigbacter-docs/).
- **Genome assembly:**
  - BigBacter now creates genome assemblies using [Shovill](https://github.com/tseemann/shovill) if one is not supplied in the samplesheet.
    This assembly is used, as is - no QC is performed! As such, it is still recommended to generate assemblies using a dedicated bacterial
    assembly workflow, like the [CDC PHoeNIx pipeline](https://github.com/CDCgov/phoenix).
  - Assemblies can now be sourced from NCBI GenBank using the [`datasets` CLI tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
    using the `genbank` column in the samplesheet.
  - Assemblies can be used in place of raw reads for variant calling (uses the Snippy `--ctgs` option).
    This method should only be used when the raw reads are not available!
- **Taxonomic classification:**
  - Bacterial samples can now be classified to the species level using [GAMBIT](https://github.com/jlumpe/gambit).
    Like genome assembly, a dedicated workflow is still recommended for this step when possible.
- **Clustering:**
  - `PopPUNK` replaced with [`Floc`](https://github.com/NW-PaGe/floc)
  - Cluster database is built on the fly from input assemblies - pre-built databases no longer necessary!
  - Clusters are now mutable via editing of sample signature files `<db>/<taxa>/sig/`
  - Clustering can by by-passed by supplying a cluster name in the samplesheet using the `cluster` column.
- **Reference Selection:**
  - Quality control steps now applied at the time of reference selection,
    including the number of contigs, genome length, and how representative
    the assembly is compare to the full cluster at the time of selection.
- **Core Genome Analysis**
  - Core genome analysis now performed using [`Polycore`](https://github.com/NW-PaGe/polycore) - replaces `snp-sites` and `snp-dist`.
  -  Supports a soft-core threshold via `min_genome_fraction` (0.8) and `min_core_fraction` (0.9), retaining more phylogenetic
  signal than strict core.
- **Microreact Report**
  - Microreact reports have been updated and include several new features
    
### `Fixed`
- see [Added](#Added)

### `Dependencies`

- Nextflow `>=25.04.0` is now required.
- IQ-TREE 2.3.4 → IQ-TREE 3.1.3. `--polytomy` replaced with `-czb`.
- Gubbins 3.3.1 → 3.4.1
- NCBI Datasets 16.15.0 → 18.18.0
- Added: Floc 1.1, Polycore 1.0.0, Shovill 1.4.2, sourmash (via Floc)
- Removed: PopPUNK 2.6.5, snp-dists 0.8.2, snp-sites
- Base image `bigbacter-base:1.0.0` → `bigbacter:2.0`
- Built on nf-core template 3.5.1 with nf-schema 2.5.1
- Snippy unchanged at 4.6.0-SC2

### `Deprecated`

- **v1 databases cannot be read by v2.0.** Layout changed from `<taxa>/pp_db/` and
  `<taxa>/clusters/<cluster>/snippy/` to `<taxa>/sig/` and
  `<taxa>/clusters/<cluster>/{asm,var,aux}/`. There is no migration path; archived
  assemblies must be re-clustered, and cluster IDs will not correspond to v1.
- Removed the `prepare_db` workflow, all 22 per-species config profiles, and the bundled
  accession lists. Floc makes them unnecessary.
- Removed merged-cluster resolution, which was specific to PopPUNK.
- Removed parameters: `ncbi`, `read_qc`, `assembly_qc`, `max_lowcov`, `max_het`, `max_ml`,
  `max_static`, `db_info`, `resolve_merged`, `run_id`, `tracedir`, and the `max_memory` /
  `max_cpus` / `max_time` caps.
- Renamed parameters:
  - `min_genfrac` (85) → `min_genome_fraction` (0.8) — **percentage to fraction**; a
    carried-over value will validate but behave incorrectly
  - `min_contig_len` → `ref_min_contig_len`
  - `strong_link_cutoff` → `strong_link_threshold`
  - `inter_link_cutoff` → `inter_link_threshold`
  - `partition_threshold` → `partition_distance`
- `--outdir` is now required. `--db` now defaults to `bigbacter_db`.
