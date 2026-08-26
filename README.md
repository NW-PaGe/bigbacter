[![](https://github.com/NW-PaGe/bigbacter-docs/raw/a0184ebe0bf6341508ca483f587a5ede6e6e0c13/assets/media/bigbacter_logo.png)](https://github.com/NW-PaGe/bigbacter-docs/blob/a0184ebe0bf6341508ca483f587a5ede6e6e0c13/assets/media/bigbacter_logo.png)

# BigBacter: Bacterial Genomic Surveillance

## Key Features

BigBacter is a bacterial genomic surveillance pipeline that can:

<div style="padding: 1em; margin: 1em 0;">

🧬 <strong>Iterative clustering</strong> - cluster assignments stay consistent across runs using a per-sample sourmash database that expands automatically with each new submission<br>
🧬 <strong>Soft-core phylogenomics</strong> - retains substantially more phylogenetic signal than strict-core approaches by tolerating a configurable level of missing data<br>
🧬 <strong>Automated reference selection</strong> - selects the most representative assembly per cluster using k-mer containment and assembly quality scoring; reuses the same reference on subsequent runs for consistent SNP distances<br>
🧬 <strong>Dual distance metrics</strong> - reports both core-genome SNP distances and whole-genome containment scores to capture both SNP-level and accessory genome variation<br>

</div>  

BigBacter accepts reads, assemblies, or **SRA and GenBank accessions**, and fills in whatever is missing — assembling, subsampling, and assigning species-level taxonomy along the way. Because clustering is k-mer based and stored per-sample, **no reference database has to be built before a new species can be run**, and a run can start with as little as a single isolate. BigBacter is designed to run downstream of a workflow with robust assembly QC, such as [PHoeNIx](https://github.com/CDCgov/phoenix), [TheiaProk](https://public-health-bacterial-genomics-theiagen.readthedocs.io/en/latest/theiaprok.html), or [Bactopia](https://github.com/bactopia/bactopia).

## Pipeline Overview

![](https://nw-page.github.io/bigbacter-docs/docs/v2.0/media/bigbacter-v2_0-nw-page.png)

## More Information:

See the [docs](https://nw-page.github.io/bigbacter-docs/docs/v2.0/) for more information.