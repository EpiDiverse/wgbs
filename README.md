EpiDiverse-WGBS Pipeline
========================

[<img width="200" align="right" src="docs/images/euflagbetter.jpg">](https://ec.europa.eu/programmes/horizon2020/en)
[<img width="200" align="right" src="docs/images/epidiverse-logo.jpg">](https://epidiverse.eu)
[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/epidiverse/wgbs.svg)](https://hub.docker.com/r/epidiverse/wgbs)
[![Release](https://img.shields.io/github/v/release/epidiverse/wgbs)]()
[![Publication](https://img.shields.io/badge/Published-bioRxiv-26af64.svg?colorB=26af64&style=popout)](https://www.biorxiv.org/content/10.1101/2020.08.28.271585v1)
[![Twitter](https://img.shields.io/twitter/follow/epidiverse?style=social)](https://twitter.com/intent/follow?screen_name=epidiverse)

---

EpiDiverse-WGBS Pipeline
========================

**EpiDiverse/wgbs** is a bioinformatics analysis pipeline for aligning whole genome bisulfite sequencing data from non-model plant species.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [cutadapt](https://github.com/marcelm/cutadapt/)), aligns the reads ([erne-bs5](http://erne.sourceforge.net/) or [segemehl](https://www.bioinf.uni-leipzig.de/Software/segemehl/)), and performs extensive quality-control on the results using custom scripts and [Picard MarkDuplicates](https://broadinstitute.github.io/picard/). Methylation calling and mbias correction is performed with [Methyldackel](https://github.com/dpryan79/MethylDackel).

> See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow NXF_VER=20.07.1 run epidiverse/wgbs -profile test,<docker|singularity|conda>
```

iv. Start running your own analysis!

```bash
nextflow NXF_VER=20.07.1 run epidiverse/wgbs -profile <docker|singularity|conda> \
--input /path/to/reads/directory --reference /path/to/reference.fasta
```

> See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.

### Wiki Documentation

The EpiDiverse/wgbs pipeline is part of the [EpiDiverse Toolkit](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/overview), a best practice suite of tools intended for the study of [Ecological Plant Epigenetics](https://app.gitbook.com/@epidiverse/s/project/). Links to general guidelines and pipeline-specific documentation can be found below:

1. [Installation](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation)
2. Pipeline configuration
    * [Local installation](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation#2-install-the-pipeline)
    * [Adding your own system config](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation#3-pipeline-configuration)
    * [EpiDiverse infrastructure](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation#appendices)
3. [Running the pipeline](docs/usage.md)
4. [Understanding the results](docs/output.md)
5. [Troubleshooting](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/troubleshooting)

### Credits

These scripts were originally written for use by the [EpiDiverse European Training Network](https://epidiverse.eu/), by Adam Nunn ([@bio15anu](https://github.com/bio15anu)) and Nilay Can ([@nilaycan](https://github.com/nilaycan)).

This project has received funding from the European Union’s Horizon 2020 research and innovation
programme under the Marie Skłodowska-Curie grant agreement No 764965

## Citation

If you use epidiverse/wgbs for your analysis, please cite it using the following doi: <placeholder>