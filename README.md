EpiDiverse-WGBS Pipeline
========================

**EpiDiverse/wgbs** is a bioinformatics analysis pipeline for aligning whole genome bisulfite sequencing data from non-model plant species.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [cutadapt](https://github.com/marcelm/cutadapt/)), aligns the reads ([erne-bs5](http://erne.sourceforge.net/) or [segemehl](https://www.bioinf.uni-leipzig.de/Software/segemehl/)), and performs extensive quality-control on the results using custom scripts and [Picard MarkDuplicates](https://broadinstitute.github.io/picard/). Methylation calling and mbias correction is performed with [Methyldackel](https://github.com/dpryan79/MethylDackel).

See the [output documentation](https://github.com/EpiDiverse/wgbs/wiki/Pipeline-Output) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run epidiverse/wgbs -profile test,<docker|singularity|conda>
```

iv. Start running your own analysis!

```bash
nextflow run epidiverse/wgbs -profile <docker|singularity|conda> --input /path/to/reads/dir --reference /path/to/ref/dir
```

See [usage docs](https://github.com/EpiDiverse/wgbs/wiki/Pipeline-Usage) for all of the available options when running the pipeline.

### Wiki Documentation

The epidiverse/wgbs pipeline comes with documentation about the pipeline, [found in the Wiki](https://github.com/EpiDiverse/wgbs/wiki):

1. [Installation](https://github.com/EpiDiverse/wgbs/wiki/Installation)
2. Pipeline configuration
    * [Local installation](https://github.com/EpiDiverse/wgbs/wiki/Installation#2-install-the-pipeline)
    * [Adding your own system config](https://github.com/EpiDiverse/wgbs/wiki/Installation#3-pipeline-configuration)
    * [EpiDiverse infrastructure](https://github.com/EpiDiverse/wgbs/wiki/Installation#appendices)
3. [Running the pipeline](https://github.com/EpiDiverse/wgbs/wiki/Pipeline-Usage)
4. [Output and how to interpret the results](https://github.com/EpiDiverse/wgbs/wiki/Pipeline-Output)
5. [Troubleshooting](https://github.com/EpiDiverse/wgbs/wiki/Troubleshooting)

### Credits

These scripts were originally written for use by the [EpiDiverse International Training Network](https://epidiverse.eu/), by Adam Nunn ([@bio15anu](https://github.com/bio15anu)) and Nilay Can ([@nilaycan](https://github.com/nilaycan)).

This project has received funding from the European Union’s Horizon 2020 research and innovation
programme under the Marie Skłodowska-Curie grant agreement No 764965

## Citation

If you use epidiverse/wgbs for your analysis, please cite it using the following doi: <placeholder>