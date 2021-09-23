# EpiDiverse-WGBS Runtime and memory usage guidelines
This document describes the default CPUs, RAM, and Time allocation specified for each pipeline process in the default configuration of the pipeline. Configuration was optimised on a HPC cluster with 64 CPUs and 256 Gb RAM, using the collection of plant population datasets provided by EpiDiverse. All values can be adjusted to suit individual needs.

|process|CPUs|RAM / Gb|Time / h|Retries|[errorStrategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)|
|-------|----|--------|--------|-------|-----------------|
|erne_bs5_indexing|2|2|8|3|finish|
|segemehl_indexing|2|2|8|3|finish|
|read_trimming|2|2|8|3|finish|
|read_merging|2|0.25|2|3|finish|
|fastqc|2|0.25|2|3|ignore|
|erne_bs5|8|1|6|3|finish|
|segemehl|8|8|12|3|finish|
|erne_bs5_processing|4|14|8|3|finish|
|segemehl_processing|2|2|8|3|finish|
|bam_merging|2|2|8|3|finish|
|bam_subsetting|2|2|10|3|finish|
|bam_statistics|2|2|10|3|ignore|
|bam_filtering|2|2|8|3|finish|
|bam_processing|2|2|8|3|finish|
|Picard_MarkDuplicates|2|10|10|3|finish|
|MethylDackel|2|2|1|3|ignore|
|conversion_rate_estimation|2|2|1|3|ignore|