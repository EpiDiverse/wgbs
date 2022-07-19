# EpiDiverse-WGBS Runtime and memory usage guidelines
This document describes the default CPUs, RAM, and Time allocation specified for each pipeline process in the default configuration of the pipeline. Configuration was optimised on a HPC cluster with 64 CPUs and 256 Gb RAM, using the collection of plant population datasets provided by EpiDiverse. All values can be adjusted to suit individual needs.

|process|CPUs|RAM / Gb|Time / h|Retries|[errorStrategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)|
|-------|----|--------|--------|-------|-----------------|
|[erne_bs5_indexing](https://github.com/EpiDiverse/wgbs/blob/master/lib/index.nf#L5)|2|2|8|3|finish|
|[segemehl_indexing](https://github.com/EpiDiverse/wgbs/blob/master/lib/index.nf#L37)|2|2|8|3|finish|
|[read_trimming](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L12)|2|2|8|3|finish|
|[read_merging](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L59)|2|0.25|2|3|finish|
|[fastqc](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L95)|2|0.25|2|3|ignore|
|[erne_bs5](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L126)|8|1|6|3|finish|
|[segemehl](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L179)|8|8|12|3|finish|
|[erne_bs5_processing](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L233)|4|14|8|3|finish|
|[segemehl_processing](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L335)|2|2|8|3|finish|
|[bam_merging](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L382)|2|2|8|3|finish|
|[bam_subsetting](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L416)|2|2|10|3|finish|
|[bam_statistics](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L468)|2|2|10|3|ignore|
|[bam_filtering](https://github.com/EpiDiverse/wgbs/blob/master/lib/wgbs.nf#L498)|2|2|8|3|finish|
|[bam_processing](https://github.com/EpiDiverse/wgbs/blob/master/lib/call.nf#L2)|2|2|8|3|finish|
|[Picard_MarkDuplicates](https://github.com/EpiDiverse/wgbs/blob/master/lib/call.nf#L29)|2|10|10|3|finish|
|[MethylDackel](https://github.com/EpiDiverse/wgbs/blob/master/lib/call.nf#L69)|2|2|1|3|ignore|
|[conversion_rate_estimation](https://github.com/EpiDiverse/wgbs/blob/master/lib/call.nf#L115)|2|2|1|3|ignore|