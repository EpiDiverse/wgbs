# EpiDiverse/wgbs
---
# Releases

---
# Prereleases

## v0.10.0dev
* Modularised CALL workflow to allow for repeat methylation calling
* Add scatter plot to check linearity of alignments vs duplicates

## v0.9.1dev
* Fix samtools merge issue with duplicated header entries
* Fix samtools sort issue with tlens correction

## v0.9.0dev - 2019.10.22
* Major changes to Nextflow DSL2
* Added capability to index reference genome in separate workflow
* Added capacity for plot-bamstats with gnuplot=5.2.6
* Added integration for conda, docker, singularity
* Added test profile to link with epidiverse/datasets
* Refactored how the config profiles are loaded