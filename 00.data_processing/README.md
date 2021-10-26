# 00.data processing

This directory contains several necessary steps, including
reads mapping, duplicates removal, matrix calculation and generation,
quality control, doublets removel, etc.

Please find useful scirpts and codes in folder `bin`.

## Requirements
* python v3.6.8
* R v3.6.0
* Perl v5.26.2

## This pipeline is build on multiple softwares and tools:
* bwa: [link](http://bio-bwa.sourceforge.net)
* samtools: [link](https://github.com/samtools/samtools)
* bedtools: [link](https://bedtools.readthedocs.io/en/latest/)
* Snaptools: [link](https://github.com/r3fang/SnapTools)
* SnapATAC: [link](https://github.com/r3fang/SnapATAC)
* Snakemake: [link](https://snakemake.readthedocs.io/en/stable/)

## Reference and annotation:
* mouse genome: [mm10](https://www.gencodegenes.org/mouse/release_M16.html)
* mouse genome annotation: [gencode vM16](https://www.gencodegenes.org/mouse/release_M16.html)
* mouse blacklist: [ENCODE blacklist](https://github.com/Boyle-Lab/Blacklist)

