# snATACutils

## The utilities for single nucleus iATAC-seq analysis.
* Please find useful scirpts and codes in folder `bin`.
* We dissect the analysis into multiple steps, detailed description can be find in corresponding folders.
* This package is still under active development.

## Table of Contents
#### - 00.data processing
> This directory contains several necessary steps, including
> reads mapping, duplicates removal, matrix calculation and generation,
> quality control, doublets removel, etc.

#### - 01.cell_clustering
> This directory contains the basic and advanced steps for cell clustering.
> for small datasets, basic clustering strategy should work.
> for large/multiple datasets, more steps need to be included.

#### - 02.cluster_analysis
> We shared additional scripts for more detailed analysis, including
> cell clustering refinement,
> integration analysis with other modalities (i.e. snRNA-seq);

#### - 03.peak_calling
> To correct potential bias due to different depth/number of cells, we optimized 
> peak calling pipeline for snATAC-seq

#### - 04.peak_analysis
> Clustering of shared and cell-type specific candidate cis-regulatory elements (cCREs)
> Identification of cell-type/regional specific cCREs

### - others
> We shared few useful scripts in folder bin/:
> - phastCons score calculation in genomic regions;
> - script for generating IGV session file, which can be directly load to genome browser
> - loomR to seurat object convertor
> - loomR to matrix convertor
> - script for calculate silhouette score for each cell clustering

#### - flagship2020:
> Supplementary tables for putative enhancers identified in [flagship paper](https://www.nature.com/articles/s41586-021-03950-0)

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

