# scMPRA

## Introduction

This repository contains code, scripts, reasonably sized data associated with scMPRA: Biorxiv link. 

The code is sufficient to generate processed single-cell core promoter library expression data underlying the reproducibility plot, DE analysis.

The scripts are sufficient for reproducing the scatterplot, UMAP, and heatmaps. 

## Prerequisite 
All dependencies are described in conda environment mascot.yml. 

## Data Processing 

### Parse read 
A stand alone Go executable program (build for linux, amd64 architecture) is build for fast extraction of the barcodes. Briefly, the first 28 bps of Read1 are directly parsed as Cell BC and UMI. Sequences around the cBC and rBC were aligned to Read2 against sequencing and PCR errors, and cBC and rBC with correct length are kept. The program returns a tsv with columns of cellBC, UMI, cBC, rBC, count. 
First, make the program executable:

```
chmod -x ./scMPRA_parsing
```

```
./scMPRA_parsing Read1.fastq.gz Read2.fastq.gz 
```

The code should finish within an hour using a single GPU with 8 Gb of memory. 

If building for different architecture is desired,the source code for the program can be found on (https://github.com/szhao045/scMPRA_parsing). 

### Processing the transcriptome using cell ranger

To run Cell Ranger 6.0.1, first download cell ranger follow the instruction from 10X genomics:
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

To process the transcriptome data associated with scMPRA data, run

```
bash ./run_cell_ranger.sh
```

### Cell barcode filtering 

We only look at the 