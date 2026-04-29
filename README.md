# miR-29 ncRNA Targetome Analysis

This repository contains analysis scripts used for the dissertation:

**Identification of MicroRNA-29 Noncoding RNA Targetome in Keratinocytes**

## Overview

These scripts analyse RNA-seq data (GSE270877) to examine expression and correlation patterns of candidate lncRNAs identified from CLIP-seq analysis.

Candidate lncRNAs:
- H19
- SNHG7
- LINC00511

## Analysis steps

- Import raw count data
- DESeq2 normalization (VST)
- Extraction of lncRNA expression
- Pearson correlation analysis
- Visualization (bubble plots)

## Input data

GSE270877_Raw_Counts.txt
in the working directory.

## Requirements

R packages:
```r
DESeq2
ggplot2
dplyr
tidyr
tibble
purrr

## Notes
Correlation was calculated using Pearson correlation on VST-normalized expression values. Genes with low variability or missing values were excluded.
