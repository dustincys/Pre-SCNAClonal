# Pre-SCNAClonal

## Introduction

Tumor samples often present heterogeneity, containing not only multiple subpopulations of cancerous cells defined by distinct somatic mutations, but also the normal, non-cancerous cells.
While using the whole-genome sequencing data of tumor sample to reconstruct the subclonal composition determined by somatic copy number alternation, the absolute copy number and its subclonal population are both required to be estimated.
However, the raw reads in WGS data present coupling GC bias in the same SCNV genome segmentation of paired tumor and normal samples, which largely affect the absolute copy number estimation conducted by existing SCNV based subclonal inferring tools.
We provide Pre-SCNVClonal, a comprehensive package for both automatically and visually correcting coupling GC bias.
Pre-SCNVClonal could be strung together with the SCNV based subclonal inferring tool as a pipeline or run individually as needed.


## Requirement

- Python2.7

- python package:
 - `scipy`
 - `numpy`
 - `pickle`
 - `matplotlib`

## Usage


