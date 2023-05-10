# Introduction

This repository contains all analysis code from the following paper:

Protein-intrinsic properties and epigenetic effects regulate pioneer-factor binding and function Tyler J Gibson, Melissa Harrison bioRxiv 2023.03.18.533281; doi: <https://doi.org/10.1101/2023.03.18.533281>

If you are interested in exploring the processed data from this paper, bigWig files and lists of peaks have been deposited in the Gene Expression Omnibus under accession number GSE227884.
Lists of class I, II, and III peaks for Zelda, Grainy head and Twist can be found as supplemental tables along with the published paper.

For anyone interested in reproducing the full analysis pipeline, follow the steps below.

# Setup

We used the workflow management system Snakemake to run all of the analysis for this manuscript.
To get started, follow the steps below to install the necessary software. You will only need to install Python3 and Snakemake itself. All other software dependencies for the pipeline are included in a Docker container that is available on Dockerhub. The Snakemake pipeline is configured to use this Docker container and will automatically retrieve it and run all of the pipeline steps within the Docker container.

1.  Install Python3 and Snakemake. Follow the instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).
2.  Install R packages.
3.  Clone this repository.

# Overview of analysis pipeline

In our manuscript, we generated data using four different genomic assays: ChIP-seq, ATAC-seq, RNA-seq, and CUT&RUN. Analysis of these data can be divided into two stages:

1. Initial processing of data for each assay (read alignment, filtering, peak calling, etc.)
2. Integration of the different genomic assays and downstream analysis. 

# Running the workflow

# Navigating the results
