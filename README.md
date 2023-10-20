# Introduction

This repository contains all analysis code from the following paper:

Protein-intrinsic properties and context-dependent effects regulate pioneer-factor binding and function Tyler J Gibson, Melissa Harrison bioRxiv 2023.03.18.533281; doi: <https://doi.org/10.1101/2023.03.18.533281>

If you are interested in exploring the processed data from this paper, bigWig files and lists of peaks have been deposited in the Gene Expression Omnibus under accession number GSE227884.
Lists of class I, II, and III peaks for Zelda, Grainy head and Twist can be found as supplemental tables along with the published paper.

For anyone interested in reproducing the full analysis pipeline, follow the steps below.

# Setup

We used the workflow management system Snakemake to run all of the analysis for this manuscript.
To get started, follow the steps below to install the necessary software.

1.  Make sure you have installed [R](https://www.r-project.org/) (version 4 or above) and a Conda-based Python3 distribution (I recommend [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)).
2.  Install Snakemake. Follow the instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). I recommend a full Snakemake installation, rather than the minimal one.
3.  Clone this repository:

```         
git clone https://github.com/tjgibson/S2_pioneers_manuscript.git
```

4.  Install R packages Most of the Snakemake pipeline uses Conda to automatically install the necessary software for each step. Because I encountered issues using Conda to install some R packages, I use the R package management system [renv](https://rstudio.github.io/renv/articles/renv.html) to handle installation of additional required R packages.

First, open a new R session and install the renv R package:

```         
install.packages("renv")
```

This repository contains a file named [renv.lock](renv.lock).
This file contains a list of all the R packages used by this workflow.
To automatically install all these packages, use the following code in your R session:

```         
renv::restore()
```

# Overview of analysis pipeline

In our manuscript, we generated data using four different genomic assays: ChIP-seq, ATAC-seq, RNA-seq, and CUT&RUN.
Analysis of these data can be divided into two stages:

1.  Initial processing of data for each assay (read alignment, filtering, peak calling, etc.)
2.  Integration of the different genomic assays and downstream analysis.

To get an overall sense of how the analysis pipeline is structured, see the [Snakefile](Snakefile).

# Running the workflow

To run the entire workflow, from raw read alignment to the generation of figures for the manscript, use the following code to invoke snakemake:

```         
snakemake --cores 1 --use-conda
```

The `--cores` parameter tells snakemake how many cores to use.
Adjust this based on your available computing resources.

The `--use-conda` parameter ensures that Snakemake will use Conda to automatically install the necessary software for each step of the workflow.
Most steps of the workflow will be run inside of a unique Conda environment, providing isolation and reproducibility.
