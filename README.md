# Mieru Splicing Project

This repository contains comprehensive RNA splicing analysis scripts and data for studying the effects of RNA-binding protein (RBP) knockouts on alternative splicing patterns in mouse embryonic stem cells. The project analyzes differential splicing events across 11 different RBP knockout lines compared to MIERU control cells.

## Project Overview

The analysis investigates how specific RBP knockouts affect:
- Alternative splicing patterns (5 event types: SE, A3SS, A5SS, MXE, RI)
- Gene expression changes
- Protein complex composition
- Functional pathway enrichment
- Splicing factor networks

**Key RBPs analyzed**: Cd2bp2, Qk, Rbm24, Rpl22l1, Spen, Strap, Tra2b, Trim71, Ubr5, Wt1, Ybx1

## Requirements

- Unix environment (WSL2 with Ubuntu or macOS recommended)
- conda package manager via [miniforge](https://github.com/conda-forge/miniforge)
- ~50GB storage space for genome indices and intermediate files

## Installation

Create and activate the conda environment with all required bioinformatics tools:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mieru -y
conda install -n mieru -y \
    fastp star samtools bedtools subread rmats rsem \
    numpy pandas matplotlib seaborn plotnine \
    r-base r-essentials r-extrafont r-janitor \
    r-ggfortify r-ggrepel r-patchwork r-ggsignif r-svglite \
    r-enrichr r-ggVennDiagram r-circlize \
    bioconductor-deseq2 \
    bioconductor-genomeinfodbdata \
    bioconductor-org.hs.eg.db \
    bioconductor-org.Mm.eg.db \
    bioconductor-clusterprofiler

conda activate mieru
```

## Dataset Access

Raw sequencing data and processed results are publicly available:

- **Fig2 data**: Control samples and marker gene analysis  
  - FASTQ files: [GSE291522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE291522)
  - Location: `Fig2/data/fastq/`

- **Fig3-6 data**: RBP knockout comparative analysis  
  - FASTQ files: [GSE291672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE291672)
  - rMATS splicing results: `GSE291672_rmats.zip`
  - Location: `Fig3-6/data/fastq/` and `Fig3-6/data/rmats/original_output/`

## Pipeline Overview

The analysis pipeline consists of:

1. **Preprocessing** (`015-preprocess/`): Quality control, genome alignment, read counting
2. **Splicing Analysis** (`025-event-frequency/`): rMATS-based alternative splicing detection
3. **Functional Analysis** (`035-compare-to-degs/`): Integration with gene expression changes
4. **Enrichment Analysis** (`045-enrichr/`): Pathway and functional annotation
5. **Complex Analysis** (`053-complex-*/`): Protein complex and cross-species analysis

## Usage

1. Download raw data from GEO repositories
2. Activate conda environment: `conda activate mieru`
3. Run preprocessing scripts in numerical order within each analysis directory
4. Scripts must be executed from their respective directories due to relative paths

## Key Tools Used

- **Alignment**: STAR, RSEM
- **Splicing**: rMATS  
- **Statistics**: DESeq2, R packages
- **Visualization**: ggplot2, circlize
