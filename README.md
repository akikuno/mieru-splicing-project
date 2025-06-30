[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/862646954.svg)](https://doi.org/10.5281/zenodo.13922578)


# Mieru Splicing Project

This repository contains comprehensive RNA splicing analysis scripts and data for studying the effects of RNA-binding protein (RBP) knockouts on alternative splicing patterns in mouse embryonic stem cells. The project analyzes differential splicing events across 11 different RBP knockout lines compared to MIERU control cells.

## Project Overview

The analysis investigates how specific RBP knockouts affect:
- Alternative splicing patterns (5 event types: SE, A3SS, A5SS, MXE, RI)
- Gene expression changes
- Protein complex composition
- Functional pathway enrichment

**Key RBPs analyzed**: Cd2bp2, Qk, Rbm24, Rpl22l1, Spen, Strap, Tra2b, Trim71, Ubr5, Wt1, Ybx1

## Requirements

- Unix environment (WSL2 with Ubuntu or macOS recommended)
- conda package manager via [miniforge](https://github.com/conda-forge/miniforge)
- ~50GB storage space for genome indices and intermediate files

## Installation

Create and activate the conda environment with all required bioinformatics tools:

```bash
conda env create -f environment.yml
conda activate mieru
```


> [!IMPORTANT]
> To ensure full reproducibility of this project's analysis results, the following files are provided:
> - **`environment.yml`**: Complete conda environment (including build numbers)
> - **`environment-no-builds.yml`**: Cross-platform compatible version
> - **`conda-packages-list.txt`**: Detailed list of installed packages
> - **`R-session-info.txt`**: R session information
> - **`REPRODUCIBILITY.md`**: Step-by-step instructions for reproducibility  
> For details, see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).


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

The analysis pipeline follows a modular, standardized organization:

### Fig2 Analysis (Control Characterization)
- **Fluorescent Analysis** (`025-fluorescents/`): Fluorescent protein expression quantification
- **Marker Gene Analysis** (`035-marker_genes/`): Validation of cell line markers

### Fig3-6 Analysis (Main Pipeline)
1. **Setup** (`00-setup/`): System dependencies, directory creation, and genome data download
2. **Preprocessing** (`01-preprocessing/`): Quality control, trimming, genome alignment, read counting, and differential splicing with rMATS
3. **Quality Control** (`02-quality-control/`): Exon characteristics analysis and quality metrics
4. **Event Analysis** (`04-event-analysis/`): Alternative splicing event frequency and ΔPSI distribution analysis (Figure 3)
5. **DEG Comparison** (`05-deg-comparison/`): Integration with differential gene expression (Figure 4)
6. **Complex Analysis** (`06-complex-analysis/`): Protein complex enrichment using CORUM/ComplexTab databases (Figure 5)
7. **Heatmap Analysis** (`07-heatmap-analysis/`): GO term and pathway visualization (Figure 6)

### Legacy Structure
- Original organization preserved in `_past/` directories for backward compatibility
- Legacy scripts (e.g., `015-preprocess/`, `Fig3-event-frequency/`) remain functional

## Usage

1. Download raw data from GEO repositories
2. Activate conda environment: `conda activate mieru`
3. Run preprocessing scripts in numerical order within each analysis directory
4. Scripts must be executed from their respective directories due to relative paths

