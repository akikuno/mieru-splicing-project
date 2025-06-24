# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a comprehensive RNA splicing analysis project studying the effects of RNA-binding protein (RBP) knockouts on alternative splicing patterns in mouse embryonic stem cells. The project analyzes differential splicing events and gene expression across 11 different RBP knockout lines compared to MIERU_CH controls.

## Environment Setup

**Required Environment**: Unix environment (WSL2/Ubuntu or macOS) with conda

**Environment Installation**:
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

## Pipeline Architecture

The analysis is organized into two main sections:
- **Fig1/**: Initial characterization studies with MIERU cell lines
- **Fig2-4/**: Main comparative analysis of RBP knockouts vs controls

### Script Organization

Scripts are numbered by execution order within each analysis phase:
- `015-`: Initial setup/preprocessing
- `025-`: Primary analysis
- `035-`: Secondary analysis
- `045-`: Visualization/output generation

### Standard Bioinformatics Pipeline

**Phase 1: Preprocessing** (`015-preprocess/`)
1. `015-download_genomes.sh` - Downloads mm39 genome and GTF from Ensembl release-111
2. `027-fastq_trimming.sh` - Quality control with fastp
3. `035-fastq_mapping.sh` - STAR alignment to genome
4. `037-featurecounts.sh` - Gene-level read counting
5. `045-rmats.sh` - Alternative splicing detection with rMATS

**Phase 2: Analysis** (various subdirectories)
- Event frequency and magnitude analysis
- Exon characteristics (length, GC content, conservation)
- Differential expression integration with DESeq2
- Functional enrichment analysis
- Protein complex analysis using CORUM/ComplexTab databases

## Data Organization

**Sample Naming Convention**:
- Controls: `MIERU_CH_{replicate}` (e.g., `MIERU_CH_1`)
- Knockouts: `{Gene}_KO_{replicate}` (e.g., `Cd2bp2_KO_5`)

**Key RBP Genes Analyzed**: Cd2bp2, Qk, Rbm24, Rpl22l1, Spen, Strap, Tra2b, Trim71, Ubr5, Wt1, Ybx1

**Data Flow**:
- Input: FASTQ files (paired-end RNA-seq)
- Intermediate: BAM files, count matrices, rMATS output
- Output: Statistical results (CSV), publication figures (PDF/SVG)

## Key File Locations

**Genome Data**: `data/genome/`
- `mm39.fa` - Mouse genome (GRCm39)
- `mm39.gtf` - Ensembl annotation
- `star_index/` - STAR genome index
- `rsem_index/` - RSEM transcriptome index

**Results Data**: `data/`
- `rmats/all_events_ko_target_fdr_dpsi.csv` - Master splicing events table
- `degs/` - Differential expression results
- `counts/featurecounts_gene_name.tsv.gz` - Gene expression counts

**Public Data**: Available at GEO accessions GSE291522 (Fig1) and GSE291672 (Fig2-4)

## Analysis Commands

**Run preprocessing pipeline**:
```bash
cd Fig2-4/scripts/015-preprocess/
./015-download_genomes.sh
./027-fastq_trimming.sh
./035-fastq_mapping.sh
./037-featurecounts.sh
./045-rmats.sh
```

**Generate event frequency analysis**:
```bash
cd Fig2-4/scripts/025-event-frequency/
./015-preprocess.sh
Rscript 025-barplot_frequency_per_event.R
```

**Run differential expression analysis**:
```bash
cd Fig2-4/scripts/035-compare-to-degs/
Rscript 015-deseq2.R
```

## Important Architecture Notes

**Splicing Event Types**: The pipeline detects 5 types of alternative splicing:
- SE (Skipped Exon)
- A3SS (Alternative 3' Splice Site)
- A5SS (Alternative 5' Splice Site) 
- MXE (Mutually Exclusive Exons)
- RI (Retained Intron)

**rMATS Integration**: The central analysis uses rMATS output, which compares each KO condition against MIERU_CH controls. The `015-preprocess.sh` script consolidates all rMATS results into a unified CSV format.

**Cross-species Analysis**: Several analysis modules incorporate human-mouse orthology mapping for comparative genomics, using databases like MGI homology and UniProt mapping.

**Statistical Framework**: Uses DESeq2 for differential expression, Fisher's exact tests for enrichment analysis, and FDR correction for multiple testing throughout the pipeline.

**Visualization**: Generates publication-ready figures using ggplot2, circlize, and other R visualization packages, with outputs in both PDF and SVG formats.

## Development Notes

When working with this codebase:
- Always activate the `mieru` conda environment before running scripts
- Scripts must be run from their respective directories due to relative path dependencies
- Large data files (BAM, FASTQ) are stored locally but not in git - use GEO accessions for data access
- R scripts expect specific data file structures - check data preprocessing steps if encountering missing file errors