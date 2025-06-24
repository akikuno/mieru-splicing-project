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
- **Fig2/**: Initial characterization studies with MIERU cell lines and marker gene analysis
- **Fig3-6/**: Main comparative analysis of RBP knockouts vs controls

### Script Organization

**Fig3-6 scripts are organized by figure and functionality**:
- `015-preprocess/`: Initial data preprocessing pipeline
- `Fig3-event-frequency/`: Alternative splicing event analysis
- `Fig4-compare-to-degs/`: Integration with differential expression
- `Fig5-enrichr/`: Functional enrichment analysis  
- `Fig6-010-complex-human-mouse/`: Cross-species protein complex analysis
- `Fig6-020-complex-heatmap/`: Complex visualization
- `Fig6-030-complex-RBP/`: RBP-specific complex analysis
- `SFig-010-exon-characteristics/`: Supplementary exon feature analysis

**Numbering within each directory**:
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

**Phase 2: Analysis** (figure-specific subdirectories)
- **Fig3**: Event frequency and magnitude analysis (`Fig3-event-frequency/`)
- **Fig4**: Differential expression integration with DESeq2 (`Fig4-compare-to-degs/`)
- **Fig5**: Functional enrichment analysis (`Fig5-enrichr/`)
- **Fig6**: Protein complex analysis using CORUM/ComplexTab databases (`Fig6-*/`)
- **SFig**: Exon characteristics (length, GC content, conservation) (`SFig-010-exon-characteristics/`)

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

**Fig2 Analysis**:
- `Fig2/data/genome/mm39_gene_name.gtf` - GTF annotation for marker genes
- `Fig2/data/counts/` - Read counts and fluorescent protein analysis
- `Fig2/data/markers.tsv` - Marker gene definitions
- `Fig2/reports/` - Visualization outputs (fluorescent analysis, marker gene violin plots)

**Fig3-6 Analysis**:
- `Fig3-6/data/genome/` - Complete genome resources (mm39.fa, GTF, indices)
- `Fig3-6/data/rmats/all_events_ko_target_fdr_dpsi.csv` - Master splicing events table
- `Fig3-6/data/degs/` - Differential expression results
- `Fig3-6/data/counts/featurecounts_gene_name.tsv.gz` - Gene expression counts
- `Fig3-6/data/exon_characteristics/` - Exon features (length, GC, conservation)
- `Fig3-6/data/Fig5/` - Protein complex databases (CORUM, ComplexTab)
- `Fig3-6/data/Fig7/` - Spliceosome analysis data

**Public Data**: Available at GEO accessions GSE291522 (Fig2) and GSE291672 (Fig3-6)

## Analysis Commands

**Fig2 Analysis Commands**:
```bash
# Marker gene analysis
cd Fig2/scripts/035-marker_genes/
./035-mapping.sh
./037-featurecounts.sh
Rscript 045-violinplot_marker_genes.R

# Fluorescent protein analysis
cd Fig2/scripts/025-fluorescents/
./033-mapping_to_fluorescents.sh
./037-count_reads_to_fluorescents.sh
Rscript 045-barplot.R
```

**Fig3-6 Analysis Commands**:
```bash
# Preprocessing pipeline
cd Fig3-6/scripts/015-preprocess/
./015-download_genomes.sh
./027-fastq_trimming.sh
./035-fastq_mapping.sh
./037-featurecounts.sh
./045-rmats.sh

# Fig3: Event frequency analysis
cd Fig3-6/scripts/Fig3-event-frequency/
./015-preprocess.sh
Rscript 025-barplot_frequency_per_event.R
Rscript 035-violinplot_dpsi_per_event.R

# Fig4: Differential expression analysis
cd Fig3-6/scripts/Fig4-compare-to-degs/
Rscript 015-deseq2.R
Rscript 025-venn_vs_expression.R

# Fig5: Enrichment analysis
cd Fig3-6/scripts/Fig5-enrichr/
Rscript 013-enrichr.R
Rscript 015-se_enrichr.R

# Fig6: Complex analysis
cd Fig3-6/scripts/Fig6-010-complex-human-mouse/
Rscript 025-fisher_complexity_complextab_human_mouse.R

# SFig: Exon characteristics
cd Fig3-6/scripts/SFig-010-exon-characteristics/
./015-preprocess.sh
Rscript 045-plot_exon_length_gc_conservation.R
```

## Important Architecture Notes

**Splicing Event Types**: The pipeline detects 5 types of alternative splicing:
- SE (Skipped Exon)
- A3SS (Alternative 3' Splice Site)
- A5SS (Alternative 5' Splice Site) 
- MXE (Mutually Exclusive Exons)
- RI (Retained Intron)

**rMATS Integration**: The central analysis uses rMATS output, which compares each KO condition against MIERU_CH controls. The `015-preprocess.sh` script in Fig3-6 consolidates all rMATS results into a unified CSV format.

**Significance Thresholds**: All splicing analyses use consistent filtering criteria:
- FDR (False Discovery Rate) < 0.05
- |ΔPSI| (absolute Delta Percent Spliced In) > 0.1 (10% change)
- Applied as: `filter(fdr < 0.05, abs(dpsi) > 0.1)`

**Cross-species Analysis**: Several analysis modules incorporate human-mouse orthology mapping for comparative genomics, using databases like MGI homology and UniProt mapping.

**Statistical Framework**: Uses DESeq2 for differential expression, Fisher's exact tests for enrichment analysis, and FDR correction for multiple testing throughout the pipeline.

**Visualization**: Generates publication-ready figures using ggplot2, circlize, and other R visualization packages, with outputs in both PDF and SVG formats.

## Development Notes

When working with this codebase:
- Always activate the `mieru` conda environment before running scripts
- Scripts must be run from their respective directories due to relative path dependencies
- Fig2 focuses on control characterization and marker genes; Fig3-6 contains the main RBP knockout analysis
- Fig3-6 scripts are organized by figure number corresponding to publication figures
- Output goes to `reports/Fig{N}/` directories matching the script directory names
- Large data files (BAM, FASTQ) are stored locally but not in git - use GEO accessions for data access
- R scripts expect specific data file structures - check data preprocessing steps if encountering missing file errors
- All splicing significance testing uses the same thresholds across all analyses for consistency