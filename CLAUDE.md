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
    r-enrichr r-ggVennDiagram r-circlize r-pheatmap r-readxl \
    r-yaml r-rcolorbrewer \
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

### Script Organization (Refactored Structure)

**Fig3-6 scripts follow a standardized, modular organization**:
- `00-setup/`: System setup, dependency checks, and genome data download
- `01-preprocessing/`: Core data preprocessing pipeline (trimming, mapping, counting, splicing)
- `02-quality-control/`: Quality metrics and exon characteristics analysis
- `03-event-analysis/`: Alternative splicing event analysis (Figure 3)
- `04-deg-comparison/`: Integration with differential expression (Figure 4)
- `05-complex-analysis/`: Protein complex enrichment analysis (Figure 5)
- `06-heatmap-analysis/`: GO term and pathway visualization (Figure 6)
- `utils/`: Shared utility functions (data loading, plotting, statistics)
- `config/`: Configuration files and parameters

**Utility Modules**:
- `utils/data-loaders.R`: Standardized data loading with consistent filtering
- `utils/plot-themes.R`: Reusable ggplot2 themes and color schemes
- `utils/statistical-tests.R`: Common statistical analysis functions
- `utils/config-loader.R`: Configuration management

**Legacy Structure** (preserved for reference):
- `015-preprocess/`, `Fig3-event-frequency/`, etc.: Original organization
- `_past/`: Archived experimental approaches

### Standard Bioinformatics Pipeline

**Refactored Pipeline Structure**:

**Phase 0: Setup** (`00-setup/`)
- `01-check-dependencies.sh` - Verify system dependencies and R packages
- `02-create-directories.sh` - Create required directory structure
- `03-download-genomes.sh` - Download mm39 genome and GTF from Ensembl release-111

**Phase 1: Preprocessing** (`01-preprocessing/`)
- `01-fastq-trimming.sh` - Quality control with fastp
- `02-read-mapping.sh` - STAR alignment to genome
- `03-feature-counts.sh` - Gene-level read counting
- `04-differential-splicing.sh` - Alternative splicing detection with rMATS

**Phase 2: Quality Control** (`02-quality-control/`)
- Exon characteristics analysis (length, GC content, conservation)

**Phase 3: Event Analysis** (`03-event-analysis/`)
- **Figure 3**: Event frequency and ΔPSI distribution analysis
- `01-event-frequency.R` - Barplot of event type percentages
- `02-dpsi-distribution.R` - Violin plots of ΔPSI distributions

**Phase 4: DEG Comparison** (`04-deg-comparison/`)
- **Figure 4**: Integration with differential gene expression
- Overlap analysis between splicing and expression changes

**Phase 5: Complex Analysis** (`05-complex-analysis/`)
- **Figure 5**: Protein complex enrichment using CORUM/ComplexTab databases

**Phase 6: Heatmap Analysis** (`06-heatmap-analysis/`)
- **Figure 6**: GO term and pathway heatmaps

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

**Refactored Analysis Commands**:

**Complete Pipeline** (Recommended):
```bash
cd Fig3-6/scripts/
conda activate mieru
./run-analysis.sh
```

**Individual Phases**:
```bash
cd Fig3-6/scripts/

# Phase 0: Setup
./00-setup/01-check-dependencies.sh
./00-setup/02-create-directories.sh
./00-setup/03-download-genomes.sh

# Phase 1: Preprocessing
./01-preprocessing/01-fastq-trimming.sh
./01-preprocessing/02-read-mapping.sh
./01-preprocessing/03-feature-counts.sh
./01-preprocessing/04-differential-splicing.sh

# Phase 2: Quality Control
./02-quality-control/015-preprocess.sh
Rscript 02-quality-control/045-plot_exon_length_gc_conservation.R

# Phase 3: Event Analysis (Figure 3)
Rscript 03-event-analysis/01-event-frequency.R
Rscript 03-event-analysis/02-dpsi-distribution.R

# Phase 4: DEG Comparison (Figure 4)
Rscript 04-deg-comparison/015-deseq2.R
Rscript 04-deg-comparison/025-venn_vs_expression.R

# Phase 5: Complex Analysis (Figure 5)
Rscript 05-complex-analysis/025-fisher_complexity_complextab_human_mouse.R

# Phase 6: Heatmap Analysis (Figure 6)
Rscript 06-heatmap-analysis/025-heatmap-go-human-mouse.R
```

**Legacy Commands** (still functional):
```bash
# Original structure (preserved for compatibility)
cd Fig3-6/scripts/015-preprocess/
./015-download_genomes.sh
# ... etc
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
- **Use refactored structure**: New modular organization in `00-setup/`, `01-preprocessing/`, etc.
- **Run complete pipeline**: Use `./run-analysis.sh` for full analysis workflow
- **Configuration-driven**: Modify parameters in `config/parameters.yaml` instead of hard-coding
- **Shared utilities**: Use functions from `utils/` for consistent data loading and plotting
- Fig2 focuses on control characterization; Fig3-6 contains main RBP knockout analysis
- Legacy structure preserved for backward compatibility in original directory names
- Output goes to `reports/Fig{N}/` directories matching publication figures
- Large data files (BAM, FASTQ) are stored locally but not in git - use GEO accessions for data access
- All splicing significance testing uses consistent thresholds: FDR < 0.05, |ΔPSI| > 0.1