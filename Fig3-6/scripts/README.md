# Fig3-6 Analysis Scripts - Refactored Structure

This directory contains the refactored and standardized analysis pipeline for Figures 3-6 of the mieru-splicing-project.

## Overview

The analysis has been reorganized into a modular, maintainable structure with:
- Standardized naming conventions
- Shared utility functions
- Configuration-driven parameters
- Clear workflow orchestration

## Directory Structure

```
scripts/
├── 00-setup/                 # System setup, dependencies, and genome download
├── 01-preprocessing/          # Core data preprocessing pipeline
├── 02-quality-control/       # Quality control and validation
├── 04-event-analysis/        # Alternative splicing event analysis (Fig 3)
├── 05-deg-comparison/         # Differential expression comparison (Fig 4)
├── 06-complex-analysis/       # Protein complex analysis (Fig 5)
├── 07-heatmap-analysis/       # Heatmap and GO analysis (Fig 6)
├── SFig-exon-characteristics/ # Supplementary exon characteristics analysis
├── utils/                     # Shared utility functions
├── config/                    # Configuration files
├── run-analysis.sh           # Main pipeline orchestration script
└── README.md                 # This file
```

### Legacy Directories (Preserved)
- `015-preprocess/`, `Fig3-event-frequency/`, etc. - Original structure (preserved for reference)
- `_past/` - Archived experimental approaches

## Key Improvements

### 1. Modular Utilities (`utils/`)
- **data-loaders.R**: Standardized data loading with consistent filtering
- **plot-themes.R**: Reusable plotting themes and color schemes
- **statistical-tests.R**: Common statistical analysis functions
- **config-loader.R**: Configuration management utilities

### 2. Configuration Management (`config/`)
- **parameters.yaml**: Centralized parameter configuration
  - Statistical thresholds (FDR, ΔPSI, etc.)
  - Plotting settings (colors, dimensions, fonts)
  - File paths and directories
  - Analysis parameters (KO genes, event types)

### 3. Standardized Workflows
Each analysis phase follows consistent patterns:
- Configuration-driven parameters
- Shared utility functions
- Standardized output formats
- Consistent error handling

## Usage

### Quick Start
Run the complete analysis pipeline:
```bash
conda activate mieru
cd Fig3-6/scripts/
./run-analysis.sh
```

### Individual Steps
Run specific analysis phases:
```bash
# Event frequency analysis (Fig 3)
Rscript 03-event-analysis/01-event-frequency.R
Rscript 03-event-analysis/02-dpsi-distribution.R

# DEG comparison (Fig 4)
Rscript 04-deg-comparison/015-deseq2.R
Rscript 04-deg-comparison/025-venn_vs_expression.R
```

### Customization
Modify analysis parameters in `config/parameters.yaml`:
```yaml
thresholds:
  fdr: 0.05      # Adjust significance threshold
  dpsi: 0.1      # Adjust effect size threshold

plotting:
  formats: ["pdf", "svg"]  # Change output formats
```

## Analysis Phases

### Phase 0: Setup (`00-setup/`)
- System dependency checks
- Directory structure creation  
- Download reference genome and annotations

### Phase 1: Preprocessing (`01-preprocessing/`)
- Quality control and read trimming with fastp
- Genome alignment with STAR
- Read counting with featureCounts
- Alternative splicing analysis with rMATS

### Phase 2: Quality Control (`02-quality-control/`)
- FastQC reports for read quality assessment
- Mapping statistics and alignment metrics
- Library size and count distribution analysis

### Phase 3: Event Analysis (`04-event-analysis/`)
- **Figure 3**: Alternative splicing event frequency and distribution
- Event type classification and quantification
- ΔPSI distribution analysis

### Phase 4: DEG Comparison (`05-deg-comparison/`)
- **Figure 4**: Integration with differential gene expression
- Overlap analysis between splicing and expression changes
- Functional enrichment of overlapping genes

### Phase 5: Complex Analysis (`06-complex-analysis/`)
- **Figure 5**: Protein complex enrichment analysis
- Fisher's exact tests for complex membership
- Cross-species (human-mouse) complex analysis

### Phase 6: Heatmap Analysis (`07-heatmap-analysis/`)
- **Figure 6**: GO term and pathway heatmaps
- RBP-specific complex analysis
- Functional annotation visualization

### Supplementary Analysis (`SFig-exon-characteristics/`)
- **Supplementary Figure**: Exon characteristics analysis
- Exon length, GC content, and conservation analysis
- Quality metrics for alternative splicing events

## Output Structure

Results are organized by figure:
```
reports/
├── Fig3/          # Event frequency and distribution plots
├── Fig4/          # DEG comparison and overlap analysis
├── Fig5/          # Complex analysis results
├── Fig6/          # Heatmaps and functional analysis
└── SFig/          # Supplementary quality control plots
```

## Configuration

### Key Parameters
- **FDR threshold**: 0.05 (adjustable in config)
- **ΔPSI threshold**: 0.1 (10% change, adjustable)
- **Event types**: A3SS, A5SS, MXE, RI, SE
- **KO genes**: 11 RBP knockouts (Cd2bp2, Qk, Rbm24, etc.)

### File Paths
All paths are configurable in `config/parameters.yaml`:
- Data directories (`data/rmats/`, `data/degs/`, etc.)
- Output directories (`reports/Fig3/`, etc.)
- External database URLs

## Dependencies

### R Packages
```r
# Core analysis
tidyverse, DESeq2, yaml

# Plotting
ggplot2, extrafont, patchwork, circlize, ggVennDiagram, pheatmap, ggsignif, RColorBrewer

# Functional analysis
enrichR, org.Mm.eg.db, clusterProfiler

# Data manipulation
readxl, janitor
```

### External Tools
- STAR (alignment)
- featureCounts (quantification) 
- rMATS (splicing analysis)
- fastp (quality control)

## Maintenance Notes

### Adding New Analyses
1. Create new script in appropriate phase directory
2. Follow naming convention: `##-descriptive-name.R`
3. Use shared utilities from `utils/`
4. Add to `run-analysis.sh` if part of main pipeline

### Modifying Parameters
1. Update `config/parameters.yaml`
2. Test with subset of data
3. Update documentation if behavior changes

## Troubleshooting

### Common Issues
1. **Missing conda environment**: Activate with `conda activate mieru`
2. **File not found errors**: Check that preprocessing steps completed successfully
3. **Memory issues**: Reduce parallelization in compute-intensive steps
4. **Plot rendering issues**: Ensure X11 forwarding or run in headless mode

### Getting Help
1. Check individual script logs for specific errors
2. Verify all dependencies are installed
3. Ensure input data files exist and are readable
4. Consult original scripts in legacy directories for comparison
