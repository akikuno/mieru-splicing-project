#!/bin/bash
# Main analysis pipeline for Fig3-6 
# Orchestrates the complete analysis workflow

set -e  # Exit on any error

echo "========================================"
echo "Starting Fig3-6 Analysis Pipeline"
echo "========================================"

# Change to script directory
cd "$(dirname "$0")"

# Check if conda environment is activated
if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "Warning: No conda environment detected. Please activate 'mieru' environment:"
    echo "conda activate mieru"
    exit 1
fi

if [[ "${CONDA_DEFAULT_ENV}" != "mieru" ]]; then
    echo "Warning: Current environment is '${CONDA_DEFAULT_ENV}'. Please activate 'mieru' environment:"
    echo "conda activate mieru"
    exit 1
fi

echo "Using conda environment: ${CONDA_DEFAULT_ENV}"

# Function to run step with logging
run_step() {
    local step_name="$1"
    local script_path="$2"
    local script_type="$3"  # "bash" or "R"
    
    echo ""
    echo "----------------------------------------"
    echo "Running: $step_name"
    echo "Script: $script_path"
    echo "----------------------------------------"
    
    if [[ "$script_type" == "bash" ]]; then
        bash "$script_path"
    elif [[ "$script_type" == "R" ]]; then
        Rscript "$script_path"
    else
        echo "Error: Unknown script type: $script_type"
        exit 1
    fi
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Completed: $step_name"
    else
        echo "✗ Failed: $step_name"
        exit 1
    fi
}

# Step 0: Setup
echo "Phase 0: Setup"
run_step "Check dependencies" "00-setup/01-check-dependencies.sh" "bash"
run_step "Create directories" "00-setup/02-create-directories.sh" "bash"
run_step "Download genome data" "00-setup/03-download-genomes.sh" "bash"

# Step 1: Preprocessing
echo ""
echo "Phase 1: Preprocessing"
run_step "FASTQ trimming" "01-preprocessing/01-fastq-trimming.sh" "bash"
run_step "Read mapping" "01-preprocessing/02-read-mapping.sh" "bash"
run_step "Feature counting" "01-preprocessing/03-feature-counts.sh" "bash"
run_step "Differential splicing analysis" "01-preprocessing/04-differential-splicing.sh" "bash"

# Step 2: Quality control
echo ""
echo "Phase 2: Quality Control"
run_step "Generate FastQC reports" "02-quality-control/01-fastqc-reports.sh" "bash"
run_step "Calculate mapping statistics" "02-quality-control/02-mapping-stats.sh" "bash"
run_step "Plot QC metrics" "02-quality-control/03-plot-qc-metrics.R" "R"

# Step 3: Event analysis (Fig 3)
echo ""
echo "Phase 3: Event Analysis (Figure 3)"
run_step "Preprocess event data" "04-event-analysis/015-preprocess.sh" "bash"
run_step "Event frequency analysis" "04-event-analysis/01-event-frequency.R" "R"
run_step "ΔPSI distribution analysis" "04-event-analysis/02-dpsi-distribution.R" "R"

# Step 4: DEG comparison (Fig 4)
echo ""
echo "Phase 4: DEG Comparison (Figure 4)"
run_step "Differential expression analysis" "05-deg-comparison/015-deseq2.R" "R"
run_step "Overlap analysis" "05-deg-comparison/025-venn_vs_expression.R" "R"
run_step "Expression visualization" "05-deg-comparison/035-expression.R" "R"

# Step 5: Complex analysis (Fig 5)
echo ""
echo "Phase 5: Complex Analysis (Figure 5)"
run_step "Download complex databases" "06-complex-analysis/013-download_uniprot.sh" "bash"
run_step "Complex enrichment analysis" "06-complex-analysis/015-download_complextab.R" "R"
run_step "Fisher's exact tests" "06-complex-analysis/025-fisher_complexity_complextab_human_mouse.R" "R"

# Step 6: Heatmap analysis (Fig 6)
echo ""
echo "Phase 6: Heatmap Analysis (Figure 6)"
run_step "Download GO annotations" "07-heatmap-analysis/013-download-go-human-mouse.sh" "bash"
run_step "Generate heatmaps" "07-heatmap-analysis/025-heatmap-go-human-mouse.R" "R"
run_step "RBP complex analysis" "07-heatmap-analysis/045-heatmap.R" "R"

# Supplementary: Exon characteristics
echo ""
echo "Supplementary: Exon Characteristics"
run_step "Extract exon data" "SFig-exon-characteristics/01-extract-exon-data.sh" "bash"
run_step "Plot exon characteristics" "SFig-exon-characteristics/04-plot-exon-characteristics.R" "R"

echo ""
echo "========================================"
echo "Analysis Pipeline Completed Successfully!"
echo "========================================"
echo ""
echo "Output files have been generated in:"
echo "  - reports/Fig3/"
echo "  - reports/Fig4/"
echo "  - reports/Fig5/"
echo "  - reports/Fig6/"
echo "  - reports/SFig/"
echo ""
echo "Next steps:"
echo "  1. Review generated plots and statistics"
echo "  2. Check log files for any warnings"
echo "  3. Run quality control checks on results"