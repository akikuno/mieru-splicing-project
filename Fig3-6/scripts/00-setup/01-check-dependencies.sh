#!/bin/bash
# Check system dependencies for Fig3-6 analysis pipeline

set -e

echo "========================================"
echo "Checking System Dependencies"
echo "========================================"

# Check if conda environment is activated
if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "❌ No conda environment detected"
    echo "Please activate the mieru environment: conda activate mieru"
    exit 1
fi

if [[ "${CONDA_DEFAULT_ENV}" != "mieru" ]]; then
    echo "❌ Wrong conda environment: ${CONDA_DEFAULT_ENV}"
    echo "Please activate the mieru environment: conda activate mieru"
    exit 1
fi

echo "✓ Conda environment: ${CONDA_DEFAULT_ENV}"

# Check required command-line tools
required_tools=("fastp" "STAR" "samtools" "bedtools" "featureCounts" "rmats.py" "rsem-calculate-expression")

for tool in "${required_tools[@]}"; do
    if command -v "$tool" &> /dev/null; then
        echo "✓ $tool found"
    else
        echo "❌ $tool not found"
        missing_tools=1
    fi
done

# Check R packages
echo ""
echo "Checking R packages..."
Rscript -e "
required_packages <- c('tidyverse', 'DESeq2', 'ggplot2', 'extrafont', 'patchwork', 'circlize', 'ggVennDiagram', 'pheatmap', 'ggsignif', 'RColorBrewer', 'enrichR', 'org.Mm.eg.db', 'clusterProfiler', 'readxl', 'janitor', 'yaml')
missing_packages <- c()

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        missing_packages <- c(missing_packages, pkg)
        cat('❌', pkg, 'not found\n')
    } else {
        cat('✓', pkg, 'found\n')
    }
}

if (length(missing_packages) > 0) {
    cat('\nMissing R packages:', paste(missing_packages, collapse = ', '), '\n')
    cat('Install with: conda install -n mieru', paste(missing_packages, collapse = ' '), '\n')
    quit(status = 1)
}
"

if [[ $? -ne 0 ]]; then
    echo "❌ Some R packages are missing"
    exit 1
fi

# Check directory structure
required_dirs=("data" "reports" "scripts/utils" "scripts/config")

echo ""
echo "Checking directory structure..."
for dir in "${required_dirs[@]}"; do
    if [[ -d "$dir" ]]; then
        echo "✓ $dir exists"
    else
        echo "⚠️ $dir missing - will be created"
        mkdir -p "$dir"
    fi
done

echo ""
echo "========================================"
echo "✓ All dependencies checked successfully!"
echo "========================================"