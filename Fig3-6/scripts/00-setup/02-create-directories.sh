#!/bin/bash
# Create directory structure for Fig3-6 analysis pipeline

set -e

echo "========================================"
echo "Creating Directory Structure"
echo "========================================"

# Change to project root
cd "$(dirname "$0")/../.."

# Create data directories
data_dirs=(
    "data/genome"
    "data/fastq_trimmed"
    "data/bam"
    "data/counts" 
    "data/degs"
    "data/rmats"
    "data/exon_characteristics"
    "data/Fig5"
    "data/Fig6"
    "data/Fig7"
    "data/rsem/bam"
)

echo "Creating data directories..."
for dir in "${data_dirs[@]}"; do
    mkdir -p "$dir"
    echo "✓ Created: $dir"
done

# Create reports directories
reports_dirs=(
    "reports/Fig3"
    "reports/Fig4"
    "reports/Fig5" 
    "reports/Fig6"
    "reports/Fig7"
    "reports/SFig"
)

echo ""
echo "Creating reports directories..."
for dir in "${reports_dirs[@]}"; do
    mkdir -p "$dir"
    echo "✓ Created: $dir"
done

# Create temporary directories
temp_dirs=(
    "tmp"
    "logs"
)

echo ""
echo "Creating temporary directories..."
for dir in "${temp_dirs[@]}"; do
    mkdir -p "$dir"
    echo "✓ Created: $dir"
done

echo ""
echo "========================================"
echo "✓ Directory structure created successfully!"
echo "========================================"