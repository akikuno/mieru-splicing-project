#!/bin/bash
# Generate FastQC reports for quality assessment

set -e

echo "========================================"
echo "Generating FastQC Quality Reports"
echo "========================================"

# Create output directory
mkdir -p "reports/QC/fastqc"

# Input and output directories
FASTQ_DIR="data/fastq_trimmed"
OUTPUT_DIR="reports/QC/fastqc"

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: FASTQ directory not found: $FASTQ_DIR"
    echo "Please run preprocessing steps first"
    exit 1
fi

# Check if FastQC is available
if ! command -v fastqc &> /dev/null; then
    echo "Error: FastQC not found. Please install FastQC"
    exit 1
fi

echo "Input directory: $FASTQ_DIR"
echo "Output directory: $OUTPUT_DIR"

# Count FASTQ files
FASTQ_COUNT=$(find "$FASTQ_DIR" -name "*.fq.gz" | wc -l)
echo "Found $FASTQ_COUNT FASTQ files to process"

if [[ $FASTQ_COUNT -eq 0 ]]; then
    echo "Warning: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

# Run FastQC on all trimmed FASTQ files
echo ""
echo "Running FastQC analysis..."

fastqc \
    --outdir "$OUTPUT_DIR" \
    --threads 4 \
    --format fastq \
    "$FASTQ_DIR"/*.fq.gz

echo ""
echo "✓ FastQC analysis completed"
echo "Reports saved to: $OUTPUT_DIR"
echo ""
echo "Summary:"
echo "- HTML reports: $OUTPUT_DIR/*.html"
echo "- ZIP archives: $OUTPUT_DIR/*.zip"