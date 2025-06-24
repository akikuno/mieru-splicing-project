#!/bin/bash
# Generate mapping statistics and quality metrics

set -e

echo "========================================"
echo "Generating Mapping Statistics"
echo "========================================"

# Create output directory
mkdir -p "reports/QC/mapping"

# Input and output directories
BAM_DIR="data/bam"
OUTPUT_DIR="reports/QC/mapping"

if [[ ! -d "$BAM_DIR" ]]; then
    echo "Error: BAM directory not found: $BAM_DIR"
    echo "Please run read mapping first"
    exit 1
fi

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found"
    exit 1
fi

echo "Input directory: $BAM_DIR"
echo "Output directory: $OUTPUT_DIR"

# Initialize summary file
SUMMARY_FILE="$OUTPUT_DIR/mapping_summary.txt"
echo "Sample,Total_Reads,Mapped_Reads,Mapping_Rate,Properly_Paired" > "$SUMMARY_FILE"

echo ""
echo "Processing BAM files..."

# Process each BAM file
for bam_file in "$BAM_DIR"/*.bam; do
    if [[ -f "$bam_file" ]]; then
        sample_name=$(basename "$bam_file" .bam)
        echo "Processing: $sample_name"
        
        # Generate detailed stats
        samtools flagstat "$bam_file" > "$OUTPUT_DIR/${sample_name}_flagstat.txt"
        samtools idxstats "$bam_file" > "$OUTPUT_DIR/${sample_name}_idxstats.txt"
        
        # Extract key metrics
        total_reads=$(samtools view -c "$bam_file")
        mapped_reads=$(samtools view -c -F 4 "$bam_file")
        properly_paired=$(samtools view -c -f 2 "$bam_file")
        
        # Calculate mapping rate
        if [[ $total_reads -gt 0 ]]; then
            mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
        else
            mapping_rate="0"
        fi
        
        # Add to summary
        echo "$sample_name,$total_reads,$mapped_reads,$mapping_rate%,$properly_paired" >> "$SUMMARY_FILE"
    fi
done

echo ""
echo "✓ Mapping statistics completed"
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Files generated:"
echo "- Summary: $SUMMARY_FILE"
echo "- Individual flagstat: $OUTPUT_DIR/*_flagstat.txt"
echo "- Individual idxstats: $OUTPUT_DIR/*_idxstats.txt"