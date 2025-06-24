#!/bin/bash
# Record detailed environment information for reproducibility

set -e

echo "========================================"
echo "Recording Environment Information"
echo "========================================"

# Create output directory
mkdir -p "logs/environment"

# Get current timestamp
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
OUTPUT_DIR="logs/environment"

# Record system information
echo "Recording system information..."
{
    echo "Environment Recording Session"
    echo "============================"
    echo "Timestamp: $(date)"
    echo "User: $(whoami)"
    echo "Working Directory: $(pwd)"
    echo ""
    
    echo "System Information:"
    echo "---------------"
    uname -a
    echo ""
    
    echo "CPU Information:"
    echo "-------------"
    lscpu | head -20
    echo ""
    
    echo "Memory Information:"
    echo "----------------"
    free -h
    echo ""
    
    echo "Disk Usage:"
    echo "----------"
    df -h
    echo ""
    
} > "$OUTPUT_DIR/system_info_${TIMESTAMP}.txt"

# Record conda environment
echo "Recording conda environment..."
if command -v conda &> /dev/null; then
    {
        echo "Conda Information:"
        echo "=================="
        conda info
        echo ""
        
        echo "Active Environment:"
        echo "=================="
        conda list
        echo ""
        
    } > "$OUTPUT_DIR/conda_env_${TIMESTAMP}.txt"
    
    # Export current environment
    conda env export > "$OUTPUT_DIR/environment_${TIMESTAMP}.yml"
else
    echo "Warning: conda not found"
fi

# Record R session information
echo "Recording R session information..."
if command -v Rscript &> /dev/null; then
    Rscript -e "
    sink('$OUTPUT_DIR/r_session_${TIMESTAMP}.txt')
    cat('R Session Information\n')
    cat('=====================\n\n')
    sessionInfo()
    cat('\n\nLoaded Packages:\n')
    cat('===============\n\n')
    search()
    cat('\n\nInstalled Packages:\n')
    cat('==================\n\n')
    ip <- installed.packages()
    write.csv(ip[,c('Package', 'Version', 'Built')], 
              file = '$OUTPUT_DIR/r_packages_${TIMESTAMP}.csv', 
              row.names = FALSE)
    sink()
    "
else
    echo "Warning: R/Rscript not found"
fi

# Record Git information
echo "Recording Git information..."
if git rev-parse --git-dir > /dev/null 2>&1; then
    {
        echo "Git Repository Information:"
        echo "=========================="
        echo "Repository: $(git config --get remote.origin.url 2>/dev/null || echo 'Local repository')"
        echo "Branch: $(git branch --show-current 2>/dev/null || echo 'Unknown')"
        echo "Commit: $(git rev-parse HEAD 2>/dev/null || echo 'Unknown')"
        echo "Status:"
        git status --porcelain
        echo ""
        
        echo "Recent Commits:"
        echo "==============" 
        git log --oneline -10
        echo ""
        
    } > "$OUTPUT_DIR/git_info_${TIMESTAMP}.txt"
else
    echo "Warning: Not a git repository"
fi

# Record environment variables
echo "Recording environment variables..."
{
    echo "Environment Variables:"
    echo "====================="
    env | sort
} > "$OUTPUT_DIR/env_vars_${TIMESTAMP}.txt"

# Record file checksums for reproducibility
echo "Recording data file checksums..."
if [ -d "data" ]; then
    find data -type f -name "*.gz" -o -name "*.fa" -o -name "*.gtf" -o -name "*.tsv" -o -name "*.csv" | \
    head -20 | \
    xargs md5sum > "$OUTPUT_DIR/data_checksums_${TIMESTAMP}.md5" 2>/dev/null || \
    echo "No data files found for checksum calculation"
fi

# Create summary file
{
    echo "Reproducibility Information Summary"
    echo "=================================="
    echo "Generated: $(date)"
    echo "Session ID: ${TIMESTAMP}"
    echo ""
    echo "Files created:"
    ls -la "$OUTPUT_DIR"/*_${TIMESTAMP}.*
    echo ""
    echo "For complete reproducibility, save all files in this directory"
    echo "and use the corresponding environment.yml file to recreate the environment."
} > "$OUTPUT_DIR/summary_${TIMESTAMP}.txt"

echo ""
echo "✓ Environment information recorded successfully!"
echo "Output directory: $OUTPUT_DIR"
echo "Session ID: $TIMESTAMP"
echo ""
echo "Files created:"
ls -la "$OUTPUT_DIR"/*_${TIMESTAMP}.*