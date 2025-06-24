# Environment Information for Reproducibility

This document describes the environment setup required to fully reproduce the analysis results of this project.

## Types of Environment Files

### 1. `environment.yml` (Recommended)
**Full reproducibility**: Complete environment specification including build numbers
```bash
conda env create -f environment.yml
```

### 2. `environment-no-builds.yml` (Cross-platform)
**Compatibility-focused**: Improved cross-platform compatibility by omitting build numbers
```bash
conda env create -f environment-no-builds.yml
```

### 3. `conda-packages-list.txt` (Reference)
A detailed list of all installed packages (for reference)

### 4. `R-session-info.txt` (R Environment Details)
R session information and a list of installed package versions

## Environment Reconstruction Steps

### Method 1: Full Reproducibility (Recommended)
```bash
# 1. Create the environment
conda env create -f environment.yml

# 2. Activate the environment
conda activate mieru

# 3. Check functionality
./Fig3-6/scripts/00-setup/01-check-dependencies.sh
```

### Method 2: Cross-Platform Compatibility
```bash
# 1. Create the environment (without build numbers)
conda env create -f environment-no-builds.yml

# 2. Activate the environment
conda activate mieru

# 3. Adjust specific packages as needed
conda install <package-name>=<version>
```

### Method 3: Manual Setup
```bash
# 1. Create the base environment
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mieru -y

# 2. Install required packages
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

# 3. Activate the environment
conda activate mieru
```

## Version Information

This project was developed and tested under the following environment:

- **Created on**: Tue Jun 24 11:50:50 JST 2025
- **Conda version**: conda 25.3.1
- **R version**: 4.3.0
- **Bioconductor version**: 3.17
- **OS**: Linux sycom-2024 6.6.87.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun  5 18:30:46 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux

## Additional Measures for Reproducibility

### Data Integrity
- Verify hash values of data files
- Ensure accurate downloads from the GEO database

### Execution Environment
- Use the same hardware configuration if possible
- Ensure sufficient memory (recommended: 32GB or more)
- Ensure sufficient storage (recommended: at least 100GB free space)

### Reproducibility Check
```bash
# Check dependencies
./Fig3-6/scripts/00-setup/01-check-dependencies.sh

# Run the full pipeline
./Fig3-6/scripts/run-analysis.sh

# Compare results
# - Check consistency of statistics
# - Visually verify figures and tables
```

## Troubleshooting

### Common Issues and Solutions

1. **Package Conflicts**
   ```bash
   conda clean --all
   conda env remove -n mieru
   conda env create -f environment.yml
   ```

2. **R Package Errors**
   ```bash
   conda run -n mieru R -e "update.packages(ask=FALSE)"
   ```

3. **Out of Memory**
   - Increase swap file size
   - Adjust batch sizes

4. **Platform-specific Issues**
   - Use `environment-no-builds.yml`
   - Manually adjust individual packages

## Contact

If you encounter issues during environment reconstruction, please contact us with the following information:
- The environment file you used
- The full error message
- Details of your execution environment (OS, conda version, etc.)

## Change Log

- Tue Jun 24 11:50:50 JST 2025: Initial version created
- Standardized environment file formats
- Added documentation to ensure reproducibility
