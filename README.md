# mieru-splicing-project

This repository contains scripts and figures related to the mieru-splicing-project.

## Environment

- Unix environment such as WSL2 with Ubuntu or macOS is required.
- Install conda via [miniforge](https://github.com/conda-forge/miniforge)

```python
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mieru -y python=3.10
conda install -n mieru -y \
    fastp star samtools bedtools subread rmats \
    r-base r-essentials r-extrafont \
    r-ggfortify r-ggrepel r-patchwork r-enrichr r-ggVennDiagram \
    bioconductor-deseq2

conda activate mieru
```

## Data storage locations

- `Fig1/data/fastq` (GSE XXXXX)
- `Fig2-4/data/fastq` (GSE YYYYY)
- `Fig2-4/data/rmats/original_output` (GSE YYYYY: contents of rMATS.zip)
