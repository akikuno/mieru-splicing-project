#!/bin/bash

########################################
# Make directories
########################################

mkdir -p data/genome

########################################
# Download genome
########################################

if ! [ -f data/genome/mm39.fa ]; then
    wget -O - https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz |
    gzip -dc >data/genome/mm39.fa
fi

if ! [ -f data/genome/mm39.gtf ]; then
    wget -O - https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz |
    gzip -dc >data/genome/mm39.gtf
fi
