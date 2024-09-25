#!/bin/bash

mkdir -p data/counts/ reports/counts

threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
threads=$((threads - 4))
[ "$threads" -lt 1 ] && threads=1

##################################
# featureCounts
##################################

# [ -s data/genome/mm39_transcript_name.gtf ] || grep transcript_name data/genome/mm39.gtf >data/genome/mm39_transcript_name.gtf
[ -s data/genome/mm39_gene_name.gtf ] || grep gene_name data/genome/mm39.gtf >data/genome/mm39_gene_name.gtf

time featureCounts -T "$threads" -p -B \
    -t exon -g gene_name -a data/genome/mm39_gene_name.gtf \
    -o data/counts/counts_gene_name.txt data/bam/*.bam

cat data/counts/counts_gene_name.txt |
    sed 1d |
    cut -f 1,6- |
    sed "s|data/bam/||g" |
    sed "s|.bam||g" |
    gzip -c >data/counts/featurecounts_gene_name.tsv.gz
