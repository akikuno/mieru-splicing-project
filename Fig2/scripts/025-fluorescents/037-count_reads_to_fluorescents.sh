#!/bin/bash

mkdir -p data/counts/

: > data/counts/total_reads.txt
find data/fastq_trimmed/*_1.fq.gz |
    sort |
    while read -r fq; do
    sample="$(basename ${fq%_1.fq.gz})"
    echo $fq
    zcat $fq | grep -c "^+$" | sed "s|^|${sample}\t|" >> data/counts/total_reads.txt
done

: > data/counts/fluorescent_reads.txt
find data/bam_fluorescents/*.bam |
    sort |
    while read -r bam; do
    sample="$(basename ${bam%.bam})"
    echo $sample
    samtools view "$bam" | cut -f 1,3 | sort -u | grep -c "EGFP" | sed "s|^|${sample}\tEGFP\t|" >> data/counts/fluorescent_reads.txt
    samtools view "$bam" | cut -f 1,3 | sort -u | grep -c "tdTomato" | sed "s|^|${sample}\ttdTomato\t|" >> data/counts/fluorescent_reads.txt
    samtools view "$bam" | cut -f 1,3 | sort -u | grep -c "TagBFP" | sed "s|^|${sample}\tTagBFP\t|" >> data/counts/fluorescent_reads.txt
done

join -a 2 data/counts/total_reads.txt data/counts/fluorescent_reads.txt |
    awk '{print $0, $4/$2 * 100}' |
    tr " " "," > data/counts/read_counts_fluorescents.csv


# sort data/counts/total_reads.txt > tmp_1
# sort data/counts/fluorescent_reads.txt > tmp_2
# join -a 1 tmp_1 tmp_2 | awk '{print $0, $4/$2 * 100}' | tr " " "," > data/counts/read_counts_fluorescents.csv
