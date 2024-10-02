#!/bin/bash

mkdir -p data/exon_characteristics

#########################
# get fasta
#########################

time bedtools getfasta -name -fi data/genome/mm39.fa -bed data/rmats/se_ko_symbol.bed -fo tmp_se.fa # 15 sec

sed "s|>||" tmp_se.fa | paste - - | sort -k 1,1 > tmp_se_sorted.txt

#########################
# Calculate Exon length and %GC
#########################

echo -e "ko_symbol\texon_length\tgc" > data/exon_characteristics/exon_length_gc.tsv

cat tmp_se_sorted.txt |
    awk '{
    sub("::.*$", "", $1)
    exon_length=length($NF)
    gc=gsub(/[C|G]/,"&", $NF) / exon_length * 100
    print $1, exon_length, gc}' |
    tr " " "\t" >> data/exon_characteristics/exon_length_gc.tsv

rm tmp_*
