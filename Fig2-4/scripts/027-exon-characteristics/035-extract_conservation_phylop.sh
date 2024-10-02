#!/bin/bash

mkdir -p data/exon_characteristics

#######################################
# Donwlonad phylop scores for mm39
#######################################

wget -P data/exon_characteristics https://hgdownload.cse.ucsc.edu/goldenPath/mm39/phyloP35way/mm39.phyloP35way.bw

#######################################
# Convert to bedgraph
#######################################

[ -f data/exon_characteristics/bigWigToBedGraph ] || wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod 755 data/exon_characteristics/bigWigToBedGraph

time ./data/exon_characteristics/bigWigToBedGraph \
    ./data/exon_characteristics/mm39.phyloP35way.bw \
    ./data/exon_characteristics/mm39.phyloP35way.bedgraph

#######################################
# Split by Chromosome to reduce memory usage
#######################################

mkdir -p data/exon_characteristics/phyloP35way_chr

for chr in {1..19} X Y; do
    echo "chr${chr}" is processing...
    time cat data/exon_characteristics/mm39.phyloP35way.bedgraph |
        awk -v chr=chr"$chr" '$1 == chr' >data/exon_characteristics/phyloP35way_chr/mm39.phyloP35way.chr"$chr".bedgraph
done

#######################################
# Subset Skipped Exon regions
#######################################

echo -e "ko_symbol\tphylop\texon_id" > data/exon_characteristics/exon_conservation.tsv

for chr in {1..19} X Y; do
    echo chr"$chr" is intersecting...
    awk -v chr="$chr" '$1 == chr' data/rmats/se_ko_symbol.bed | sed "s|^|chr|" > tmp_se.bed
    time bedtools intersect \
        -a tmp_se.bed \
        -b data/exon_characteristics/phyloP35way_chr/mm39.phyloP35way.chr"$chr".bedgraph \
        -wa -wb |
    awk '{print $4, $NF, $1"_"$2"_"$3"_"$5"_"$6}' |
    tr " " "\t" >> data/exon_characteristics/exon_conservation.tsv
done

rm tmp_se.bed
