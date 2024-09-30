#!/bin/bash

mkdir -p data/exon_conservations

#######################################
# Donwlonad phylop scores for mm39
#######################################

wget -P data/exon_conservations https://hgdownload.cse.ucsc.edu/goldenPath/mm39/phyloP35way/mm39.phyloP35way.bw

#######################################
# Convert to bedgraph
#######################################

[ -f data/exon_conservations/bigWigToBedGraph ] || wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod 755 data/exon_conservations/bigWigToBedGraph

time ./data/exon_conservations/bigWigToBedGraph \
    ./data/exon_conservations/mm39.phyloP35way.bw \
    ./data/exon_conservations/mm39.phyloP35way.bedgraph

#######################################
# Split by Chromosome
#######################################

for chr in {1..19} X Y; do
    echo "$chr" is processing...
    time grep chr"$chr" data/exon_conservations/mm39.phyloP35way.bedgraph >data/exon_conservations/mm39.phyloP35way.chr"$chr".bedgraph
done

time awk '{ print > "data/exon_conservations/mm39.phyloP35way."$1".bedgraph" }' data/exon_conservations/mm39.phyloP35way.bedgraph

#######################################
# Subset Skipped Exon regions
#######################################

find data/rmats/original_output -type f |
    grep JC.txt |
    grep SE |
    while read -r file; do
        ko_symbol=$(echo "$file" | cut -d/ -f 4)

        cat "$file" |
            awk '{print $3, $4, $6, $7, $20, $NF}' |
            # Filter by FDR < 0.05 and IncLevelDifference > 0.1 or IncLevelDifference < -0.1
            awk '$(NF-1) < 0.05 && ($NF > 0.1 || $NF < -0.1)' |
            tr -d '"' |
            awk '{print $2,$3,$4,$1,$5,$6}' |
            tr " " "\t" |
            sort -k1,1 -k2,2n >data/exon_conservations/tmp_SE.bed

        : > data/exon_conservations/"$ko_symbol".bed
        for chr in {1..19} X Y; do
            echo "$ko_symbol" of chr"$chr" is intersecting...
            time bedtools intersect \
                -a data/exon_conservations/tmp_SE.bed \
                -b data/exon_conservations/mm39.phyloP35way.chr"$chr".bedgraph \
                -wa -wb |
            awk '{print $1,$2,$3,$4,$5,$6,$NF}' |
            tr " " "\t" >> data/exon_conservations/"$ko_symbol".bed
        done
    done

rm data/exon_conservations/tmp_SE.bed

#######################################
# Merge files
#######################################

echo "ko_symbol,chr,start,end,target_symbol,fdr,dpsi,phylop" > data/exon_conservations/exon_conservation.csv

find data/exon_conservations/ -type f |
    grep ".bed$" |
    while read -r file; do
        ko_symbol=$(basename "$file" | cut -d. -f1)
        echo "$ko_symbol" is merging...
        sed "s|^|${ko_symbol},|" "$file" | tr "\t" "," >> data/exon_conservations/exon_conservation.csv
    done
