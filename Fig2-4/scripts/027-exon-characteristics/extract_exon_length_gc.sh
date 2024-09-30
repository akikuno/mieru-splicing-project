#!/bin/bash

rm -rf data/rmats/rmats_jc
mkdir -p data/rmats/rmats_jc data/exon_conservations

find original_output |
    grep -v tmp |
    grep "JC.txt$" |
    while read -r file; do
        echo $file
        ko_symbol=$(echo $file | cut -d "/" -f 2)
        file_basename=$(basename "$file")
        output="$ko_symbol"_"$file_basename"
        cp "$file" data/rmats/rmats_jc/"$output"
    done

#########################
# Extract SE events with FDR < 0.05 and abs(dPSI) > 0.1
#########################

: > tmp_se.bed
find data/rmats/rmats_jc -type f |
    grep SE |
    while read -r file; do
    ko_symbol=$(basename $file | cut -d "_" -f 1)
    echo $ko_symbol
    cat "$file" |
        awk '$20 < 0.05 && ($NF > 0.1 || $NF < -0.1)' |
        awk '{print $4,$6,$7,$3,$7-$6, $20, $NF}' |
        sed "s|$| ${ko_symbol}|" |
        tr -d '"' |
        tr " " "\t" >> tmp_se.bed
    done

sed "s|chr||" tmp_se.bed | sort -k 1,1 -k 2,2n > tmp_se_sorted.bed

#########################
# get fasta
#########################

time bedtools getfasta -fi data/genome/mm39.fa -bed tmp_se_sorted.bed -fo tmp_se.fa # 15 sec

sed "s|>||" tmp_se.fa | paste - - | sort -k 1,1 > tmp_se_sorted.txt

awk '{print $1":"$2"-"$3,$4,$5,$6,$7,$8}' tmp_se_sorted.bed | sort > tmp_se_sorted_regions.txt

join -1 1 -2 1 tmp_se_sorted_regions.txt tmp_se_sorted.txt | sort -u | tr " " "\t" > tmp_se_sorted_final.txt

#########################
# Calculate %GC
#########################

cut -f 7 tmp_se_sorted_final.txt |
    awk '{print (gsub(/[C|G]/,"&") / length($0)) * 100}' |
    paste tmp_se_sorted_final.txt - > tmp_se_sorted_final_gc.txt


#########################
# Output raw data
#########################
echo "ko_symbol,target_symbol,fdr,dpsi,mm39_coordinate,exon_length,percent_gc,exon_seq" > data/exon_conservations/exon_length_gc_all.csv
awk '{print $6, $2, $4, $5, $1, $3, $8, $7}' tmp_se_sorted_final_gc.txt |
    tr " " "," |
    sort >> data/exon_conservations/exon_length_gc_all.csv

#########################
# Summarize data
#########################

echo "ko_symbol, exon_length, gc" > data/exon_conservations/exon_length_gc.csv
awk '{print $6,$3,$NF}' tmp_se_sorted_final_gc.txt |
    tr " " "," |
    sort >> data/exon_conservations/exon_length_gc.csv

rm tmp_*
