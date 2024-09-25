#!/bin/bash

mkdir -p data/rmats/original_output

threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
threads=$((threads - 4))
[ "$threads" -lt 1 ] && threads=1

find data/bam -type f |
    grep "bam$" |
    grep MIERU_CH |
    tr "\n" "," |
    sed "s|,$||" |
    grep ^ >tmp_rmats_control.txt

samples=$(
    find data/bam -type f |
        grep "bam$" |
        grep -v MIERU_CH |
        sed "s|data/bam/||" |
        sed "s|_KO.*||" |
        sort -u
)

for sample in $samples; do
    if [ -d data/rmats/original_output/"$sample" ]; then
        echo "$sample" already exists. Skipping...
        continue
    fi
    echo $sample
    find data/bam -type f |
        grep "bam$" |
        grep "$sample" |
        tr "\n" "," |
        sed "s|,$||" |
        grep ^ >tmp_rmats_sample.txt
    mkdir -p data/rmats/original_output/"$sample"/tmp
    time RNASeq-MATS.py \
        --b1 tmp_rmats_control.txt \
        --b2 tmp_rmats_sample.txt \
        --gtf data/genome/mm39.gtf \
        --od data/rmats/original_output/"$sample" \
        --tmp data/rmats/original_output/"$sample"/tmp \
        --variable-read-length --readLength 150 \
        --allow-clipping \
        -t paired \
        --nthread "$threads"
done

rm tmp_rmats_*
