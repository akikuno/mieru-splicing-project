#!/bin/bash

fastp --version 2>&1 | tr " " "," >>reports/versions.csv

mkdir -p data/fastq_trimmed

threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
threads=$((threads - 4))
[ "$threads" -lt 1 ] && threads=1

find data/fastq -type f |
grep fq.gz |
sort |
paste - - |
while read -r R1 R2; do
    [ -s data/fastq_trimmed/"$(basename "$R1")" ] && continue
    basename "$R1"
    basename "$R2"
    echo "=============="
    filename_html=$(basename "$R1" | cut -d "_" -f 1-3).html
    time fastp \
    -i "$R1" -I "$R2" \
    -o data/fastq_trimmed/"$(basename "$R1")" \
    -O data/fastq_trimmed/"$(basename "$R2")" \
    -h data/fastq_trimmed/"$filename_html" \
    --thread "$threads"
done
