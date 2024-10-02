#!/bin/bash

tmp_file="data/rmats/tmp_se.bed"
output_file="data/rmats/se_ko_symbol.bed"

find data/rmats/original_output -type f |
    grep "SE.MATS.JC.txt" |
    while read -r file; do
    SYMBOL="$(echo "$file" | cut -d "/" -f 4)"
    echo "$SYMBOL"
    sed 1d "$file" |
        # Filter by FDR < 0.05 and IncLevelDifference > 0.1 or IncLevelDifference < -0.1
        awk '$20 < 0.05 && ($NF > 0.1 || $NF < -0.1)' |
        awk -v symbol="$SYMBOL" '{print $4,$6,$7,symbol,$NF,$5}' |
        tr " " "\t" |
        sed "s|^chr||g" |
        sort -u >> "$tmp_file"
    done

sort -k1,1 -k2,2n "$tmp_file" > "$output_file"

rm "$tmp_file"
