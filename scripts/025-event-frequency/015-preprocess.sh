#!/bin/bash

output_file="data/rmats/all_events_ko_target_fdr_dpsi.csv"
echo "event, ko_symbol, target_symbol, fdr, dpsi" > "$output_file"

find data/rmats/original_output -type f |
    grep "JC.txt" |
    while read -r file; do
    SYMBOL="$(echo "$file" | cut -d "/" -f 4)"
    EVENT="$(echo "$file" | cut -d "/" -f 5 | cut -d "." -f 1)"
    echo "$SYMBOL" "$EVENT"
    sed 1d "$file" |
        awk '{print $3,$(NF-3), $NF}' |
        tr " " "," |
        sed "s|^|${SYMBOL},|g" |
        sed "s|^|${EVENT},|g" >> "$output_file"
    done

