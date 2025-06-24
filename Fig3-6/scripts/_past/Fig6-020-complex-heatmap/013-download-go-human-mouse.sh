#!/bin/bash

mkdir -p data/Fig6

echo "taxon_id,symbol,go" > data/Fig6/go_human_mouse.csv
wget -qO- https://current.geneontology.org/annotations/goa_human.gaf.gz |
    gzip -d | grep TRBV20OR9-2 | head
    grep -v "^!" |
    cut -f 3,5 |
    tr "\t" "," |
    sed "s|^|9606,|" >> data/Fig6/go_human_mouse.csv

wget -qO- https://current.geneontology.org/annotations/goa_mouse.gaf.gz |
    gzip -d |
    grep -v "^!" |
    cut -f 3,5 |
    tr "\t" "," |
    sed "s|^|10090,|" >> data/Fig6/go_human_mouse.csv
