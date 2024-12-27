#!/bin/bash

mkdir -p data/Fig5

curl https://current.geneontology.org/annotations/mgi.gaf.gz |
    gunzip > data/Fig5/mgi.gaf

cat data/Fig5/mgi.gaf |
    grep -v '^!' |
    cut -f 3,5 |
    sort -u > data/Fig5/mgi_symbol_go.txt
