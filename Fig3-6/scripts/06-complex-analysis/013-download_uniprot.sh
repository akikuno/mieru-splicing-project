#!/bin/bash

mkdir -p data/Fig5

# time curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz > data/Fig5/uniprot_idmapping.txt.gz # 30min
# Uniprot ID
time curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz > data/Fig5/uniprot_idmapping_selected.txt.gz # 231m
time zcat data/Fig5/uniprot_idmapping_selected.txt.gz |
    cut -f 1,3,13 |
    awk '$3 == 9606 || $3 == 10090' | # Human and Mouse
    awk 'BEGIN{OFS="\t"} {print $3, $2, $1}' |
    sort > data/Fig5/uniprot_idmapping_selected_flatten.txt # 20 min

# NCBI Gene ID
time curl https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz > data/Fig5/All_Mammalia.gene_info.gz # 4 min
zcat data/Fig5/All_Mammalia.gene_info.gz | cut -f 1-3 | awk '$1 == 9606 || $1 == 10090' | sort > data/Fig5/gene_info_mouse_human.txt

awk '{print $1"-"$2 "\t" $3}' data/Fig5/uniprot_idmapping_selected_flatten.txt | sort > tmp_uniprot.txt

awk '{print $1"-"$2 "\t" $3}' data/Fig5/gene_info_mouse_human.txt | sort > tmp_gene_info.txt

wc -l tmp_uniprot.txt tmp_gene_info.txt

echo "taxon_id,uniprot,symbol" > data/Fig5/uniprot_gene_symbol_mouse_human.csv
join -t $'\t' tmp_uniprot.txt tmp_gene_info.txt | tr "-" "\t" | cut -f 1,3,4 | sort | tr "\t" "," >> data/Fig5/uniprot_gene_symbol_mouse_human.csv
