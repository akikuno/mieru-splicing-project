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

# time zcat data/Fig5/uniprot_idmapping.txt.gz |
#     awk -F '\t' '
#     $1 != prev {
#         if (NR > 1) print line;
#         prev = $1;
#         line = $1 "\t" $2 ":" $3;
#         next;
#     }
#     {
#         line = line "\t" $2 ":" $3;
#     }
#     END { print line; }
#     ' |
#     gzip -c > data/Fig5/uniprot_idmapping_flatten.txt.gz # 5m14.102s

# zcat data/Fig5/uniprot_idmapping_flatten.txt.gz | wc -l
# time zcat data/Fig5/uniprot_idmapping_flatten.txt.gz |
#     grep Gene_Name |
#     grep -e NCBI_TaxID:9606 -e NCBI_TaxID:10090 |
#     awk '{print $1; for (i=1; i<=NF; i++) if ($i ~ /^Gene_Name:/ || $i ~ /^NCBI_TaxID/) print $i}' |
#     paste - - - |
#     sed -e "s|Gene_Name:||" -e "s|NCBI_TaxID:||" > data/Fig5/uniprot_idmapping_mouse_human.txt

# wc -l data/Fig5/uniprot_idmapping_mouse_human.txt
# ls -lh data/Fig5/uniprot_idmapping_mouse_human.txt

# zcat data/Fig5/uniprot_idmapping_flatten.txt.gz | grep -e Q9UK58 -e Q9UQ88 | tee tmp_uniprot.txt
# curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.example |
# awk -F'\t' '
#     $1 != prev {
#         if (NR > 1) print line;
#         prev = $1;
#         line = $1 "\t" $2 ":" $3;
#         next;
#     }
#     {
#         line = line "\t" $2 ":" $3;
#     }
#     END { print line; }
#     ' |
#     gzip -c > tmp.txt.gz

# zcat tmp.txt.gz | wc -l
# zcat tmp.txt.gz |
#     grep Gene_Name |
#     grep -e NCBI_TaxID:9606 -e NCBI_TaxID:10090 |
#     awk '{print $1; for (i=1; i<=NF; i++) if ($i ~ /^Gene_Name:/ || $i ~ /^NCBI_TaxID/) print $i}' |
#     paste - - - |
#     sed -e "s|Gene_Name:||" -e "s|NCBI_TaxID:||" > data/Fig5/uniprot_idmapping_mouse_human.txt

# cat << EOF |
# A	ID	hoge
# A	Num	3
# B	ID	fuga
# B	Num	4
# B	Other	XXX
# EOF
# awk -F'\t' '
# $1 != prev {
#     if (NR > 1) print line;
#     prev = $1;
#     line = $1 "\t" $2 ":" $3;
#     next;
# }
# {
#     line = line "\t" $2 ":" $3;
# }
# END { print line; }
# ' |
# gzip -c > data/Fig5/uniprot_idmapping_flatten.txt.gz
