#!/bin/bash

mkdir -p data/Fig5

time curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz > data/Fig5/uniprot_idmapping.txt.gz

time awk -F'\t' '
    $1 != prev {
        if (NR > 1) print line;
        prev = $1;
        line = $1 "\t" $2 ":" $3;
        next;
    }
    {
        line = line "\t" $2 ":" $3;
    }
    END { print line; }
    ' < data/Fig5/uniprot_idmapping.txt.gz |
    gzip -c > data/Fig5/uniprot_idmapping_flatten.txt.gz

zcat data/Fig5/uniprot_idmapping_flatten.txt.gz | wc -l
zcat data/Fig5/uniprot_idmapping_flatten.txt.gz |
    grep Gene_Name |
    grep -e NCBI_TaxID:9606 -e NCBI_TaxID:10090 |
    awk '{print $1; for (i=1; i<=NF; i++) if ($i ~ /^Gene_Name:/ || $i ~ /^NCBI_TaxID/) print $i}' |
    paste - - - |
    sed -e "s|Gene_Name:||" -e "s|NCBI_TaxID:||" > data/Fig5/uniprot_idmapping_mouse_human.txt

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
