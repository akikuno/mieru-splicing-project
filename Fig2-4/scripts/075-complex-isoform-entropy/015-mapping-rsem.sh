#!/bin/bash

###############################################################################
# Index
###############################################################################

mkdir -p data/genome/rsem_index

time rsem-prepare-reference --gtf data/genome/mm39.gtf --star -p 12 data/genome/mm39.fa data/genome/rsem_index/mm39
# real    50m37.548s

###############################################################################
# BAM
###############################################################################
# rm -rf data/rsem/bam
mkdir -p data/rsem/bam

find data/fastq_trimmed -type f |
grep fq.gz$ |
sort |
paste - - |
sort -u |
tail -n 12 |
while read -r R1 R2; do
    filename=$(basename "${R1%_R1_*}" | cut -d "_" -f 1-3)
    
    if [ -f data/rsem/bam/"$filename".isoforms.results ]; then
        echo "$filename" already processed. Skipping...
        continue
    fi
    
    echo "======================================"
    echo "$filename" is now processing...
    echo "======================================"
    
    zcat "$R1" >tmp_R1_"$filename".fq &
    zcat "$R2" >tmp_R2_"$filename".fq &
    time wait # 6 minutes
    
    time rsem-calculate-expression \
    --star \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --paired-end \
    -p 12 \
    --append-names \
    tmp_R1_"$filename".fq \
    tmp_R2_"$filename".fq \
    data/genome/rsem_index/mm39 \
    data/rsem/bam/"$filename"
    # 6 hour / sample...
done

rm tmp_R1.fq tmp_R2.fq
## STAR version: 2.7.11b

###############################################################################
# Expression
###############################################################################
mkdir -p data/rsem/expression

find data/rsem/bam -type f |
grep bam$ |
while read -r bam; do
    filename=$(basename "${bam%.bam}")
    echo "$filename"
    time rsem-calculate-expression \
    -p 12 \
    --paired-end \
    --alignments \
    --append-names \
    --estimate-rspd \
    --no-bam-output \
    SRR22571458_Aligned.toTranscriptome.out.bam \
    data/genome/rsem_index/mm39 \
    data/rsem/expression/"$filename"
done

###############################################################################
# Merge
###############################################################################
