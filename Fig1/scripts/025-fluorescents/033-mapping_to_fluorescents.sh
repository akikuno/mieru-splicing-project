#!/bin/bash

mkdir -p data/genome/star_index data/bam_fluorescents

STAR --version | sed "s|^|STAR,|" >>reports/versions.csv
samtools --version | head -n 1 | tr " " "," >>reports/versions.csv

threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
threads=$((threads - 4))
[ "$threads" -lt 1 ] && threads=1

if ! [ -f data/genome/star_index_fluorescents/genomeParameters.txt ]; then
    time STAR \
    --runThreadN "$threads" \
    --runMode genomeGenerate \
    --genomeDir data/genome/star_index_fluorescents \
    --genomeFastaFiles data/fluorescents.fasta
fi

find data/fastq_trimmed -type f |
grep fq.gz$ |
sort |
paste - - |
sort -u |
while read -r R1 R2; do
    filename=$(basename "${R1%_1.fq.gz*}" | cut -d "_" -f 1-3)

    if [ -s data/bam_fluorescents/"$filename".bam ]; then
        echo "$filename".bam already exists. Skipping...
        continue
    fi

    echo "======================================"
    echo "$filename" is now processing...
    echo "======================================"

    gzip -dc "$R1" >tmp_R1.fq
    gzip -dc "$R2" >tmp_R2.fq
    time STAR \
    --runThreadN "$threads" \
    --genomeDir data/genome/star_index_fluorescents \
    --readFilesIn tmp_R1.fq tmp_R2.fq \
    --outSAMtype BAM SortedByCoordinate
    mv Aligned.sortedByCoord.out.bam data/bam_fluorescents/"$filename".bam
    time samtools index -@ "$threads" data/bam_fluorescents/"$filename".bam
    cat Log.final.out >data/bam_fluorescents/"$filename"_log.txt
    rm tmp* Log* SJ.out.tab
done

