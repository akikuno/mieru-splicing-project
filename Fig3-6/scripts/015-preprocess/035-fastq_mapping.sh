#!/bin/bash

mkdir -p data/genome/star_index data/bam data/bam_to_transcriptome

threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
threads=$((threads - 4))
[ "$threads" -lt 1 ] && threads=1

if ! [ -f data/genome/star_index/genomeParameters.txt ]; then
    time STAR \
        --runThreadN "$threads" \
        --runMode genomeGenerate \
        --genomeDir data/genome/star_index \
        --genomeFastaFiles data/genome/mm39.fa \
        --sjdbGTFfile data/genome/mm39.gtf
fi

find data/fastq_trimmed -type f |
    grep fq.gz$ |
    sort |
    paste - - |
    sort -u |
    while read -r R1 R2; do
        filename=$(basename "${R1%_R1_*}" | cut -d "_" -f 1-3)

        if [ -f data/bam/"$filename".bam ]; then
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
            --genomeDir data/genome/star_index \
            --readFilesIn tmp_R1.fq tmp_R2.fq \
            --quantMode TranscriptomeSAM \
            --outFilterMultimapNmax 1 \
            --outSAMtype BAM SortedByCoordinate
        mv Aligned.sortedByCoord.out.bam data/bam/"$filename".bam
        mv Aligned.toTranscriptome.out.bam data/bam_to_transcriptome/"$filename"_totranscriptome.bam
        time samtools index -@ "$threads" data/bam/"$filename".bam
        cat Log.final.out >data/bam/"$filename"_log.txt
        rm tmp* Log* SJ.out.tab
    done

: >reports/reads_mapped.csv

find data/bam -type f |
    grep bam$ |
    while read -r bam; do
        echo "$bam"
        sample=$(basename "${bam%.bam}")
        samtools view -F 4 "$bam" |
            cut -f 1 |
            sort -u |
            wc -l |
            sed "s|^|${sample},|" >> reports/reads_mapped.csv
    done

# samtools view data/bam/Ybx1_KO_7.bam | wc -l
# samtools view -F 4 data/bam/Ybx1_KO_7.bam | wc -l

###########################################################
# merge BAM files
###########################################################

mkdir -p data/bam_merged

samples=$(
    find data/bam -type f |
        grep bam$ |
        xargs -I{} basename {} |
        cut -d_ -f 1 |
        sort -u
)

for sample in $samples; do
    if ! [ -f data/bam_merged/"$sample".bam ]; then
        echo data/bam/"$sample"_*.bam
        time samtools merge -@ "$threads" data/bam_merged/"$sample".bam data/bam/"$sample"_*.bam
        time samtools index -@ "$threads" data/bam_merged/"$sample".bam
    fi
done
