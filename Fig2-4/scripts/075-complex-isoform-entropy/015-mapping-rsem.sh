#!/bin/bash

###############################################################################
# Index
###############################################################################

mkdir -p data/genome/rsem_index

[ -f data/genome/rsem_index/SAindex ] || time rsem-prepare-reference --gtf data/genome/mm39.gtf --star -p 12 data/genome/mm39.fa data/genome/rsem_index/mm39 # real    50m37.548s

###############################################################################
# BAM
###############################################################################
# rm -rf data/rsem/bam
mkdir -p data/rsem/bam
num_threads=12

find data/fastq_trimmed -type f |
    grep fq.gz$ |
    sort |
    paste - - |
    sort -u |
    while read -r R1 R2; do
        filename=$(basename "${R1%_R1_*}" | cut -d "_" -f 1-3)

        if [ -f data/rsem/bam/"$filename".isoforms.results ]; then
            echo "$filename" already processed. Skipping...
            continue
        fi

        echo "======================================"
        echo "$filename" is now processing...
        echo "======================================"

        if ! [ -f tmp_R1_"$filename".fq ]; then
            echo "$filename" is unzipping...
            zcat "$R1" >tmp_R1_"$filename".fq &
            zcat "$R2" >tmp_R2_"$filename".fq &
            time wait # 6 minutes
        fi

        time rsem-calculate-expression \
            --star \
            --output-genome-bam \
            --sort-bam-by-coordinate \
            --paired-end \
            -p "$num_threads" \
            --append-names \
            tmp_R1_"$filename".fq \
            tmp_R2_"$filename".fq \
            data/genome/rsem_index/mm39 \
            data/rsem/bam/"$filename"
        # 6 hour / sample...
    done

# rm tmp_R1.fq tmp_R2.fq
## STAR version: 2.7.11b

# find data/fastq_trimmed -type f |
#     grep fq.gz$ |
#     sort |
#     paste - - |
#     sort -u |
#     while read -r R1 R2; do
#         filename=$(basename "${R1%_R1_*}" | cut -d "_" -f 1-3)
#         # filename="Trim71_KO_2"
#         if [ -f data/rsem/bam/"$filename".temp/"$filename".bam ]; then
#             echo "$filename".bam already exists. Skipping...
#             continue
#         fi

#         time STAR --genomeDir data/genome/rsem_index --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN "$num_threads" --genomeLoad NoSharedMemory --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outSAMheaderHD @HD VN:1.4 SO:unsorted --outFileNamePrefix data/rsem/bam/"$filename".temp/"$filename" --readFilesIn tmp_R1_"$filename".fq tmp_R2_"$filename".fq
#     done

# find data/rsem/bam -type f |
#     grep isoforms.results$ |
#     sed "s|data/rsem/bam/||" |
#     sed "s|\.isoforms.results||" |
#     sort >tmp_finished_samples.txt

# find data/fastq_trimmed -type f |
#     grep fq.gz$ |
#     sort |
#     paste - - |
#     cut -f 1 |
#     sed "s|data/fastq_trimmed/||" |
#     sed "s|_1.fq.gz||" |
#     grep -v -f tmp_finished_samples.txt |
#     paste - - - - - |
#     while read -r sample1 sample2 sample3 sample4 sample5; do
#         rsem-parse-alignments data/genome/rsem_index/mm39 data/rsem/bam/"$sample1".temp/"$sample1" data/rsem/bam/"$sample1".stat/"$sample1" data/rsem/bam/"$sample1".temp/"$sample1".bam 3 -tag XM &
#         rsem-parse-alignments data/genome/rsem_index/mm39 data/rsem/bam/"$sample2".temp/"$sample2" data/rsem/bam/"$sample2".stat/"$sample2" data/rsem/bam/"$sample2".temp/"$sample2".bam 3 -tag XM &
#         rsem-parse-alignments data/genome/rsem_index/mm39 data/rsem/bam/"$sample3".temp/"$sample3" data/rsem/bam/"$sample3".stat/"$sample3" data/rsem/bam/"$sample3".temp/"$sample3".bam 3 -tag XM &
#         rsem-parse-alignments data/genome/rsem_index/mm39 data/rsem/bam/"$sample4".temp/"$sample4" data/rsem/bam/"$sample4".stat/"$sample4" data/rsem/bam/"$sample4".temp/"$sample4".bam 3 -tag XM &
#         rsem-parse-alignments data/genome/rsem_index/mm39 data/rsem/bam/"$sample5".temp/"$sample5" data/rsem/bam/"$sample5".stat/"$sample5" data/rsem/bam/"$sample5".temp/"$sample5".bam 3 -tag XM &
#         time wait
#     done
