# mieru-splicing-project
 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mieru -y python=3.10
conda install -n mieru -y \
    fastqc fastp \
    star samtools bedtools deeptools subread rmats \
    r-base r-essentials r-writexl r-janitor r-extrafont \
    r-ggfortify r-ggrepel r-patchwork r-enrichr \
    bioconductor-deseq2 \
    bioconductor-goseq \
    bioconductor-org.Mm.eg.db \
    bioconductor-TxDb.Mmusculus.UCSC.mm10.ensGene \

conda activate mieru
