FROM condaforge/miniforge3:latest

LABEL maintainer="Mark"
LABEL description="FSHD D4Z4 analysis pipeline (v3) for ONT long-read data"

RUN mamba install -y -n base \
    bioconda::samtools=1.21 \
    bioconda::bcftools=1.21 \
    bioconda::minimap2=2.28 \
    bioconda::winnowmap=2.03 \
    bioconda::meryl=1.4.1 \
    bioconda::blast=2.16.0 \
    bioconda::modkit=0.4.3 \
    bioconda::whatshap=2.3 \
    bioconda::sniffles=2.4 \
    bioconda::flye=2.9.5 \
    bioconda::last=1542 \
    bioconda::tandem-genotypes=1.9.2 \
    bioconda::ensembl-vep=113.0 \
    bioconda::pysam=0.22.1 \
    conda-forge::pandas=2.2.3 \
    conda-forge::plotnine=0.14.3 \
    && mamba clean -afy

WORKDIR /data
