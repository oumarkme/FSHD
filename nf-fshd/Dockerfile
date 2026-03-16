FROM condaforge/miniforge3:latest

LABEL maintainer="Mark"
LABEL description="FSHD D4Z4 analysis pipeline (v3) for ONT long-read data"

RUN mamba install -y -n base \
    bioconda::samtools=1.23 \
    bioconda::bcftools=1.23 \
    bioconda::minimap2=2.30 \
    bioconda::blast=2.17.0 \
    bioconda::ont-modkit=0.6.1 \
    bioconda::pysam=0.23.3 \
    conda-forge::pandas=3.0.1 \
    conda-forge::plotnine=0.15.3 \
    && mamba clean -afy

WORKDIR /data
