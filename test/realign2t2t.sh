#!/bin/bash
PREFIX='WGS-ONT-A0002-7298388-D26478'

export PATH=/Everythings/misc/samtools-1.21:$PATH

BAMS='/nanopore_nas/nanopore/'${PREFIX}
TMP='/raid/tmp/playground'

# List all the bam files
find ${BAMS} -name '*.bam' > ${TMP}/bam.list

# Create Fastq file
REF='/lab_nas_home_user/Mark/reference_genome/hg38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set_modified.fasta'
samtools view -b -T ${REF} -@ 50 /lab_nas_home_user/Mark/ONT/${PREFIX}/DUX4.cram | \
samtools fastq -@ 50 - | \
gzip > ${TMP}/${PREFIX}.fastq.gz

# Align to reference
REF='/lab_nas_home_user/Mark/reference_genome/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna'
MMI='/lab_nas_home_user/Mark/reference_genome/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic.mmi'

sudo docker run --rm --gpus all \
    -v '/lab_nas_home_user':'/lab_nas_home_user' \
    -v '/raid':'/raid' \
    nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1 pbrun minimap2 \
    --ref ${REF} --index ${MMI} \
    --interval chr4:193000000-193574945 \
    --num-threads 50 \
    --preset map-ont \
    --gpuwrite \
    --gpusort \
    --tmp-dir ${TMP} \
    --in-fq ${TMP}/${PREFIX}.fastq.gz \
    --out-bam /lab_nas_home_user/Mark/ONT/${PREFIX}/DUX4.t2t.cram

