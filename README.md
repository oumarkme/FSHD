# FSHD

## Usage:
```{bash}
python pipeline/fshd_v2.py -p PREFIX
```

## Required python packages
- os
- argparse
- pandas
- plotnine

## Required system installed softwares
- samtools
- makeblasdb
- blastn
- pb-CpG-tools

## Protocol
- From BAM files, extract reads mapped to chr4:190000000-190175845

- Create FASTA file of these reads

- Create BLAST database using `makeblastdb`

- Grab sequence ID with tail sequence
    - Exon 2 (version 2)
    - pLAM (version 1)

- Grab sequence ID with head sequence
    - DBET (version 2)
    - P13E11 (version 1)

- Count D4Z4 count of all available reads
    - Count by MIR8078 sequence

- Check if permissive (poly-A tail)

- Methylation level
    - Extract methylation level using `aligned_bam_to_cpg_scores` from `pb-CpG-tools`
    - Plot by region `chr4:190175600-190175834`