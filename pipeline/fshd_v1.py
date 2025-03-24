import os
import argparse
import sys
import pandas as pd
from plotnine import *

### Basic infos ###
# Get arguments
parser = argparse.ArgumentParser(description='FSHD Pipeline - Analyze FSHD D4Z4 counts and methylation level')
parser.add_argument('-p', '--prefix', required=True, help='Sample Prefix')
PREFIX = parser.parse_args().prefix

# working directory
WD = f'/lab_nas_home_user/Mark/ONT/{PREFIX}'
os.makedirs(WD, exist_ok=True)

# Set threads
CORES = int(round(os.cpu_count() * 0.8))

# Softwares
SAMTOOLS = '/Everythings/misc/samtools-1.21/samtools'
MAKEBLASTDB = "/Everythings/misc/ncbi-blast-2.13.0+/bin/makeblastdb"
BLASTN = "/Everythings/misc/ncbi-blast-2.13.0+/bin/blastn"
BAM2CPG = '/Everythings/misc/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores'

# Path
TMP = '/raid/tmp/FSHD'
os.makedirs(TMP, exist_ok=True)
REF = '/nanopore_nas/pipeline/data/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set_modified.fasta'
QUERY = "/lab_nas_home_user/Mark/FSHD/querySeq"

# Other infos
blastTitle = ['QuerySeq-id', 'Seq-id', 'Identity%', 'Alignment-length', 'Mismatches', 'Gap', 'Q.Start', 'Q.End', 'S.Start', 'S.End', 'E-value', 'Bit-score']



### Merge BAM files ###
BAMS = f'/nanopore_nas/nanopore/{PREFIX}'
os.system(f'find {BAMS} -type f -name "*.bam" > {TMP}/bams.txt')
os.system(f'{SAMTOOLS} merge -b {TMP}/bams.txt -R chr4:190000000-190175845 -o - -@ {CORES} | {SAMTOOLS} sort -m 10G -u -@ {CORES} | {SAMTOOLS} view -b -@ {CORES} -o {TMP}/tmp.bam')
os.system(f'{SAMTOOLS} index -@ {CORES} {TMP}/tmp.bam')

### Create fasta file from BAM file
os.system(f'{SAMTOOLS} fasta -@ {CORES} {TMP}/tmp.bam > {TMP}/tmp_chr4.fasta')

### Create BLAST database
DB = f'{TMP}/blastdb'
os.system(f'{MAKEBLASTDB} -in {TMP}/tmp_chr4.fasta -dbtype nucl -out {DB} -title {PREFIX}')


### Sequence ID with tail sequence (pLAM)
os.system(f'{BLASTN} -db {DB} -query {QUERY}/pLAM.fasta -out {TMP}/tail.blast -num_threads {CORES} -outfmt 6')

tail = (
    pd.read_table(f'{TMP}/tail.blast', header=None, sep='\t', names=blastTitle)
        .drop(columns=['QuerySeq-id'])
        .sort_values(by='Identity%', ascending=False)
)


### Sequence ID with head sequence
os.system(f'{BLASTN} -db {DB} -query {QUERY}/P13E11.fasta -out {TMP}/head.blast -num_threads {CORES} -outfmt 6')

head = (
    pd.read_table(f'{TMP}/head.blast', header=None, sep='\t', names=blastTitle)
        .drop(columns=['QuerySeq-id'])
        .sort_values(by='Identity%', ascending=False)
)


### Count D4Z4 counts of all available reads
os.system(f'{BLASTN} -db {DB} -query {QUERY}/MIR8078.fasta -out {TMP}/d4z4.blast -num_threads {CORES} -outfmt 6')

D4Z4 = (
    pd.read_table(f'{TMP}/d4z4.blast', header=None, sep='\t', names=blastTitle)
        .drop(columns=['QuerySeq-id'])
        .sort_values(by='Identity%', ascending=False)
)

D4Z4 = D4Z4.groupby('Seq-id').size().reset_index(name='Count').sort_values(by='Count', ascending=False)
D4Z4.head()

# Add a new column "Target" to the D4Z4 dataframe
D4Z4['Target'] = D4Z4['Seq-id'].apply(
    lambda x:   'both' if x in head['Seq-id'].values and x in tail['Seq-id'].values else
                'head' if x in head['Seq-id'].values else
                'tail' if x in tail['Seq-id'].values else
                'none'
)


### Check if permissive (poly-A)
# Sequence ID with poly-A sequence
os.system(f'{BLASTN} -db {DB} -query {QUERY}/pLAM-polyA.fasta -out {TMP}/polyA.blast -num_threads {CORES} -outfmt 6')
polyA = (
    pd.read_table(f'{TMP}/polyA.blast', header=None, sep="\t", names=blastTitle)
        .sort_values(by='Identity%', ascending=False)
).iloc[:, [1, 2, 4, 3, 0]]

# Add permissive information to D4Z4 dataframe
def determine_permissive(seq):
    tmp = polyA[polyA['Seq-id'] == seq]
    if tmp.shape[0] == 2:
        diff = tmp.loc[tmp['QuerySeq-id'] == 'hapA', 'Mismatches'].values[0] - tmp.loc[tmp['QuerySeq-id'] == 'hapB', 'Mismatches'].values[0]
        if diff == 1:
            return 'Permissive (4qA)'
        elif diff == -1:
            return 'Non-permissive (4qB)'
    return 'Unknown'

D4Z4['Permissive'] = D4Z4['Seq-id'].apply(determine_permissive)


### Save results
D4Z4.to_csv(f'{WD}/FSHD-{PREFIX}.tsv', sep='\t', index=False)

# BAM file
pd.DataFrame({'Seq-id': list(set(head['Seq-id']).union(set(tail['Seq-id'])))}).to_csv(f'{TMP}/headandtail.txt', index=False, header=False)
os.system(f'{SAMTOOLS} view -C -T {REF} -@ {CORES} -N {TMP}/headandtail.txt -o {WD}/DUX4.cram {TMP}/tmp.bam')
os.system(f'{SAMTOOLS} index -@ {CORES} {WD}/DUX4.cram')


### Calculate methylation level
os.system(f'{BAM2CPG} --bam {WD}/DUX4.cram --ref {REF} --model /nanopore_nas/pipeline/data/pileup_calling_model.v1.tflite --output-prefix {TMP}/tmp --threads {CORES}')

# Load CpG scores
meth = (
    pd.read_table(f'{TMP}/tmp.combined.bed', header=None, sep='\t', 
                    usecols=[0, 2, 3, 5], names=['chrom', 'pos', 'meth', 'coverage'])
        .query('coverage > 1 and pos >= 190175600 and pos <= 190175834')
)

# Plot
plt = (
    ggplot(meth, aes(y='meth')) +
        geom_boxplot() +
        labs(y="Methylation frequency (%)") +
        ylim(0, 100) +
        theme_bw() + 
        theme(figure_size=(3, 3), dpi=200,  axis_ticks_x=None, axis_text_x=None, axis_title_x=None)
)
plt.save(f'{WD}/methylation.jpeg')

# Remove temporary files
os.system(f'rm -rf {TMP}/*')
sys.exit(0)