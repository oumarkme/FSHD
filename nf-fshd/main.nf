#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process EXTRACT_FSHD_READS {

    input:
    path input_file
    path input_index
    path ref_hg38

    output:
    path 'subset.bam'

    script:
    """
    cat > regions.bed << 'EOF'
chr4	189900000	190215000	4q35_D4Z4
chr10	133500000	133798000	10q26_D4Z4
EOF

    REF_OPT=""
    if [[ "${input_file}" == *.cram ]]; then
        REF_OPT="--reference ${ref_hg38}"
    fi

    # Tolerate CRAM decode errors on edge slices (common with ONT long reads)
    samtools view -b -h -M -L regions.bed \$REF_OPT -@ ${task.cpus} ${input_file} -o target.bam 2>&1 | \
        grep -v 'cram_decode_seq\\|cram_next_slice' >&2 || true
    samtools view -b -f 4 \$REF_OPT -@ ${task.cpus} ${input_file} -o unmapped.bam 2>&1 | \
        grep -v 'cram_decode_seq\\|cram_next_slice' >&2 || true
    samtools merge -@ ${task.cpus} -o subset.bam target.bam unmapped.bam
    """
}


process REALIGN_T2T {

    input:
    path subset_bam
    path ref_t2t
    val  sample_id

    output:
    path "${sample_id}.t2t.bam",     emit: bam
    path "${sample_id}.t2t.bam.bai", emit: bai

    script:
    """
    samtools fastq -T MM,ML -@ ${task.cpus} ${subset_bam} \
        | minimap2 -ax lr:hq -y --eqx --secondary=no -t ${task.cpus} \
            ${ref_t2t} - \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.t2t.bam

    samtools index -@ ${task.cpus} ${sample_id}.t2t.bam
    """
}


process EXTRACT_D4Z4 {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    path bam
    path bai
    val  sample_id

    output:
    path "${sample_id}.d4z4.bam",     emit: bam
    path "${sample_id}.d4z4.bam.bai", emit: bai

    script:
    """
    cat > d4z4_regions.bed << 'EOF'
chr4	193230000	193610000	chr4_D4Z4
chr10	134560000	134800000	chr10_D4Z4
EOF

    samtools view -b -L d4z4_regions.bed -@ ${task.cpus} ${bam} -o ${sample_id}.d4z4.bam
    samtools index -@ ${task.cpus} ${sample_id}.d4z4.bam
    """
}


process BLAST_ANALYSIS {

    input:
    path bam
    path bai
    val  sample_id
    path query_d4z4
    path query_plam_polya
    path query_p13e11

    output:
    path "${sample_id}.d4z4_counts.tsv"

    script:
    """
    samtools fasta -@ ${task.cpus} ${bam} > reads.fasta
    makeblastdb -in reads.fasta -dbtype nucl -out blastdb
    blastn -db blastdb -query ${query_d4z4} -out d4z4.blast -num_threads ${task.cpus} -outfmt 6
    blastn -db blastdb -query ${query_plam_polya} -out polya.blast -num_threads ${task.cpus} -outfmt 6
    blastn -db blastdb -query ${query_p13e11} -out p13e11.blast -num_threads ${task.cpus} -outfmt 6

    python3 - "${sample_id}" << 'PYEOF'
import sys
import pysam
import pandas as pd

sample_id = sys.argv[1]

blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# --- Chromosome assignment from BAM alignments ---
read_chr = {}
with pysam.AlignmentFile('${bam}', 'rb') as af:
    for read in af:
        if not read.is_unmapped and read.query_name:
            read_chr[read.query_name] = read.reference_name

# --- D4Z4 RU counting (filter alignment length > 1500) ---
d4z4 = pd.read_csv('d4z4.blast', sep='\\t', header=None, names=blast_cols)
d4z4 = d4z4[d4z4['length'] > 1500]
ru_counts = d4z4.groupby('sseqid').size().reset_index(name='RU_Count')
ru_counts = ru_counts.sort_values('RU_Count', ascending=False)

# --- Completeness: proximal (P13E11) + distal (pLAM) per DUCKS4 ---
proximal_reads = set()  # P13E11 = D4F104S1 probe (proximal anchor)
distal_reads = set()    # pLAM (distal anchor)
polya = pd.DataFrame(columns=blast_cols)
try:
    p13e11_df = pd.read_csv('p13e11.blast', sep='\\t', header=None, names=blast_cols)
    proximal_reads = set(p13e11_df['sseqid'].unique())
except pd.errors.EmptyDataError:
    pass
try:
    polya = pd.read_csv('polya.blast', sep='\\t', header=None, names=blast_cols)
    distal_reads = set(polya['sseqid'].unique())
except pd.errors.EmptyDataError:
    pass

def classify_completeness(read_id):
    has_proximal = read_id in proximal_reads
    has_distal = read_id in distal_reads
    if has_proximal and has_distal:
        return 'complete'
    elif has_distal:
        return 'partial-distal'
    elif has_proximal:
        return 'partial-proximal'
    return 'incomplete'

# --- Permissive check (pLAM-polyA, chr4 only) ---
# hapA = ATTAAA (4qA, permissive), hapB = ATCAAA (4qB/10q)
# diff = hapA.mismatch - hapB.mismatch:
#   < 0 → read matches hapA better → Permissive (4qA)
#   > 0 → read matches hapB better → Non-permissive (4qB)
permissive_map = {}
for seq_id, grp in polya.groupby('sseqid'):
    if read_chr.get(seq_id) != 'chr4':
        continue
    hap_a = grp[grp['qseqid'] == 'hapA']
    hap_b = grp[grp['qseqid'] == 'hapB']
    if len(hap_a) >= 1 and len(hap_b) >= 1:
        diff = hap_a.iloc[0]['mismatch'] - hap_b.iloc[0]['mismatch']
        if diff < 0:
            permissive_map[seq_id] = 'Permissive (4qA)'
        elif diff > 0:
            permissive_map[seq_id] = 'Non-permissive (4qB)'
        else:
            permissive_map[seq_id] = 'Unknown'

# --- Build output ---
ru_counts['Chromosome'] = ru_counts['sseqid'].map(read_chr).fillna('Unknown')
ru_counts['Completeness'] = ru_counts['sseqid'].apply(classify_completeness)
ru_counts['Permissive'] = ru_counts['sseqid'].map(permissive_map).fillna('-')

result = ru_counts.rename(columns={'sseqid': 'Read_ID'})[
    ['Read_ID', 'Chromosome', 'RU_Count', 'Completeness', 'Permissive']
]
result.to_csv(f'{sample_id}.d4z4_counts.tsv', sep='\\t', index=False)
PYEOF
    """
}


process METHYLATION {

    input:
    path bam
    path bai
    path counts_tsv
    path ref_t2t
    val  sample_id

    output:
    path "${sample_id}.methylation.bed"

    script:
    """
    cat > meth_targets.bed << 'EOF'
chr4	193540172	193543634
EOF

    # Build read selection tiers: most specific → broadest fallback
    python3 - << 'PYEOF'
import pandas as pd

counts = pd.read_csv('${counts_tsv}', sep='\\t')
chr4 = counts[counts['Chromosome'] == 'chr4']

qA = chr4[chr4['Permissive'].str.startswith('Permissive')]
qB = chr4[chr4['Permissive'].str.startswith('Non-permissive')]

tiers = []
if len(qA) > 0 and len(qB) > 0:
    # Two haplotypes: shorter (contracted) allele first
    shorter = qA if qA['RU_Count'].median() <= qB['RU_Count'].median() else qB
    tiers.append(shorter)
if len(qA) > 0:
    complete_qA = qA[qA['Completeness'] == 'complete']
    if len(complete_qA) > 0:
        tiers.append(complete_qA)
    tiers.append(qA)
# Final fallback: all chr4 reads
tiers.append(chr4)

for i, tier in enumerate(tiers):
    tier['Read_ID'].to_csv(f'tier_{i}.txt', index=False, header=False)

with open('n_tiers', 'w') as f:
    f.write(str(len(tiers)))
PYEOF

    N_TIERS=\$(cat n_tiers)

    # Try each tier until modkit produces output
    for i in \$(seq 0 \$((N_TIERS - 1))); do
        samtools view -b -N "tier_\${i}.txt" -@ ${task.cpus} ${bam} -o selected.bam
        samtools index -@ ${task.cpus} selected.bam

        modkit pileup selected.bam ${sample_id}.methylation.bed \
            --ref ${ref_t2t} \
            --cpg \
            --combine-strands \
            --combine-mods \
            --modified-bases 5mC \
            --include-bed meth_targets.bed \
            --threads ${task.cpus}

        if [ -s "${sample_id}.methylation.bed" ]; then
            break
        fi
    done

    # Final fallback: all reads in the D4Z4 BAM
    if [ ! -s "${sample_id}.methylation.bed" ]; then
        modkit pileup ${bam} ${sample_id}.methylation.bed \
            --ref ${ref_t2t} \
            --cpg \
            --combine-strands \
            --combine-mods \
            --modified-bases 5mC \
            --include-bed meth_targets.bed \
            --threads ${task.cpus}
    fi
    """
}


process REPORT {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    path counts_tsv
    path methylation_bed
    val  sample_id

    output:
    path "${sample_id}.repeats.tsv"
    path "${sample_id}.methylation.tsv"

    script:
    """
    python3 - "${sample_id}" << 'PYEOF'
import sys
import pandas as pd

sample_id = sys.argv[1]

# --- Table 1: D4Z4 Repeat Analysis ---
counts = pd.read_csv(f'{sample_id}.d4z4_counts.tsv', sep='\\t')

# Derive haplotype: only reads with distal anchor (pLAM) can be typed
def get_haplotype(row):
    if row['Completeness'] in ('partial-proximal', 'incomplete'):
        return '-'
    if row['Permissive'].startswith('Permissive'):
        return '4qA'
    elif row['Permissive'].startswith('Non-permissive'):
        return '4qB'
    return '-'

counts['Haplotype'] = counts.apply(get_haplotype, axis=1)
counts['Complete'] = (counts['Completeness'] == 'complete')

rows = []
for (chrom, hap, complete), grp in counts.groupby(['Chromosome', 'Haplotype', 'Complete']):
    ru_min = int(grp['RU_Count'].min())
    ru_max = int(grp['RU_Count'].max())
    ru_counts = f'{ru_min}~{ru_max}' if ru_min != ru_max else str(ru_min)
    rows.append({
        'Chromosome': chrom,
        'Haplotype': hap,
        'RU_Counts': ru_counts,
        'Read_Count': len(grp),
        'Complete': str(complete).upper(),
    })

pd.DataFrame(rows).to_csv(f'{sample_id}.repeats.tsv', sep='\\t', index=False)

# --- Table 2: Methylation Analysis ---
# Methylation sub-regions (T2T-CHM13v2.0 chr4 coordinates)
REGIONS = {
    'DR1':      (193541000, 193541600),   # 29 CpGs, D4Z4 promoter
    'DUX4-PAS': (193543430, 193543634),   # 10 CpGs, polyadenylation signal
}
# Control IQR from healthy reference cohort
CONTROL_IQR = {
    'DR1':      [49.91, 52.29, 54.65],
    'DUX4-PAS': [53.13, 55.33, 59.20],
}

try:
    meth = pd.read_csv(f'{sample_id}.methylation.bed', sep='\\t', header=None, comment='#')
    # bedMethyl: col 0 = chrom, col 1 = start, col 10 = pct_modified
    meth = meth.rename(columns={0: 'chrom', 1: 'start', 10: 'pct_modified'})
except (pd.errors.EmptyDataError, pd.errors.ParserError):
    meth = pd.DataFrame(columns=['chrom', 'start', 'pct_modified'])

meth_rows = []
for region, (rstart, rend) in REGIONS.items():
    sub = meth[(meth['chrom'] == 'chr4') & (meth['start'] >= rstart) & (meth['start'] < rend)]
    mean_pct = round(float(sub['pct_modified'].mean()), 2) if len(sub) > 0 else 'NA'
    q25, q50, q75 = CONTROL_IQR[region]
    meth_rows.append({
        'Region': region,
        'Methylation_Pct': mean_pct,
        'Control_Methylation_IQR [25%, 50%, 75%]': f'[{q25}, {q50}, {q75}]',
    })

pd.DataFrame(meth_rows).to_csv(f'{sample_id}.methylation.tsv', sep='\\t', index=False)
PYEOF
    """
}


workflow {
    if (!params.input) {
        error "Please provide an input BAM/CRAM file with --input"
    }

    def input_file = file(params.input)
    def sample_id  = input_file.baseName

    // Auto-detect index: .cram.crai or .bam.bai
    def idx_ext    = input_file.name.endsWith('.cram') ? '.crai' : '.bai'
    def input_idx  = file("${params.input}${idx_ext}")

    EXTRACT_FSHD_READS(input_file, input_idx, file(params.ref_hg38))
    REALIGN_T2T(EXTRACT_FSHD_READS.out, file(params.ref_t2t), sample_id)
    EXTRACT_D4Z4(REALIGN_T2T.out.bam, REALIGN_T2T.out.bai, sample_id)

    BLAST_ANALYSIS(
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai, sample_id,
        file(params.query_d4z4), file(params.query_plam_polya),
        file(params.query_p13e11)
    )

    METHYLATION(
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai,
        BLAST_ANALYSIS.out, file(params.ref_t2t), sample_id
    )

    REPORT(BLAST_ANALYSIS.out, METHYLATION.out, sample_id)
}
