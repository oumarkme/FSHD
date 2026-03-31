#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// ============================================================================
// Step 1: Extract FSHD-region reads from hg38 BAM/CRAM
// ============================================================================
process EXTRACT_FSHD_READS {

    input:
    path input_file
    path input_index
    path ref_hg38
    val  sample_id

    output:
    path "${sample_id}.subset.bam", emit: bam

    script:
    """
    cat > regions.bed << 'EOF'
chr4	189000000	190500000	4q35_D4Z4_flanking
chr10	133000000	134000000	10q26_D4Z4_flanking
EOF

    REF_OPT=""
    if [[ "${input_file}" == *.cram ]]; then
        REF_OPT="--reference ${ref_hg38}"
    fi

    samtools view -b -h -M -L regions.bed \$REF_OPT -@ ${task.cpus} ${input_file} -o target.bam 2>&1 | \
        grep -v 'cram_decode_seq\\|cram_next_slice' >&2 || true
    samtools view -b -f 4 \$REF_OPT -@ ${task.cpus} ${input_file} -o unmapped.bam 2>&1 | \
        grep -v 'cram_decode_seq\\|cram_next_slice' >&2 || true
    samtools merge -@ ${task.cpus} -o ${sample_id}.subset.bam target.bam unmapped.bam
    """
}


// ============================================================================
// Step 2: Realign to T2T-CHM13v2.0, preserving MM/ML and HP/PS tags
// ============================================================================
process ALIGN_T2T {

    input:
    path subset_bam
    path ref_t2t
    val  sample_id

    output:
    path "${sample_id}.t2t.bam",     emit: bam
    path "${sample_id}.t2t.bam.bai", emit: bai

    script:
    def tags = 'MM,ML,HP,PS'
    """
    samtools fastq -T ${tags} -@ ${task.cpus} ${subset_bam} \
        | minimap2 -ax lr:hq -y --eqx --secondary=no -t ${task.cpus} \
            ${ref_t2t} - \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.t2t.bam

    samtools index -@ ${task.cpus} ${sample_id}.t2t.bam
    """
}


// ============================================================================
// Step 3: Extract D4Z4 candidate reads from T2T BAM
// ============================================================================
process EXTRACT_D4Z4 {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: '*.d4z4.bam*'

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

    samtools view -b -h -M -L d4z4_regions.bed -@ ${task.cpus} ${bam} -o target.bam
    samtools view -b -f 4 -@ ${task.cpus} ${bam} -o unmapped.bam
    samtools merge -@ ${task.cpus} -o ${sample_id}.d4z4.bam target.bam unmapped.bam
    samtools index -@ ${task.cpus} ${sample_id}.d4z4.bam
    """
}


// ============================================================================
// Step 4: Chromosome assignment and feature annotation (BLAST)
//
// Chromosome assignment (three-tier, strongest last):
//   Tier 1 — BAM alignment to T2T (initial, weakest — 98.5% identity)
//   Tier 2 — Restriction motif voting: BlnI + XapI (overrides Tier 1)
//   Tier 3 — Haplotype marker concordance (strongest when both agree):
//            concordant proximal+distal → override all
//            single marker → tiebreaker when motif evidence is LOW
//
// Haplotype determination (DUCKS4-style marker comparison):
//   6 haplotype-specific markers: 3 proximal (D4F104S1) + 3 distal (pLAM)
//   Proximal: 4qA_D4F104S1, 4qB_D4F104S1, 10qA_D4F104S1
//   Distal:   4qA_pLAM (ATTAAA PAS), 4qB_pLAM (no PAS), 10qA_pLAM (ATCAAA PAS)
//   Best bitscore wins; chimeric if proximal ≠ distal type.
//
// Haplotype classification (7 values):
//   4qA-S  — 4qA_pLAM hit, terminal D4Z4 unit ≤2800 bp
//   4qA-L  — 4qA_pLAM hit, terminal D4Z4 unit >2800 bp
//   4qB    — 4qB_pLAM hit (or 4qB_D4F104S1 only)
//   10qA   — 10qA_pLAM hit
//   10qB   — no distal hit, 10qA_D4F104S1 proximal, chr10
//   4q     — no marker hits, chr4
//   10q    — no marker hits, chr10
// ============================================================================
process BLAST_ANNOTATE {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy',
        pattern: '{*.validation.txt,*.complete_reads.bam,*.complete_reads.bam.bai}'

    input:
    path bam
    path bai
    val  sample_id
    path query_d4z4
    path query_haplotype_markers
    path query_p13e11

    output:
    path "${sample_id}.reads_annotated.tsv",      emit: tsv
    path "${sample_id}.validation.txt",            emit: validation
    path "${sample_id}.complete_reads.bam",        emit: complete_bam
    path "${sample_id}.complete_reads.bam.bai",    emit: complete_bai

    script:
    """
    samtools fasta -@ ${task.cpus} ${bam} > reads.fasta
    makeblastdb -in reads.fasta -dbtype nucl -out blastdb
    blastn -db blastdb -query ${query_d4z4}               -out d4z4.blast      -num_threads ${task.cpus} -outfmt 6
    blastn -db blastdb -query ${query_haplotype_markers}   -out haplotype.blast -num_threads ${task.cpus} -outfmt 6
    blastn -db blastdb -query ${query_p13e11}              -out p13e11.blast    -num_threads ${task.cpus} -outfmt 6

    python3 - "${sample_id}" << 'PYEOF'
import re
import sys
import pysam
import pandas as pd

sample_id = sys.argv[1]

blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# --- Tier 1: Chromosome assignment from BAM alignments ---
read_chr = {}
read_len = {}
with pysam.AlignmentFile('${bam}', 'rb') as af:
    for read in af:
        if read.query_name:
            if not read.is_unmapped:
                read_chr[read.query_name] = read.reference_name
            if read.query_length:
                read_len[read.query_name] = read.query_length

# --- D4Z4 RU counting (filter alignment length > 1500 bp) ---
d4z4 = pd.read_csv('d4z4.blast', sep='\\t', header=None, names=blast_cols)
d4z4 = d4z4[d4z4['length'] > 1500]
ru_counts = d4z4.groupby('sseqid').size().reset_index(name='RU_Count')

# --- Haplotype marker BLAST hits ---
PROXIMAL_MARKERS = {'4qA_D4F104S1', '4qB_D4F104S1', '10qA_D4F104S1'}
DISTAL_MARKERS   = {'4qA_pLAM', '4qB_pLAM', '10qA_pLAM'}

hap_markers = pd.DataFrame(columns=blast_cols)
try:
    hap_markers = pd.read_csv('haplotype.blast', sep='\\t', header=None, names=blast_cols)
    hap_markers = hap_markers[hap_markers['pident'] >= 85]
except pd.errors.EmptyDataError:
    pass

# --- P13E11 proximal anchor (additional evidence) ---
proximal_reads = set()
try:
    p13e11_df = pd.read_csv('p13e11.blast', sep='\\t', header=None, names=blast_cols)
    p13e11_df = p13e11_df[(p13e11_df['pident'] >= 85) & (p13e11_df['length'] >= 200)]
    proximal_reads = set(p13e11_df['sseqid'].unique())
except pd.errors.EmptyDataError:
    pass

# --- Per-read best marker lookup ---
def get_best_marker(read_id, marker_set):
    hits = hap_markers[(hap_markers['sseqid'] == read_id) &
                       (hap_markers['qseqid'].isin(marker_set))]
    if len(hits) == 0:
        return None
    return hits.sort_values('bitscore', ascending=False).iloc[0]['qseqid']

# --- Populate Has_Head / Has_Tail from marker hits ---
d4z4_read_ids = set(ru_counts['sseqid'])

read_proximal_marker = {}
read_distal_marker   = {}
distal_reads = set()

for read_id in d4z4_read_ids:
    best_prox = get_best_marker(read_id, PROXIMAL_MARKERS)
    best_dist = get_best_marker(read_id, DISTAL_MARKERS)
    if best_prox:
        read_proximal_marker[read_id] = best_prox
        proximal_reads.add(read_id)
    if best_dist:
        read_distal_marker[read_id] = best_dist
        distal_reads.add(read_id)

# --- Tier 2: Restriction motif voting ---
CHR4_MOTIF  = 'CTAGGTCTAGG'
CHR10_MOTIF = 'CTAGGCCTAGG'
CHR4_RC     = CHR4_MOTIF.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
CHR10_RC    = CHR10_MOTIF.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
XAPI_RE     = re.compile(r'[AG]AATT[CT]')

blni_chr4_counts  = {}
blni_chr10_counts = {}
xapi_counts       = {}

with pysam.AlignmentFile('${bam}', 'rb') as af:
    for read in af:
        if read.query_name not in d4z4_read_ids:
            continue
        if not read.query_sequence:
            continue
        seq = read.query_sequence.upper()
        n_chr4_blni  = seq.count(CHR4_MOTIF)  + seq.count(CHR4_RC)
        n_chr10_blni = seq.count(CHR10_MOTIF) + seq.count(CHR10_RC)
        n_xapi       = len(XAPI_RE.findall(seq))
        blni_chr4_counts[read.query_name]  = n_chr4_blni
        blni_chr10_counts[read.query_name] = n_chr10_blni
        xapi_counts[read.query_name]       = n_xapi
        chr4_score  = n_chr4_blni + n_xapi
        chr10_score = n_chr10_blni
        if chr10_score > chr4_score:
            read_chr[read.query_name] = 'chr10'
        elif chr4_score > chr10_score:
            read_chr[read.query_name] = 'chr4'

# --- Tier 3: Chromosome refinement from haplotype markers ---
# Haplotype markers imply a chromosome of origin. When both proximal and
# distal markers agree, this is the strongest chromosome evidence — it
# overrides BAM alignment and motif voting. When only one marker is
# present, it is used as a tiebreaker when motif evidence is LOW.
MARKER_CHR = {
    '4qA_D4F104S1': 'chr4',  '4qB_D4F104S1': 'chr4',
    '10qA_D4F104S1': 'chr10',
    '4qA_pLAM': 'chr4',      '4qB_pLAM': 'chr4',
    '10qA_pLAM': 'chr10',
}

def _motif_total(read_id):
    return (blni_chr4_counts.get(read_id, 0) +
            blni_chr10_counts.get(read_id, 0) +
            xapi_counts.get(read_id, 0))

for read_id in d4z4_read_ids:
    prox = read_proximal_marker.get(read_id)
    dist = read_distal_marker.get(read_id)
    prox_chr = MARKER_CHR.get(prox)
    dist_chr = MARKER_CHR.get(dist)

    if prox_chr and dist_chr and prox_chr == dist_chr:
        # Concordant markers — strongest evidence, always override
        read_chr[read_id] = prox_chr
    elif prox_chr and dist_chr:
        # Markers disagree (chimeric/translocation) — keep Tier 1/2 assignment
        pass
    elif _motif_total(read_id) == 0:
        # Motif evidence absent — single marker as tiebreaker
        read_chr[read_id] = prox_chr or dist_chr or read_chr.get(read_id, 'Unknown')

# --- Chr confidence ---
def motif_confidence(read_id):
    n4    = blni_chr4_counts.get(read_id, 0)
    n10   = blni_chr10_counts.get(read_id, 0)
    n_xap = xapi_counts.get(read_id, 0)
    total = n4 + n10 + n_xap
    if total == 0:
        return 'LOW'
    chr4_score  = n4 + n_xap
    chr10_score = n10
    if min(chr4_score, chr10_score) == 0:
        return 'HIGH'
    ratio = max(chr4_score, chr10_score) / total
    if ratio >= 0.8:
        return 'HIGH'
    elif ratio >= 0.6:
        return 'MEDIUM'
    return 'LOW'

# --- 4qA-L/S subtyping ---
# pLAM position on read (from 4qA_pLAM BLAST hit)
plam_sstart = {}
if len(hap_markers) > 0:
    plam_hits = hap_markers[hap_markers['qseqid'] == '4qA_pLAM']
    for seq_id, grp in plam_hits.groupby('sseqid'):
        best = grp.sort_values('evalue').iloc[0]
        plam_sstart[seq_id] = min(int(best['sstart']), int(best['send']))

d4z4_hits = d4z4[['sseqid', 'sstart', 'send', 'length']].copy()
d4z4_hits['smin'] = d4z4_hits[['sstart', 'send']].min(axis=1)

def determine_subtype(read_id):
    if read_id not in plam_sstart:
        return '4qA'
    pLAM_pos = plam_sstart[read_id]
    read_d4z4 = d4z4_hits[d4z4_hits['sseqid'] == read_id]
    if len(read_d4z4) == 0:
        return '4qA'
    read_d4z4 = read_d4z4.copy()
    read_d4z4['dist_to_plam'] = (read_d4z4['smin'] - pLAM_pos).abs()
    closest = read_d4z4.sort_values('dist_to_plam').iloc[0]
    return '4qA-L' if int(closest['length']) > 2800 else '4qA-S'

# --- Permissive haplotype check ---
# Permissive = functional PAS (4qA_pLAM with ATTAAA) + non-chimeric.
# Permissiveness is sequence-dependent, not chromosome-dependent:
#   - chr4  + 4qA_pLAM → permissive (standard FSHD1)
#   - chr10 + 4qA_pLAM → permissive (rare translocation; Lemmers et al., 2021)
#   - any   + 10qA_pLAM → non-permissive (ATCAAA, dysfunctional PAS)
#   - any   + 4qB_pLAM  → non-permissive (no PAS)
# 4A166 has ATTAAA PAS but is non-permissive; it appears as chimeric
# (proximal=10qA, distal=4qA) so chimeric reads are always non-permissive.
def is_permissive(haplotype, chimeric, distal_marker):
    if chimeric:
        return False
    return distal_marker == '4qA_pLAM'

# --- Haplotype classification ---
# The haplotype prefix always matches the assigned chromosome:
#   chr4  → 4qA / 4qA-S / 4qA-L / 4qB / 4q
#   chr10 → 10qA / 10qB / 10q
# The Distal_Marker column records which pLAM variant actually matched
# (4qA_pLAM, 4qB_pLAM, 10qA_pLAM), revealing translocation cases.
def classify_haplotype(read_id, chromosome):
    best_dist = read_distal_marker.get(read_id)
    best_prox = read_proximal_marker.get(read_id)
    prefix = '4q' if chromosome == 'chr4' else ('10q' if chromosome == 'chr10' else '')

    # Chimeric: proximal and distal marker types disagree
    chimeric = False
    if best_prox and best_dist:
        prox_type = best_prox.split('_')[0]
        dist_type = best_dist.split('_')[0]
        if prox_type != dist_type:
            chimeric = True

    # Haplotype from distal marker, prefix forced to match chromosome
    if best_dist in ('4qA_pLAM', '10qA_pLAM'):
        # A-type pLAM (either variant)
        if prefix == '4q':
            hap = determine_subtype(read_id) if best_dist == '4qA_pLAM' else '4qA'
        elif prefix == '10q':
            hap = '10qA'
        else:
            hap = 'Unknown'
    elif best_dist == '4qB_pLAM':
        hap = prefix + 'B' if prefix else '4qB'
    elif best_prox:
        # No distal hit — use proximal marker, but force prefix to chromosome
        prox_suffix = 'A' if 'A_D4F104S1' in best_prox else 'B'
        hap = prefix + prox_suffix if prefix else 'Unknown'
    elif read_id in proximal_reads:
        hap = prefix if prefix else 'Unknown'
    else:
        hap = prefix if prefix else 'Unknown'

    return hap, chimeric

# --- Build output ---
rows = []
for _, row in ru_counts.iterrows():
    read_id    = row['sseqid']
    chromosome = read_chr.get(read_id, 'Unknown')
    haplotype, chimeric = classify_haplotype(read_id, chromosome)

    rows.append({
        'Read_ID':          read_id,
        'Chromosome':       chromosome,
        'Chr_Confidence':   motif_confidence(read_id),
        'RU_Count':         int(row['RU_Count']),
        'Haplotype':        haplotype,
        'Permissive':       is_permissive(haplotype, chimeric, read_distal_marker.get(read_id, '')),
        'Read_Length':      int(read_len.get(read_id, 0)),
        'Has_Head':         read_id in proximal_reads,
        'Has_Tail':         read_id in distal_reads,
        'Proximal_Marker':  read_proximal_marker.get(read_id, ''),
        'Distal_Marker':    read_distal_marker.get(read_id, ''),
        'Chimeric':         chimeric,
        'BlnI_Chr4':        blni_chr4_counts.get(read_id, 0),
        'BlnI_Chr10':       blni_chr10_counts.get(read_id, 0),
    })

result = pd.DataFrame(rows, columns=[
    'Read_ID', 'Chromosome', 'Chr_Confidence', 'RU_Count', 'Haplotype',
    'Permissive', 'Read_Length', 'Has_Head', 'Has_Tail', 'Proximal_Marker',
    'Distal_Marker', 'Chimeric', 'BlnI_Chr4', 'BlnI_Chr10',
])

# --- Biological validation ---
A_TYPES = {'4qA-L', '4qA-S', '4qA', '10qA'}
B_TYPES = {'4qB', '10qB'}
GENERIC = {'4q', '10q', 'Unknown'}

result['Flag'] = ''
warnings = []

for chrom in ['chr4', 'chr10']:
    chrom_df = result[result['Chromosome'] == chrom]
    if len(chrom_df) == 0:
        continue

    # Check 1: At most 2 typed haplotypes per chromosome
    typed = chrom_df[~chrom_df['Haplotype'].isin(GENERIC)]
    distinct = sorted(typed['Haplotype'].unique())
    if len(distinct) > 2:
        warnings.append(f'WARN: {chrom} has {len(distinct)} typed haplotypes: {distinct}')
    elif len(distinct) <= 2:
        warnings.append(f'PASS: {chrom} has {len(distinct)} typed haplotype(s): {distinct}')

    for hap in chrom_df['Haplotype'].unique():
        grp = chrom_df[chrom_df['Haplotype'] == hap]

        # Check 2: Completeness consistency for A-type alleles
        if hap in A_TYPES:
            complete = grp[grp['Has_Head'] & grp['Has_Tail']]
            if len(complete) > 0:
                min_complete_len = int(complete['Read_Length'].min())
                max_complete_ru  = int(complete['RU_Count'].max())
                for idx in grp.index:
                    r = grp.loc[idx]
                    if r['Has_Tail'] and not r['Has_Head']:
                        flags = []
                        if r['Read_Length'] >= min_complete_len:
                            flags.append('longer_than_complete')
                        if r['RU_Count'] > max_complete_ru:
                            flags.append('RU_exceeds_complete')
                        if flags:
                            result.at[idx, 'Flag'] = '; '.join(flags)

        # Check 3: Completeness consistency for B-type alleles
        if hap in B_TYPES:
            complete = grp[grp['Has_Head'] & grp['Has_Tail']]
            if len(complete) > 0:
                min_complete_len = int(complete['Read_Length'].min())
                max_complete_ru  = int(complete['RU_Count'].max())
                for idx in grp.index:
                    r = grp.loc[idx]
                    if not (r['Has_Head'] and r['Has_Tail']):
                        flags = []
                        if r['Read_Length'] >= min_complete_len:
                            flags.append('longer_than_complete')
                        if r['RU_Count'] > max_complete_ru:
                            flags.append('RU_exceeds_complete')
                        if flags:
                            existing = result.at[idx, 'Flag']
                            result.at[idx, 'Flag'] = (existing + '; ' if existing else '') + '; '.join(flags)

    # Check 4: Fragment RU sanity
    all_complete = chrom_df[chrom_df['Has_Head'] & chrom_df['Has_Tail'] &
                            (~chrom_df['Haplotype'].isin(GENERIC))]
    if len(all_complete) > 0:
        max_allele_ru = int(all_complete['RU_Count'].max())
        for idx in chrom_df.index:
            r = chrom_df.loc[idx]
            if r['Haplotype'] in GENERIC and r['RU_Count'] > max_allele_ru:
                existing = result.at[idx, 'Flag']
                flag = 'RU_exceeds_all_alleles'
                result.at[idx, 'Flag'] = (existing + '; ' + flag) if existing else flag

# Check 5: Chimeric reads
for idx in result.index:
    if result.loc[idx, 'Chimeric']:
        existing = result.at[idx, 'Flag']
        flag = 'chimeric'
        result.at[idx, 'Flag'] = (existing + '; ' + flag) if existing else flag

result = result.sort_values('RU_Count', ascending=False)
result.to_csv(f'{sample_id}.reads_annotated.tsv', sep='\\t', index=False)

# --- Validation report ---
with open(f'{sample_id}.validation.txt', 'w') as vf:
    vf.write(f'=== Biological Validation: {sample_id} ===\\n\\n')

    # --- Section 1: Summary ---
    total_reads = len(result)
    total_complete = len(result[result['Has_Head'] & result['Has_Tail']])
    total_head_only = len(result[result['Has_Head'] & ~result['Has_Tail']])
    total_tail_only = len(result[~result['Has_Head'] & result['Has_Tail']])
    total_fragment = len(result[~result['Has_Head'] & ~result['Has_Tail']])
    vf.write(f'Summary:\\n')
    vf.write(f'  Total D4Z4 reads: {total_reads}\\n')
    vf.write(f'  Complete (head+tail): {total_complete}\\n')
    vf.write(f'  Head only: {total_head_only}\\n')
    vf.write(f'  Tail only: {total_tail_only}\\n')
    vf.write(f'  Fragment (no anchors): {total_fragment}\\n')
    vf.write('\\n')

    # --- Section 2: Haplotype counts per chromosome ---
    vf.write('Haplotype counts per chromosome:\\n')
    for w in warnings:
        vf.write(f'  [{w.split(":")[0]}] {w.split(":", 1)[1].strip()}\\n')
    vf.write('\\n')

    # --- Section 3: Per-allele summary ---
    vf.write('Per-allele summary:\\n')
    for chrom in ['chr4', 'chr10']:
        chrom_df = result[result['Chromosome'] == chrom]
        if len(chrom_df) == 0:
            continue
        for hap in sorted(chrom_df['Haplotype'].unique()):
            if hap in GENERIC:
                continue
            grp = chrom_df[chrom_df['Haplotype'] == hap]
            complete = grp[grp['Has_Head'] & grp['Has_Tail']]
            incomplete = grp[~(grp['Has_Head'] & grp['Has_Tail'])]
            perm = grp['Permissive'].iloc[0] if 'Permissive' in grp.columns else ''
            chimeric_n = int(grp['Chimeric'].sum()) if 'Chimeric' in grp.columns else 0
            max_ru = int(complete['RU_Count'].max()) if len(complete) > 0 else int(grp['RU_Count'].max())
            ru_source = 'complete' if len(complete) > 0 else 'incomplete'
            vf.write(f'\\n  {chrom}/{hap}:\\n')
            vf.write(f'    Permissive: {perm}\\n')
            vf.write(f'    Total reads: {len(grp)} ({len(complete)} complete, {len(incomplete)} incomplete)\\n')
            vf.write(f'    Max RU: {max_ru} (from {ru_source} reads)\\n')
            if chimeric_n > 0:
                vf.write(f'    Chimeric reads: {chimeric_n}\\n')
    vf.write('\\n')

    # --- Section 4: Complete reads detail ---
    vf.write('Complete reads (Has_Head AND Has_Tail):\\n')
    any_complete = False
    for chrom in ['chr4', 'chr10']:
        chrom_df = result[result['Chromosome'] == chrom]
        for hap in sorted(chrom_df['Haplotype'].unique()):
            if hap in GENERIC:
                continue
            grp = chrom_df[chrom_df['Haplotype'] == hap]
            complete = grp[grp['Has_Head'] & grp['Has_Tail']]
            if len(complete) == 0:
                continue
            any_complete = True
            perm = complete['Permissive'].iloc[0] if 'Permissive' in complete.columns else ''
            vf.write(f'\\n  {chrom}/{hap} (Permissive={perm}): {len(complete)} complete read(s)\\n')
            for _, r in complete.sort_values('RU_Count', ascending=False).iterrows():
                chim_str = '  CHIMERIC' if r.get('Chimeric', False) else ''
                vf.write(f'    {r["Read_ID"]}\\n')
                vf.write(f'      RU={r["RU_Count"]}  Read_Length={r["Read_Length"]}  Chr_Confidence={r["Chr_Confidence"]}\\n')
                vf.write(f'      Proximal_Marker={r["Proximal_Marker"] if r["Proximal_Marker"] else "-"}  '
                         f'Distal_Marker={r["Distal_Marker"] if r["Distal_Marker"] else "-"}{chim_str}\\n')
                vf.write(f'      BlnI_Chr4={r["BlnI_Chr4"]}  BlnI_Chr10={r["BlnI_Chr10"]}\\n')
    if not any_complete:
        vf.write('  No complete reads found.\\n')
    vf.write('\\n')

    # --- Section 5: Low confidence reads ---
    low_conf = result[result['Flag'].astype(str).str.len() > 0]
    if len(low_conf) > 0:
        vf.write(f'Low confidence reads ({len(low_conf)}):\\n')
        vf.write(f'  Possible reasons: longer_than_complete, RU_exceeds_complete,\\n')
        vf.write(f'                    RU_exceeds_all_alleles, chimeric\\n\\n')
        for _, r in low_conf.iterrows():
            vf.write(f'  {r["Read_ID"]}\\n')
            vf.write(f'    {r["Chromosome"]}  {r["Haplotype"]}  RU={r["RU_Count"]}  len={r["Read_Length"]}\\n')
            vf.write(f'    reason: {r["Flag"]}\\n')
    else:
        vf.write('No low confidence reads.\\n')

# --- Extract complete reads to BAM ---
complete_ids = set(result[result['Has_Head'] & result['Has_Tail']]['Read_ID'])
with open('complete_read_ids.txt', 'w') as f:
    for rid in complete_ids:
        f.write(rid + '\\n')
PYEOF

    # Extract complete reads into a separate BAM (empty if none)
    if [ -s complete_read_ids.txt ]; then
        samtools view -b -h -N complete_read_ids.txt -@ ${task.cpus} ${bam} -o ${sample_id}.complete_reads.bam
        samtools index -@ ${task.cpus} ${sample_id}.complete_reads.bam
    else
        samtools view -b -h -@ ${task.cpus} ${bam} -o ${sample_id}.complete_reads.bam -e 'False'
        samtools index -@ ${task.cpus} ${sample_id}.complete_reads.bam
    fi
    """
}


// ============================================================================
// Step 5: Per-allele CpG methylation
//   Uses HP tags carried from the input through ALIGN_T2T.
// ============================================================================
process METHYLATION {

    input:
    path t2t_bam
    path t2t_bai
    path ref_t2t
    path ref_t2t_fai
    val  sample_id

    output:
    path "${sample_id}.methylation_HP1.bed",       emit: hp1_bed
    path "${sample_id}.methylation_HP2.bed",       emit: hp2_bed
    path "${sample_id}.methylation_ungrouped.bed", emit: ungrouped_bed

    script:
    """
    cat > meth_targets.bed << 'EOF'
chr4	193540172	193543634
chr10	134738000	134741148
EOF

    modkit pileup ${t2t_bam} meth_out/ \
        --ref ${ref_t2t} \
        --cpg \
        --combine-strands \
        --modified-bases 5mC \
        --phased \
        --include-bed meth_targets.bed \
        --prefix ${sample_id} \
        --threads ${task.cpus}

    if ls meth_out/${sample_id}*hp1* 1>/dev/null 2>&1; then
        cp meth_out/${sample_id}*hp1* ${sample_id}.methylation_HP1.bed
    else
        touch ${sample_id}.methylation_HP1.bed
    fi

    if ls meth_out/${sample_id}*hp2* 1>/dev/null 2>&1; then
        cp meth_out/${sample_id}*hp2* ${sample_id}.methylation_HP2.bed
    else
        touch ${sample_id}.methylation_HP2.bed
    fi

    # Modkit names the non-phased partition *_combined.bedmethyl
    if ls meth_out/${sample_id}*combined* 1>/dev/null 2>&1; then
        cp meth_out/${sample_id}*combined* ${sample_id}.methylation_ungrouped.bed
    elif ls meth_out/${sample_id}*ungrouped* 1>/dev/null 2>&1; then
        cp meth_out/${sample_id}*ungrouped* ${sample_id}.methylation_ungrouped.bed
    else
        touch ${sample_id}.methylation_ungrouped.bed
    fi
    """
}


// ============================================================================
// Step 6: Generate report
//   repeats.tsv    — grouped by (Chromosome, Haplotype, Has_Head, Has_Tail)
//   methylation.tsv — DR1 and DUX4-PAS means with atlas/literature references
// ============================================================================
process REPORT {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    path annotated_tsv
    path hp1_bed
    path hp2_bed
    path ungrouped_bed
    val  sample_id

    output:
    path "${sample_id}.repeats.tsv",     emit: repeats
    path "${sample_id}.methylation.tsv", emit: methylation

    script:
    """
    python3 - "${sample_id}" << 'PYEOF'
import sys
import pandas as pd

sample_id = sys.argv[1]

df = pd.read_csv(f'{sample_id}.reads_annotated.tsv', sep='\\t')

# --- repeats.tsv ---
# Group by (Chromosome, Haplotype, Permissive, Has_Head, Has_Tail).
# Report the maximum RU_Count, read count, and flagged count per group.
group_cols = ['Chromosome', 'Haplotype', 'Has_Head', 'Has_Tail']
repeats = (
    df.groupby(group_cols, dropna=False)
    .agg(
        RU_Count=('RU_Count', 'max'),
        Read_Count=('Read_ID', 'count'),
    )
    .reset_index()
    .rename(columns={'Has_Head': 'Head', 'Has_Tail': 'Tail'})
)
repeats = repeats[['Chromosome', 'Haplotype', 'RU_Count', 'Head', 'Tail', 'Read_Count']]
# Add ≥ prefix for incomplete reads (RU count is a lower bound)
repeats['RU_Count'] = repeats.apply(
    lambda r: f'>={int(r["RU_Count"])}' if not (r['Head'] and r['Tail']) else str(int(r['RU_Count'])),
    axis=1
)
chr_order = {'chr4': 0, 'chr10': 1}
repeats['_chr_sort'] = repeats['Chromosome'].map(chr_order).fillna(99)
repeats = repeats.sort_values(
    ['_chr_sort', 'Haplotype', 'Head', 'Tail'],
    ascending=[True, True, False, False],
).drop(columns=['_chr_sort'])
repeats.to_csv(f'{sample_id}.repeats.tsv', sep='\\t', index=False)

# --- methylation.tsv ---
# DR1 and DUX4-PAS mean methylation per haplotype partition.
# Reference atlas: ONT WGS, ~45-49 healthy controls, Gaussian-smoothed (sigma=150bp)
# Literature: Caputo et al. 2022 (Cells 11:4114) — thresholds from 150 controls
REGIONS = {
    'DR1':      (193541000, 193541600),
    'DUX4-PAS': (193543430, 193543634),
}
ATLAS_REF = {
    'DR1':      (54.70, 1.97),
    'DUX4-PAS': (51.81, 0.23),
}
# IQR computed from all CpG REF_METH values in the reference atlas (ref_atlas_150.tsv.gz)
# hg38 coordinates via UCSC liftover: DR1 chr4:190173187-190173795, DUX4-PAS chr4:190175625-190175829
ATLAS_IQR = {
    'DR1':      (53.07, 54.89, 56.83),
    'DUX4-PAS': (51.63, 51.75, 51.99),
}
LITERATURE_REF = {
    'DR1':      'CpG1 <= 53%',
    'DUX4-PAS': 'CpG6 <= 78%',
}

def read_bedmethyl(path):
    try:
        bm = pd.read_csv(path, sep='\\t', header=None, comment='#')
        return bm.rename(columns={
            0: 'chrom', 1: 'start', 9: 'N_valid', 10: 'pct_modified'
        })
    except (pd.errors.EmptyDataError, pd.errors.ParserError, FileNotFoundError):
        return pd.DataFrame(columns=['chrom', 'start', 'N_valid', 'pct_modified'])

bed_hp1 = read_bedmethyl(f'{sample_id}.methylation_HP1.bed')
bed_hp2 = read_bedmethyl(f'{sample_id}.methylation_HP2.bed')
bed_ung = read_bedmethyl(f'{sample_id}.methylation_ungrouped.bed')

# Combined: weighted average across HP1 + HP2 + ungrouped at each CpG position
all_beds = pd.concat([bed_hp1, bed_hp2, bed_ung], ignore_index=True)
if len(all_beds) > 0 and 'N_valid' in all_beds.columns and all_beds['N_valid'].notna().any():
    all_beds['n_modified'] = all_beds['N_valid'] * all_beds['pct_modified'] / 100.0
    bed_combined = (
        all_beds.groupby(['chrom', 'start'])
        .agg(N_valid=('N_valid', 'sum'), n_modified=('n_modified', 'sum'))
        .reset_index()
    )
    bed_combined['pct_modified'] = (
        bed_combined['n_modified'] / bed_combined['N_valid'] * 100.0
    )
else:
    bed_combined = pd.DataFrame(
        columns=['chrom', 'start', 'N_valid', 'pct_modified']
    )

meth_rows = []
for hp_label, meth in [('HP1', bed_hp1), ('HP2', bed_hp2), ('Combined', bed_combined)]:
    for region, (rstart, rend) in REGIONS.items():
        sub = meth[
            (meth['chrom'] == 'chr4') &
            (meth['start'] >= rstart) &
            (meth['start'] < rend)
        ]
        mean_pct = round(float(sub['pct_modified'].mean()), 2) if len(sub) > 0 else 'NA'
        n_cpgs   = len(sub)
        atlas_mean, atlas_sd = ATLAS_REF[region]
        q1, q2, q3 = ATLAS_IQR[region]
        meth_rows.append({
            'Haplotype':             hp_label,
            'Region':                region,
            'Methylation':           mean_pct,
            'N_CpG':                 n_cpgs,
            'Ref_Atlas (Mean (SD) %)': f'{atlas_mean}% ({atlas_sd}%)',
            'Ref_Atlas (IQR %)':     f'[{q1}%, {q2}%, {q3}%]',
            'Lit_Caputo (%)':        LITERATURE_REF[region],
        })

pd.DataFrame(meth_rows, columns=[
    'Haplotype', 'Region', 'Methylation', 'N_CpG',
    'Ref_Atlas (Mean (SD) %)', 'Ref_Atlas (IQR %)', 'Lit_Caputo (%)',
]).to_csv(f'{sample_id}.methylation.tsv', sep='\\t', index=False)
PYEOF
    """
}


// ============================================================================
// Workflow
//   Input BAM/CRAM must already have HP/PS tags from upstream phasing.
//   Tags are carried through realignment via samtools fastq -T and minimap2 -y.
// ============================================================================
workflow {
    if (!params.input) {
        error "Please provide an input BAM/CRAM file with --input"
    }

    def input_file = file(params.input)
    def sample_id  = input_file.baseName
    def ref_t2t     = file(params.ref_t2t)
    def ref_t2t_fai = file("${params.ref_t2t}.fai")
    def ref_hg38    = file(params.ref_hg38)

    def idx_ext   = input_file.name.endsWith('.cram') ? '.crai' : '.bai'
    def input_idx = file("${params.input}${idx_ext}")

    // Step 1: Extract FSHD-region reads
    EXTRACT_FSHD_READS(input_file, input_idx, ref_hg38, sample_id)

    // Step 2: Realign to T2T (HP/PS tags preserved)
    ALIGN_T2T(EXTRACT_FSHD_READS.out.bam, ref_t2t, sample_id)

    // Step 3: Extract D4Z4 candidates
    EXTRACT_D4Z4(ALIGN_T2T.out.bam, ALIGN_T2T.out.bai, sample_id)

    // Step 4: BLAST annotation with haplotype-specific markers
    BLAST_ANNOTATE(
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai, sample_id,
        file(params.query_d4z4), file(params.query_haplotype_markers),
        file(params.query_p13e11)
    )

    // Step 5: CpG methylation (HP tags from T2T BAM)
    METHYLATION(
        ALIGN_T2T.out.bam, ALIGN_T2T.out.bai,
        ref_t2t, ref_t2t_fai, sample_id
    )

    // Step 6: Report (repeats + methylation summary)
    REPORT(
        BLAST_ANNOTATE.out.tsv,
        METHYLATION.out.hp1_bed,
        METHYLATION.out.hp2_bed,
        METHYLATION.out.ungrouped_bed,
        sample_id
    )
}
