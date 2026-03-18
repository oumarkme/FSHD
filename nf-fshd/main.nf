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
// ============================================================================
process BLAST_ANNOTATE {

    input:
    path bam
    path bai
    val  sample_id
    path query_d4z4
    path query_plam_polya
    path query_p13e11

    output:
    path "${sample_id}.reads_annotated.tsv", emit: tsv

    script:
    """
    samtools fasta -@ ${task.cpus} ${bam} > reads.fasta
    makeblastdb -in reads.fasta -dbtype nucl -out blastdb
    blastn -db blastdb -query ${query_d4z4}       -out d4z4.blast   -num_threads ${task.cpus} -outfmt 6
    blastn -db blastdb -query ${query_plam_polya}  -out polya.blast  -num_threads ${task.cpus} -outfmt 6
    blastn -db blastdb -query ${query_p13e11}      -out p13e11.blast -num_threads ${task.cpus} -outfmt 6

    python3 - "${sample_id}" << 'PYEOF'
import sys
import pysam
import pandas as pd

sample_id = sys.argv[1]

blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# --- Chromosome assignment from BAM alignments ---
read_chr = {}
read_len = {}
with pysam.AlignmentFile('${bam}', 'rb') as af:
    for read in af:
        if read.query_name:
            if not read.is_unmapped:
                read_chr[read.query_name] = read.reference_name
            if read.query_length:
                read_len[read.query_name] = read.query_length

# --- D4Z4 RU counting (filter alignment length > 1500) ---
d4z4 = pd.read_csv('d4z4.blast', sep='\\t', header=None, names=blast_cols)
d4z4 = d4z4[d4z4['length'] > 1500]
ru_counts = d4z4.groupby('sseqid').size().reset_index(name='RU_Count')

# --- Completeness: proximal (P13E11) + distal (pLAM) ---
proximal_reads = set()
distal_reads = set()
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

# --- Permissive haplotype check (ALL reads with pLAM hit) ---
# pLAM-polyA.fa contains two sequences:
#   hapA = 4qA PAS (ATTAAA) -- functional, permissive, chr4-specific
#   hapB = 10qA PAS (ATCAAA) -- non-functional, chr10
# 4qB lacks pLAM entirely and produces no hit to either query.
haplotype_map = {}
for seq_id, grp in polya.groupby('sseqid'):
    hap_a = grp[grp['qseqid'] == 'hapA']
    hap_b = grp[grp['qseqid'] == 'hapB']
    if len(hap_a) >= 1 and len(hap_b) >= 1:
        diff = hap_a.iloc[0]['mismatch'] - hap_b.iloc[0]['mismatch']
        if diff < 0:
            haplotype_map[seq_id] = '4qA'    # ATTAAA PAS -> chr4 permissive
        elif diff > 0:
            haplotype_map[seq_id] = '10qA'   # ATCAAA PAS -> chr10
        else:
            haplotype_map[seq_id] = 'Unknown'

# --- BlnI-based chr4/chr10 correction ---
# Chr4 D4Z4 units lack BlnI site (B-): ...CTAGGTCTAGG...
# Chr10 D4Z4 units have BlnI site (B+): ...CTAGGCCTAGG...
# Count diagnostic motifs per read; majority vote determines chromosome.
def rev_comp(seq):
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]

CHR4_MOTIF  = 'CTAGGTCTAGG'
CHR10_MOTIF = 'CTAGGCCTAGG'
CHR4_RC     = rev_comp(CHR4_MOTIF)
CHR10_RC    = rev_comp(CHR10_MOTIF)

d4z4_read_ids = set(ru_counts['sseqid'])
with pysam.AlignmentFile('${bam}', 'rb') as af:
    for read in af:
        if read.query_name in d4z4_read_ids and read.query_sequence:
            seq = read.query_sequence.upper()
            n4  = seq.count(CHR4_MOTIF)  + seq.count(CHR4_RC)
            n10 = seq.count(CHR10_MOTIF) + seq.count(CHR10_RC)
            if n10 > n4:
                read_chr[read.query_name] = 'chr10'
            elif n4 > n10:
                read_chr[read.query_name] = 'chr4'
            # if equal, keep BAM alignment

# --- Haplotype-based correction ---
# 4qA reads must be chr4; 10qA reads must be chr10
for seq_id, hap in haplotype_map.items():
    if hap == '4qA' and read_chr.get(seq_id) != 'chr4':
        read_chr[seq_id] = 'chr4'
    elif hap == '10qA' and read_chr.get(seq_id) != 'chr10':
        read_chr[seq_id] = 'chr10'

ru_counts['Chromosome'] = ru_counts['sseqid'].map(read_chr).fillna('Unknown')
ru_counts['Completeness'] = ru_counts['sseqid'].apply(classify_completeness)
ru_counts['Haplotype'] = ru_counts['sseqid'].map(haplotype_map).fillna('-')
ru_counts['Read_Length'] = ru_counts['sseqid'].map(read_len).fillna(0).astype(int)
ru_counts['Has_Head'] = ru_counts['sseqid'].isin(proximal_reads)
ru_counts['Has_Tail'] = ru_counts['sseqid'].isin(distal_reads)

result = ru_counts.rename(columns={'sseqid': 'Read_ID'})[
    ['Read_ID', 'Chromosome', 'RU_Count', 'Completeness', 'Haplotype', 'Read_Length', 'Has_Head', 'Has_Tail']
]
result = result.sort_values('RU_Count', ascending=False)
result.to_csv(f'{sample_id}.reads_annotated.tsv', sep='\\t', index=False)
PYEOF
    """
}


// ============================================================================
// Step 5: Three-phase allele assignment
//   Phase 1: Separate alleles (HP → RU-anchor → Jenks → single)
//   Phase 2: Assign unphased reads (RU-constraint → minimap2 → RU-proximity)
//   Phase 3: Per-allele counting strategy (read-level vs assembly)
// ============================================================================
process ALLELE_ASSIGN {

    input:
    path annotated_tsv
    path t2t_bam
    path t2t_bai
    path d4z4_bam
    path d4z4_bai
    val  sample_id

    output:
    path "${sample_id}.allele_assignments.tsv", emit: tsv
    path "${sample_id}.allele_strategy.json",   emit: strategy

    script:
    """
    # ==== Phase 1: Separate alleles by HP / RU-anchor / Jenks ====
    python3 - "${sample_id}" << 'PYEOF'
import sys
import json
import pandas as pd
import pysam
import jenkspy

sample_id = sys.argv[1]

df = pd.read_csv(f'{sample_id}.reads_annotated.tsv', sep='\\t')

hp_tags = {}
with pysam.AlignmentFile(f'{sample_id}.t2t.bam', 'rb') as af:
    for read in af:
        if read.query_name and read.has_tag('HP'):
            hp_tags[read.query_name] = read.get_tag('HP')

assignments = []
phase1_info = {}

for chrom in ['chr4', 'chr10']:
    chrom_reads = df[(df['Chromosome'] == chrom) & (df['RU_Count'] >= 1)].copy()

    if len(chrom_reads) == 0:
        phase1_info[chrom] = {'separation_method': 'empty'}
        continue

    chrom_reads['HP_Tag'] = chrom_reads['Read_ID'].map(hp_tags)
    has_any_hp = chrom_reads['HP_Tag'].notna().any()

    if has_any_hp:
        chrom_reads['Allele'] = chrom_reads['HP_Tag'].apply(
            lambda x: 'HAP1' if x == 1 else ('HAP2' if x == 2 else 'unphased')
        )
        chrom_reads['Assignment_Method'] = chrom_reads['HP_Tag'].apply(
            lambda x: 'hp' if pd.notna(x) else 'unphased'
        )
        separation_method = 'hp'

    elif (chrom_reads['Has_Head'] & chrom_reads['Has_Tail']).any():
        complete_ru = chrom_reads.loc[
            chrom_reads['Has_Head'] & chrom_reads['Has_Tail'], 'RU_Count'].median()
        chrom_reads['Allele'] = chrom_reads['RU_Count'].apply(
            lambda x: 'HAP1' if x <= complete_ru + 2 else 'HAP2'
        )
        chrom_reads['Assignment_Method'] = 'ru-anchor'
        separation_method = 'ru-anchor'

    else:
        ru_vals = chrom_reads['RU_Count'].values
        if len(chrom_reads) >= 2 and ru_vals.std() >= 1.5:
            breaks = jenkspy.jenks_breaks(ru_vals.tolist(), n_classes=2)
            threshold = breaks[1]
            chrom_reads['Allele'] = chrom_reads['RU_Count'].apply(
                lambda x: 'HAP1' if x <= threshold else 'HAP2'
            )
            chrom_reads['Assignment_Method'] = 'jenks'
            separation_method = 'jenks'
        else:
            chrom_reads['Allele'] = 'HAP1'
            chrom_reads['Assignment_Method'] = 'single'
            separation_method = 'single'

    # Write per-chromosome read lists for minimap2 alignment
    for allele_label in ['HAP1', 'HAP2', 'unphased']:
        ids = chrom_reads.loc[chrom_reads['Allele'] == allele_label, 'Read_ID']
        ids.to_csv(f'{chrom}_{allele_label}.txt', index=False, header=False)

    chrom_reads = chrom_reads.drop(columns=['HP_Tag'])
    assignments.append(chrom_reads)
    phase1_info[chrom] = {'separation_method': separation_method}

if assignments:
    result = pd.concat(assignments, ignore_index=True)
else:
    result = pd.DataFrame(columns=['Read_ID', 'Chromosome', 'RU_Count', 'Completeness',
                                    'Haplotype', 'Read_Length', 'Has_Head', 'Has_Tail',
                                    'Allele', 'Assignment_Method'])

result.to_csv('phase1_assignments.tsv', sep='\\t', index=False)

with open('phase1_info.json', 'w') as f:
    json.dump(phase1_info, f)
PYEOF

    # ==== Alignment-based unphased read assignment ====
    for C in chr4 chr10; do
        if [[ ! -s "\${C}_unphased.txt" ]]; then
            continue
        fi

        N_UNPHASED=\$(wc -l < "\${C}_unphased.txt")
        N_HAP1=\$(wc -l < "\${C}_HAP1.txt" 2>/dev/null || echo 0)
        N_HAP2=\$(wc -l < "\${C}_HAP2.txt" 2>/dev/null || echo 0)

        if [[ \$N_UNPHASED -eq 0 ]]; then
            continue
        fi

        # Extract FASTA for each group
        samtools view -b -N "\${C}_unphased.txt" -@ ${task.cpus} ${d4z4_bam} \
            | samtools fasta -@ ${task.cpus} > "\${C}_unphased.fasta"

        if [[ \$N_HAP1 -gt 0 ]]; then
            samtools view -b -N "\${C}_HAP1.txt" -@ ${task.cpus} ${d4z4_bam} \
                | samtools fasta -@ ${task.cpus} > "\${C}_HAP1.fasta"
        else
            touch "\${C}_HAP1.fasta"
        fi

        if [[ \$N_HAP2 -gt 0 ]]; then
            samtools view -b -N "\${C}_HAP2.txt" -@ ${task.cpus} ${d4z4_bam} \
                | samtools fasta -@ ${task.cpus} > "\${C}_HAP2.fasta"
        else
            touch "\${C}_HAP2.fasta"
        fi

        # Align unphased reads against each allele pool
        if [[ -s "\${C}_HAP1.fasta" ]]; then
            minimap2 -x map-ont --secondary=no -t ${task.cpus} \
                "\${C}_HAP1.fasta" "\${C}_unphased.fasta" 2>/dev/null \
                > "\${C}_vs_hap1.paf" || true
        else
            touch "\${C}_vs_hap1.paf"
        fi

        if [[ -s "\${C}_HAP2.fasta" ]]; then
            minimap2 -x map-ont --secondary=no -t ${task.cpus} \
                "\${C}_HAP2.fasta" "\${C}_unphased.fasta" 2>/dev/null \
                > "\${C}_vs_hap2.paf" || true
        else
            touch "\${C}_vs_hap2.paf"
        fi
    done

    # ==== Phase 2 + 3: Reassign unphased reads + determine strategy ====
    python3 - "${sample_id}" << 'PYEOF2'
import sys
import json
import pandas as pd

sample_id = sys.argv[1]

df = pd.read_csv('phase1_assignments.tsv', sep='\\t')
with open('phase1_info.json') as f:
    phase1_info = json.load(f)

def parse_paf_scores(paf_path):
    scores = {}
    try:
        with open(paf_path) as f:
            for line in f:
                if not line.strip():
                    continue
                fields = line.strip().split('\\t')
                read_id = fields[0]
                # matching bases (col 10) as score
                match_bases = int(fields[9])
                if read_id not in scores or match_bases > scores[read_id]:
                    scores[read_id] = match_bases
    except FileNotFoundError:
        pass
    return scores

strategy = {}

for chrom in ['chr4', 'chr10']:
    info = phase1_info.get(chrom, {})
    sep_method = info.get('separation_method', 'empty')

    if sep_method == 'empty':
        strategy[chrom] = {'HAP1': None, 'HAP2': None}
        continue

    chrom_mask = df['Chromosome'] == chrom

    # Phase 2: Reassign unphased reads
    if sep_method == 'hp':
        unphased_mask = chrom_mask & (df['Allele'] == 'unphased')

        if unphased_mask.any():
            hap1_reads = df[chrom_mask & (df['Allele'] == 'HAP1')]
            hap2_reads = df[chrom_mask & (df['Allele'] == 'HAP2')]

            # Biological constraint: if an allele has complete reads,
            # its array is fully spanned — unphased reads with RU far
            # above the complete read's RU cannot be from that allele.
            hap1_complete = hap1_reads[hap1_reads['Has_Head'] & hap1_reads['Has_Tail']]
            hap2_complete = hap2_reads[hap2_reads['Has_Head'] & hap2_reads['Has_Tail']]
            hap1_max_complete_ru = int(hap1_complete['RU_Count'].max()) if len(hap1_complete) > 0 else None
            hap2_max_complete_ru = int(hap2_complete['RU_Count'].max()) if len(hap2_complete) > 0 else None

            hap1_scores = parse_paf_scores(f'{chrom}_vs_hap1.paf')
            hap2_scores = parse_paf_scores(f'{chrom}_vs_hap2.paf')
            hap1_median = hap1_reads['RU_Count'].median() if len(hap1_reads) > 0 else 0
            hap2_median = hap2_reads['RU_Count'].median() if len(hap2_reads) > 0 else 0

            for idx in df.index[unphased_mask]:
                ru = df.at[idx, 'RU_Count']
                read_id = df.at[idx, 'Read_ID']

                # Force-assign if RU exceeds a complete allele's span
                if hap1_max_complete_ru is not None and ru > hap1_max_complete_ru + 2:
                    df.at[idx, 'Allele'] = 'HAP2'
                    df.at[idx, 'Assignment_Method'] = 'ru-constraint'
                    continue
                if hap2_max_complete_ru is not None and ru > hap2_max_complete_ru + 2:
                    df.at[idx, 'Allele'] = 'HAP1'
                    df.at[idx, 'Assignment_Method'] = 'ru-constraint'
                    continue

                # Alignment-based assignment for ambiguous reads
                s1 = hap1_scores.get(read_id, 0)
                s2 = hap2_scores.get(read_id, 0)

                if s1 > s2:
                    target = 'HAP1'
                    method = 'alignment'
                elif s2 > s1:
                    target = 'HAP2'
                    method = 'alignment'
                else:
                    target = 'HAP1' if abs(ru - hap1_median) <= abs(ru - hap2_median) else 'HAP2'
                    method = 'ru-proximity'

                df.at[idx, 'Allele'] = target
                df.at[idx, 'Assignment_Method'] = method

    # Phase 3: Determine strategy per allele
    strat = {}
    for allele in ['HAP1', 'HAP2']:
        allele_reads = df[chrom_mask & (df['Allele'] == allele)]
        if len(allele_reads) == 0:
            strat[allele] = None
        elif (allele_reads['Has_Head'] & allele_reads['Has_Tail']).any():
            strat[allele] = 'read-level'
        else:
            strat[allele] = 'assembly'
    strategy[chrom] = strat

result = df[['Read_ID', 'Chromosome', 'RU_Count', 'Completeness', 'Haplotype',
              'Read_Length', 'Has_Head', 'Has_Tail', 'Allele', 'Assignment_Method']]
result.to_csv(f'{sample_id}.allele_assignments.tsv', sep='\\t', index=False)

with open(f'{sample_id}.allele_strategy.json', 'w') as f:
    json.dump(strategy, f, indent=2)
PYEOF2
    """
}


// ============================================================================
// Step 6: Read-level counting via consensus BLAST
// ============================================================================
process READLEVEL_COUNT {

    input:
    path d4z4_bam
    path d4z4_bai
    path allele_assignments
    path query_d4z4
    path query_plam_polya
    path query_p13e11
    val  sample_id
    val  chromosome
    val  allele

    output:
    path "${sample_id}.${chromosome}.${allele}.readlevel_ru_count.txt", emit: ru_count

    script:
    """
    python3 - "${sample_id}" "${chromosome}" "${allele}" << 'PYEOF'
import sys
import pandas as pd

sample_id, chrom, allele = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.read_csv(f'{sample_id}.allele_assignments.tsv', sep='\\t')
reads = df[(df['Chromosome'] == chrom) & (df['Allele'] == allele)]

reads['Read_ID'].to_csv('allele_reads.txt', index=False, header=False)

spanning = reads[reads['Completeness'].isin(['complete', 'partial-distal'])]
spanning_sorted = spanning.sort_values('Read_Length', ascending=False)
spanning_sorted['Read_ID'].to_csv('spanning_reads.txt', index=False, header=False)

if len(spanning_sorted) > 0:
    with open('ref_read.txt', 'w') as f:
        f.write(spanning_sorted.iloc[0]['Read_ID'])
    with open('n_spanning', 'w') as f:
        f.write(str(len(spanning_sorted)))
else:
    with open('n_spanning', 'w') as f:
        f.write('0')
PYEOF

    N_SPANNING=\$(cat n_spanning)

    samtools view -b -N allele_reads.txt -@ ${task.cpus} ${d4z4_bam} -o allele.bam
    samtools fasta -@ ${task.cpus} allele.bam > allele.fasta

    if [[ \$N_SPANNING -ge 1 ]]; then
        REF_READ=\$(cat ref_read.txt)
        seqkit grep -p "\$REF_READ" allele.fasta > raw_ref.fasta

        if [[ \$N_SPANNING -ge 3 ]]; then
            samtools view -b -N spanning_reads.txt -@ ${task.cpus} ${d4z4_bam} -o spanning.bam
            samtools fasta -@ ${task.cpus} spanning.bam > spanning.fasta
            seqkit grep -v -p "\$REF_READ" spanning.fasta > rest.fasta

            minimap2 -ax map-ont -t ${task.cpus} raw_ref.fasta rest.fasta > aln.sam
            racon --no-trimming -t ${task.cpus} rest.fasta aln.sam raw_ref.fasta \
                > consensus.fasta
        else
            cp raw_ref.fasta consensus.fasta
        fi
    else
        echo ">${chromosome}_${allele}" > consensus.fasta
        echo "" >> consensus.fasta
    fi

    HAS_SEQ=false
    if grep -v '^>' consensus.fasta | grep -q '[ACGTacgt]'; then
        HAS_SEQ=true
    fi

    if [[ "\$HAS_SEQ" == "true" ]]; then
        makeblastdb -in consensus.fasta -dbtype nucl -out consensus_db
        blastn -db consensus_db -query ${query_d4z4} -outfmt 6 > consensus_d4z4.blast
        blastn -db consensus_db -query ${query_plam_polya} -outfmt 6 > anchor_polya.blast
        blastn -db consensus_db -query ${query_p13e11} -outfmt 6 > anchor_p13e11.blast
    else
        touch consensus_d4z4.blast anchor_polya.blast anchor_p13e11.blast
    fi

    python3 - "${sample_id}" "${chromosome}" "${allele}" << 'PYEOF2'
import sys
import pandas as pd

sample_id, chrom, allele = sys.argv[1], sys.argv[2], sys.argv[3]
blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
try:
    d4z4 = pd.read_csv('consensus_d4z4.blast', sep='\\t', header=None, names=blast_cols)
    ru_count = len(d4z4[d4z4['length'] > 1500])
except pd.errors.EmptyDataError:
    ru_count = 0

has_proximal = False
has_distal = False
try:
    p = pd.read_csv('anchor_p13e11.blast', sep='\\t', header=None, names=blast_cols)
    has_proximal = len(p) > 0
except pd.errors.EmptyDataError:
    pass
try:
    d = pd.read_csv('anchor_polya.blast', sep='\\t', header=None, names=blast_cols)
    has_distal = len(d) > 0
except pd.errors.EmptyDataError:
    pass

confidence = 'High' if (has_proximal and has_distal) else 'Medium'

with open(f'{sample_id}.{chrom}.{allele}.readlevel_ru_count.txt', 'w') as f:
    f.write(f'{ru_count}\\t{confidence}\\tconsensus\\n')
PYEOF2
    """
}


// ============================================================================
// Step 7: Assembly-based counting (Flye vs longest raw read)
// ============================================================================
process ASSEMBLY_COUNT {

    input:
    path d4z4_bam
    path d4z4_bai
    path allele_assignments
    path query_d4z4
    path query_plam_polya
    path query_p13e11
    val  sample_id
    val  chromosome
    val  allele

    output:
    path "${sample_id}.${chromosome}.${allele}.ru_count.txt", emit: ru_count

    script:
    """
    python3 - "${sample_id}" "${chromosome}" "${allele}" << 'PYEOF'
import sys
import pandas as pd

sample_id, chrom, allele = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.read_csv(f'{sample_id}.allele_assignments.tsv', sep='\\t')
reads = df[(df['Chromosome'] == chrom) & (df['Allele'] == allele)]
reads['Read_ID'].to_csv('allele_reads.txt', index=False, header=False)
max_ru = int(reads['RU_Count'].max()) if len(reads) > 0 else 0
with open('n_reads', 'w') as f:
    f.write(str(len(reads)))
with open('max_raw_ru', 'w') as f:
    f.write(str(max_ru))
PYEOF

    N_READS=\$(cat n_reads)

    samtools view -b -N allele_reads.txt -@ ${task.cpus} ${d4z4_bam} -o allele.bam
    samtools fastq -@ ${task.cpus} allele.bam > allele.fastq

    ASSEMBLY_OK=false

    if [[ \$N_READS -ge 3 ]]; then
        flye --nano-hq allele.fastq \
            --genome-size 300k \
            --min-overlap 5000 \
            --threads ${task.cpus} \
            --out-dir flye_out \
        && ASSEMBLY_OK=true || true
    fi

    MAX_RAW_RU=\$(cat max_raw_ru)
    USE_FLYE=false

    if [[ "\$ASSEMBLY_OK" == "true" && -s flye_out/assembly.fasta ]]; then
        seqkit sort -l -r flye_out/assembly.fasta | seqkit head -n 1 > contig.fasta

        makeblastdb -in contig.fasta -dbtype nucl -out contig_db
        blastn -db contig_db -query ${query_d4z4} -outfmt 6 > contig_d4z4.blast
        blastn -db contig_db -query ${query_p13e11} -outfmt 6 > contig_p13e11.blast
        blastn -db contig_db -query ${query_plam_polya} -outfmt 6 > contig_polya.blast

        python3 - "${sample_id}" "${chromosome}" "${allele}" << 'PYEOF2'
import sys
import pandas as pd

sample_id, chrom, allele = sys.argv[1], sys.argv[2], sys.argv[3]
blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

try:
    d4z4 = pd.read_csv('contig_d4z4.blast', sep='\\t', header=None, names=blast_cols)
    ru_blast = len(d4z4[d4z4['length'] > 1500])
except pd.errors.EmptyDataError:
    ru_blast = 0

ru_length = None
try:
    p13 = pd.read_csv('contig_p13e11.blast', sep='\\t', header=None, names=blast_cols)
    plam = pd.read_csv('contig_polya.blast', sep='\\t', header=None, names=blast_cols)
    if len(p13) > 0 and len(plam) > 0:
        prox_start = min(p13['sstart'].min(), p13['send'].min())
        dist_end = max(plam['sstart'].max(), plam['send'].max())
        array_length = abs(dist_end - prox_start)
        ru_length = round(array_length / 3300)
except (pd.errors.EmptyDataError, ValueError):
    pass

if ru_blast > 0 and ru_length is not None and abs(ru_blast - ru_length) <= 2:
    confidence = 'Medium'
elif ru_blast > 0:
    confidence = 'Low'
else:
    confidence = 'Low'

with open('flye_ru', 'w') as f:
    f.write(str(ru_blast))
with open('flye_confidence', 'w') as f:
    f.write(confidence)
PYEOF2

        FLYE_RU=\$(cat flye_ru)
        if [[ \$FLYE_RU -ge \$MAX_RAW_RU ]]; then
            USE_FLYE=true
        fi
    fi

    if [[ "\$USE_FLYE" == "true" ]]; then
        FLYE_CONF=\$(cat flye_confidence)
        echo -e "\${FLYE_RU}\\t\${FLYE_CONF}\\tflye" > ${sample_id}.${chromosome}.${allele}.ru_count.txt
    else
        # Use longest raw read as reference
        samtools fasta -@ ${task.cpus} allele.bam \
            | seqkit sort -l -r | seqkit head -n 1 > raw_read_ref.fasta

        HAS_RAW=false
        if grep -v '^>' raw_read_ref.fasta | grep -q '[ACGTacgt]'; then
            HAS_RAW=true
        fi

        if [[ "\$HAS_RAW" == "true" ]]; then
            makeblastdb -in raw_read_ref.fasta -dbtype nucl -out rawread_db
            blastn -db rawread_db -query ${query_d4z4} -outfmt 6 > rawread_d4z4.blast

            python3 - "${sample_id}" "${chromosome}" "${allele}" << 'PYEOF3'
import sys
import pandas as pd

sample_id, chrom, allele = sys.argv[1], sys.argv[2], sys.argv[3]
blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

try:
    d4z4 = pd.read_csv('rawread_d4z4.blast', sep='\\t', header=None, names=blast_cols)
    ru_count = len(d4z4[d4z4['length'] > 1500])
except pd.errors.EmptyDataError:
    ru_count = 0

with open(f'{sample_id}.{chrom}.{allele}.ru_count.txt', 'w') as f:
    f.write(f'{ru_count}\\tLow\\traw-read\\n')
PYEOF3
        else
            echo -e ">=\${MAX_RAW_RU}\\tLow\\traw-read" > ${sample_id}.${chromosome}.${allele}.ru_count.txt
        fi
    fi
    """
}


// ============================================================================
// Step 8: Per-allele CpG methylation
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
        --combine-mods \
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
// Step 9: Generate report
// ============================================================================
process REPORT {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    path allele_assignments
    path strategy_json
    path hp1_bed
    path hp2_bed
    path ungrouped_bed
    path readlevel_ru_files
    path assembly_ru_files
    val  sample_id

    output:
    path "${sample_id}.repeats.tsv",      emit: repeats
    path "${sample_id}.methylation.tsv",  emit: methylation

    script:
    """
    python3 - "${sample_id}" << 'PYEOF'
import sys
import json
import pandas as pd

sample_id = sys.argv[1]

df = pd.read_csv(f'{sample_id}.allele_assignments.tsv', sep='\\t')
with open(f'{sample_id}.allele_strategy.json') as f:
    strategy = json.load(f)

# --- Helpers ---
def read_ru_file(path):
    try:
        with open(path) as rf:
            parts = rf.read().strip().split('\\t')
            ru = parts[0] if len(parts) > 0 else 'NA'
            confidence = parts[1] if len(parts) > 1 else 'Low'
            ref_source = parts[2] if len(parts) > 2 else 'unknown'
            return ru, confidence, ref_source
    except (FileNotFoundError, IndexError):
        return 'NA', 'NA', 'unknown'

def majority_haplotype(series):
    for h in series.value_counts().index:
        if h != '-':
            return h
    return '-'

# --- repeats.tsv (grouped by Haplotype and Completeness within each allele) ---
def completeness_label(row):
    if row['Has_Head'] and row['Has_Tail']:
        return 'complete'
    elif row['Has_Head']:
        return 'head-only'
    elif row['Has_Tail']:
        return 'tail-only'
    return 'incomplete'

rows = []
for chrom in ['chr4', 'chr10']:
    chrom_strat = strategy.get(chrom, {})
    chrom_reads = df[df['Chromosome'] == chrom]

    if len(chrom_reads) == 0:
        continue

    for allele in sorted(chrom_reads['Allele'].unique()):
        allele_reads = chrom_reads[chrom_reads['Allele'] == allele].copy()
        if len(allele_reads) == 0:
            continue

        allele_strat = chrom_strat.get(allele)
        if allele_strat == 'read-level':
            _, _, ref_source = read_ru_file(
                f'{sample_id}.{chrom}.{allele}.readlevel_ru_count.txt')
            count_method = 'read-level'
        elif allele_strat == 'assembly':
            _, _, ref_source = read_ru_file(
                f'{sample_id}.{chrom}.{allele}.ru_count.txt')
            count_method = 'assembly'
        else:
            count_method = '-'
            ref_source = '-'

        allele_label = '-' if allele == 'unphased' else allele
        allele_reads['Complete'] = allele_reads.apply(completeness_label, axis=1)

        for (haplotype, complete), grp in allele_reads.groupby(['Haplotype', 'Complete']):
            rows.append({
                'Chromosome': chrom,
                'Allele': allele_label,
                'Haplotype': haplotype,
                'RU_Count': int(grp['RU_Count'].max()),
                'Complete': complete,
                'Count_Method': count_method,
                'Reference_Source': ref_source,
                'Read_Count': len(grp),
            })

repeats = pd.DataFrame(rows, columns=['Chromosome', 'Allele', 'Haplotype', 'RU_Count',
                                        'Complete', 'Count_Method', 'Reference_Source',
                                        'Read_Count'])
repeats.to_csv(f'{sample_id}.repeats.tsv', sep='\\t', index=False)

# --- methylation.tsv ---
REGIONS = {
    'DR1':      (193541000, 193541600),
    'DUX4-PAS': (193543430, 193543634),
}
# Reference atlas: ONT WGS, ~45-49 healthy controls, Gaussian-smoothed (sigma=150bp)
ATLAS_REF = {
    'DR1':      (54.70, 1.97),
    'DUX4-PAS': (51.81, 0.23),
}
# Caputo et al. 2022 (Cells 11:4114) — diagnostic thresholds from 150 controls
# DR1_CpG1 <= 53% and DUX4-PAS_CpG6 <= 78% indicate FSHD-compatible hypomethylation
LITERATURE_REF = {
    'DR1':      'CpG1 <= 53%',
    'DUX4-PAS': 'CpG6 <= 78%',
}

def read_bedmethyl(path):
    try:
        bm = pd.read_csv(path, sep='\\t', header=None, comment='#')
        return bm.rename(columns={0: 'chrom', 1: 'start', 9: 'N_valid', 10: 'pct_modified'})
    except (pd.errors.EmptyDataError, pd.errors.ParserError, FileNotFoundError):
        return pd.DataFrame(columns=['chrom', 'start', 'N_valid', 'pct_modified'])

bed_hp1 = read_bedmethyl(f'{sample_id}.methylation_HP1.bed')
bed_hp2 = read_bedmethyl(f'{sample_id}.methylation_HP2.bed')
bed_ung = read_bedmethyl(f'{sample_id}.methylation_ungrouped.bed')

# Combined: weighted average across HP1 + HP2 + ungrouped at each CpG
all_beds = pd.concat([bed_hp1, bed_hp2, bed_ung], ignore_index=True)
if len(all_beds) > 0 and 'N_valid' in all_beds.columns and all_beds['N_valid'].notna().any():
    all_beds['n_modified'] = all_beds['N_valid'] * all_beds['pct_modified'] / 100.0
    bed_combined = all_beds.groupby(['chrom', 'start']).agg(
        N_valid=('N_valid', 'sum'),
        n_modified=('n_modified', 'sum')
    ).reset_index()
    bed_combined['pct_modified'] = bed_combined['n_modified'] / bed_combined['N_valid'] * 100.0
else:
    bed_combined = pd.DataFrame(columns=['chrom', 'start', 'N_valid', 'pct_modified'])

meth_rows = []
for hp_label, meth in [('HP1', bed_hp1), ('HP2', bed_hp2), ('Combined', bed_combined)]:
    for region, (rstart, rend) in REGIONS.items():
        sub = meth[(meth['chrom'] == 'chr4') &
                   (meth['start'] >= rstart) &
                   (meth['start'] < rend)]
        mean_pct = round(float(sub['pct_modified'].mean()), 2) if len(sub) > 0 else 'NA'
        n_cpgs = len(sub)
        atlas_mean, atlas_sd = ATLAS_REF[region]
        meth_rows.append({
            'Haplotype': hp_label,
            'Region': region,
            'Methylation': mean_pct,
            'N_CpG': n_cpgs,
            'Ref_Atlas (Mean (SD) %)': f'{atlas_mean}% ({atlas_sd}%)',
            'Lit_Caputo (%)': LITERATURE_REF[region],
        })

pd.DataFrame(meth_rows, columns=['Haplotype', 'Region', 'Methylation', 'N_CpG',
                                   'Ref_Atlas (Mean (SD) %)', 'Lit_Caputo (%)']).to_csv(
    f'{sample_id}.methylation.tsv', sep='\\t', index=False)
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

    def idx_ext    = input_file.name.endsWith('.cram') ? '.crai' : '.bai'
    def input_idx  = file("${params.input}${idx_ext}")

    // Step 1: Extract FSHD-region reads
    EXTRACT_FSHD_READS(input_file, input_idx, ref_hg38, sample_id)

    // Step 2: Realign to T2T (HP/PS tags preserved)
    ALIGN_T2T(EXTRACT_FSHD_READS.out.bam, ref_t2t, sample_id)

    // Step 3: Extract D4Z4 candidates
    EXTRACT_D4Z4(ALIGN_T2T.out.bam, ALIGN_T2T.out.bai, sample_id)

    // Step 4: BLAST annotation
    BLAST_ANNOTATE(
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai, sample_id,
        file(params.query_d4z4), file(params.query_plam_polya),
        file(params.query_p13e11)
    )

    // Step 5: Allele assignment (HP tags from T2T BAM, alignment from D4Z4 BAM)
    ALLELE_ASSIGN(
        BLAST_ANNOTATE.out.tsv,
        ALIGN_T2T.out.bam, ALIGN_T2T.out.bai,
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai,
        sample_id
    )

    // Parse strategy JSON → route alleles to read-level or assembly
    allele_routing = ALLELE_ASSIGN.out.strategy
        .map { json_file ->
            new groovy.json.JsonSlurper().parse(json_file)
        }
        .flatMap { strat ->
            def tuples = []
            strat.each { chr, alleles ->
                alleles.each { allele, strategy ->
                    if (strategy) tuples << [chr, allele, strategy]
                }
            }
            return tuples
        }
        .branch { tuple ->
            read_level: tuple[2] == 'read-level'
            assembly:   tuple[2] == 'assembly'
        }

    // Step 6: Read-level counting (conditional — only alleles with complete reads)
    READLEVEL_COUNT(
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai,
        ALLELE_ASSIGN.out.tsv,
        file(params.query_d4z4), file(params.query_plam_polya),
        file(params.query_p13e11),
        sample_id,
        allele_routing.read_level.map { tuple -> tuple[0] },
        allele_routing.read_level.map { tuple -> tuple[1] }
    )

    // Step 7: Assembly-based counting (conditional — only alleles without complete reads)
    ASSEMBLY_COUNT(
        EXTRACT_D4Z4.out.bam, EXTRACT_D4Z4.out.bai,
        ALLELE_ASSIGN.out.tsv,
        file(params.query_d4z4), file(params.query_plam_polya),
        file(params.query_p13e11),
        sample_id,
        allele_routing.assembly.map { tuple -> tuple[0] },
        allele_routing.assembly.map { tuple -> tuple[1] }
    )

    // Collect outputs — .ifEmpty([]) prevents REPORT from hanging
    readlevel_ru_collected = READLEVEL_COUNT.out.ru_count.collect().ifEmpty([])
    assembly_ru_collected  = ASSEMBLY_COUNT.out.ru_count.collect().ifEmpty([])

    // Step 8: Methylation (HP tags from T2T BAM)
    METHYLATION(
        ALIGN_T2T.out.bam, ALIGN_T2T.out.bai,
        ref_t2t, ref_t2t_fai, sample_id
    )

    // Step 9: Report
    REPORT(
        ALLELE_ASSIGN.out.tsv,
        ALLELE_ASSIGN.out.strategy,
        METHYLATION.out.hp1_bed,
        METHYLATION.out.hp2_bed,
        METHYLATION.out.ungrouped_bed,
        readlevel_ru_collected,
        assembly_ru_collected,
        sample_id
    )
}
