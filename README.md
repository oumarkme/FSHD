# nf-fshd: ONT D4Z4 Repeat Unit Counting, Haplotyping, and Methylation Pipeline

A Nextflow DSL2 pipeline for FSHD diagnosis from hg38-aligned ONT BAM/CRAM files with HP tags from upstream phasing. The pipeline resolves all four D4Z4 alleles (2× chr4, 2× chr10), determines 4qA/4qB haplotype, counts D4Z4 repeat units per allele, and reports allele-specific CpG methylation at the distal DUX4 region. All analysis coordinates reference **T2T-CHM13v2.0**.

---

## Biological background

### D4Z4 and FSHD

The D4Z4 macrosatellite array on chromosome 4q35 consists of tandem 3.3 kb repeat units. Healthy individuals carry 11–100+ RU per allele. FSHD1 occurs when one allele contracts to 1–10 RU **on a permissive haplotype background** (4qA or 4qA-L), which stabilizes DUX4 mRNA via a polyadenylation signal (PAS) in the pLAM region distal to the array (Lemmers et al., 2012).

### Haplotype backgrounds

| Haplotype | Permissive? | Key feature | Population freq |
|-----------|-------------|-------------|-----------------|
| **4qA** (includes 4A161) | Yes | Contains pLAM with functional DUX4 PAS (ATTAAA) | ~50% of chr4 alleles |
| **4qA-L** (4A161L) | Yes | 1.6 kb longer distal partial D4Z4 unit; functional PAS | ~20% of European 4A161 alleles |
| **4qB** | No | Lacks pLAM entirely; DUX4 mRNA unstable | ~50% of chr4 alleles |
| **10qA** | No | 98.5% identity to 4qA but non-functional PAS (SNPs) | Most chr10 alleles |

Only contractions on 4qA or 4qA-L cause FSHD1. Contractions on 4qB are clinically silent (Lemmers et al., 2012).

### D4Z4 array sizes

| RU count | Array size | Clinical significance |
|----------|-----------|----------------------|
| 1–3 | 3.3–10 kb | Severe/infantile FSHD1 (high penetrance) |
| 4–7 | 13–23 kb | Classic FSHD1 (variable expressivity) |
| 8–10 | 26–33 kb | Borderline; reduced penetrance (~40–50%); requires methylation |
| 11–20 | 36–66 kb | Normal; may contribute to FSHD2 if modifier present |
| 21–100+ | 69–330+ kb | Normal/healthy |

---

## Technical constraints (ONT sequencing)

### Read length vs D4Z4 spanning

| Parameter | Standard ligation (LSK114, with size selection) | Ultra-long (ULK114) |
|-----------|------------------------------------------------|---------------------|
| Typical N50 | 25–30 kb | 50–100+ kb |
| Max reliable spanning | ~7–8 RU (~26 kb) | ~15–20 RU |
| Marginal spanning | 9–10 RU (~33 kb) | 20–30 RU |
| Never spanned | ≥11 RU (≥36 kb) | ≥30 RU |

A read is **"complete"** when it contains both the proximal p13E-11 anchor AND the distal pLAM anchor. With standard N50 ~25–30 kb, complete reads can only come from alleles with ≤~10 RU. Healthy alleles (≥11 RU) will never produce complete reads with standard prep.

### Coverage per allele

| Genome coverage | Reads per allele (theoretical) | Effective (subtelomeric penalty) |
|-----------------|-------------------------------|----------------------------------|
| 30× | ~7–8× | ~4–6× |
| 60× | ~15× | ~10–12× |

### HP tag reliability

- HP tags are reliable in unique flanking regions adjacent to D4Z4
- Within the D4Z4 array: very few heterozygous SNPs → phasing may fail
- Most D4Z4 reads carry HP tags inherited from flanking-region phasing, but some reads (especially those entirely within the array) may lack HP tags

### Chr4 vs chr10 discrimination

- D4Z4 arrays share ~98.5–99% identity between chr4 and chr10
- Discrimination relies on unique flanking sequences and a chr10-specific BlnI restriction site within D4Z4 units
- Reads must extend into flanking regions for reliable chromosome assignment
- Haplotype detection (4qA vs 10qA) requires reads extending past the last D4Z4 unit to detect pLAM

### Flye assembly caveats

- Minimum ~10× coverage; 30–40× recommended. Per-allele coverage at D4Z4 is typically 4–12× (marginal)
- Flye may collapse identical D4Z4 copies, underestimating RU count (Kolmogorov et al., 2019)
- Resolution: take whichever has MORE RUs — assembly or longest raw read

---

## Algorithm: complete-reads-first allele assignment

The pipeline processes chr4 and chr10 independently. After BLAST annotation, every read has: `RU_Count`, `Has_Head` (p13E-11), `Has_Tail` (pLAM), `Haplotype` (4qA/10qA/−), and `HP` tag (1/2/none).

The algorithm has three sequential phases:

```
Phase 1: SEPARATE ALLELES     → which reads belong together?
Phase 2: FIND COMPLETE READS  → which allele has spanning reads?
Phase 3: COUNT per allele     → read-level (consensus) or assembly?
```

### Phase 1: Separate alleles

Group reads into two alleles using the best available signal:

```
IF any reads on this chromosome have HP tags:
    Allele_1 = reads with HP == 1
    Allele_2 = reads with HP == 2
    unphased = reads without HP tag
    separation_method = 'hp'

ELSE (no HP tags at all):
    IF complete reads exist:
        complete_ru = median(RU_Count of complete reads)
        Allele_1 = reads with RU_Count <= complete_ru + 2
        Allele_2 = reads with RU_Count > complete_ru + 2
        separation_method = 'ru-anchor'
    ELSE:
        IF std(RU_counts) >= 1.5 AND len(reads) >= 2:
            Jenks 2-class split -> Allele_1 (shorter), Allele_2 (longer)
            separation_method = 'jenks'
        ELSE:
            Allele_1 = all reads
            Allele_2 = empty
            separation_method = 'single'
```

HP is always preferred when available — it uses flanking-region SNPs and doesn't depend on RU count differences.

### Phase 2: Assign unphased reads

Complete reads (Has_Head AND Has_Tail) span the full array and **can only come from the shorter allele** (biological constraint: healthy alleles ≥11 RU exceed ONT N50; Huang et al., 2024).

```
IF no complete reads:
    Skip to Phase 3 (all alleles use assembly)

IF complete reads already have HP tags:
    Already in the correct allele from Phase 1

IF complete reads lack HP tags (separation was by HP):
    Assign to the allele whose median RU is closest

All remaining unphased reads:
    IF allele has complete reads (max_RU = X) AND read RU > X + 2:
        Force-assign to the OTHER allele (ru-constraint)
    ELSE:
        minimap2 align read against HAP1 pool and HAP2 pool
        Assign to allele with more matching bases (alignment)
        Tie: assign to allele with closer median RU (ru-proximity)
```

### Phase 3: Per-allele strategy

Each allele independently gets one of two strategies:

```
For each allele:

    IF allele has complete reads -> STRATEGY: READ-LEVEL
        spanning = complete + partial-distal reads
        raw_ref = longest spanning read

        IF >= 3 spanning reads:
            Racon consensus (raw_ref as seed, polish with others)
        ELSE:
            consensus = raw_ref as-is

        reference = consensus
        RU_count = BLAST d4z4.fa against consensus (hits > 1500 bp)

    ELSE -> STRATEGY: ASSEMBLY
        Run Flye on allele reads (--nano-hq --genome-size 300k)
        max_raw_ru = max(RU_Count) among reads in this allele

        IF Flye succeeded AND Flye_RU >= max_raw_ru:
            reference = Flye contig
        ELSE:
            reference = longest raw read

        RU_count = BLAST d4z4.fa against reference (hits > 1500 bp)
```

### Example walkthroughs

**Classic FSHD1 (most common):**

| Phase | What happens |
|-------|-------------|
| Phase 1 | HP tags exist → Allele_1 (HP=1), Allele_2 (HP=2) |
| Phase 2 | Complete reads found with HP=1 → already in Allele_1 |
| Phase 3 | Allele_1: read-level (Racon consensus). Allele_2: assembly (Flye vs max raw read) |

**No HP tags, complete reads exist:**

| Phase | What happens |
|-------|-------------|
| Phase 1 | No HP → use complete reads as RU anchor. Allele_1 (≤RU+2), Allele_2 (>RU+2) |
| Phase 2 | Complete reads already in Allele_1 by construction |
| Phase 3 | Allele_1: read-level. Allele_2: assembly (or empty → single allele) |

**Both alleles short (homozygous contraction):**

| Phase | What happens |
|-------|-------------|
| Phase 1 | HP tags → Allele_1 (HP=1), Allele_2 (HP=2) |
| Phase 2 | Complete reads in BOTH alleles |
| Phase 3 | Both alleles: read-level (Racon consensus each) |

**No complete reads (chr10, or healthy chr4):**

| Phase | What happens |
|-------|-------------|
| Phase 1 | HP tags → Allele_1 (HP=1), Allele_2 (HP=2) |
| Phase 2 | No complete reads → skip |
| Phase 3 | Both alleles: assembly (Flye vs max raw read) |

---

## Pipeline steps

### Step 1: Extract FSHD-region reads (EXTRACT_FSHD_READS)

Extract reads from FSHD-associated regions in hg38 coordinates plus unmapped reads.

```
BED (hg38):
    chr4:189000000-190500000     # 4q35 D4Z4 flanking
    chr10:133000000-134000000    # 10q26 D4Z4 flanking

samtools view -b -h -M -L regions.bed -> target.bam
samtools view -b -f 4 -> unmapped.bam
samtools merge -> subset.bam
```

### Step 2: Realign to T2T-CHM13v2.0 (ALIGN_T2T)

Realign to T2T, preserving MM, ML, HP, PS tags via `samtools fastq -T` and `minimap2 -y`.

`minimap2 -ax lr:hq -y --eqx --secondary=no | samtools sort`

### Step 3: Extract D4Z4 candidate reads (EXTRACT_D4Z4)

```
BED (T2T):
    chr4:193230000-193610000
    chr10:134560000-134800000
```

Published output: `{SAMPLE_ID}.d4z4.bam` + `.bai`

### Step 4: Chromosome assignment and feature annotation (BLAST_ANNOTATE)

BLAST-based feature annotation following the approach of Lowenstern et al. (2026) and Huang et al. (2024):

```
blastn d4z4.fa       -> RU counting (hits > 1500 bp)
blastn pLAM-polyA.fa -> distal anchor + 4qA/10qA haplotype
blastn P13E11.fasta  -> proximal anchor
```

**Chromosome assignment (three-tier):**

1. **BAM alignment** (initial): minimap2 alignment to T2T reference
2. **BlnI motif counting** (correction): chr4-type `CTAGGTCTAGG` vs chr10-type `CTAGGCCTAGG` — majority wins
3. **Haplotype override** (final): 4qA (ATTAAA) → force chr4; 10qA (ATCAAA) → force chr10

### Step 5: Allele assignment (ALLELE_ASSIGN)

Three-phase algorithm: Phase 1 (separate by HP/RU-anchor/Jenks), Phase 2 (assign unphased reads via RU-constraint + alignment + proximity), Phase 3 (determine read-level vs assembly strategy per allele).

### Step 6: Read-level counting (READLEVEL_COUNT)

For alleles with complete reads: extract spanning reads → Racon consensus (Vaser et al., 2017; if ≥3 reads) → BLAST RU count. This approach follows D4Z4End2End (Xiao et al., 2025).

### Step 7: Assembly-based counting (ASSEMBLY_COUNT)

For alleles without complete reads: Flye assembly → compare Flye RU with max raw-read RU → select reference with more RUs → BLAST RU count.

### Step 8: Per-allele CpG methylation (METHYLATION)

```
modkit pileup --phased --cpg --combine-strands --combine-mods
    --modified-bases 5mC
    --include-bed chr4:193540172-193543634, chr10:134738000-134741148
```

Runs independently; only needs T2T BAM from Step 2.

### Step 9: Generate report (REPORT)

Produces `repeats.tsv` and `methylation.tsv` from allele assignments, RU counts, and methylation data.

---

## Pipeline flow diagram

```
Input: hg38 BAM/CRAM (HP/PS + MM/ML tags)
    |
    v
Step 1: EXTRACT_FSHD_READS (hg38 regions)
    |
    v
Step 2: ALIGN_T2T (preserve MM/ML/HP/PS) -------.
    |                                            |
    v                                            |
Step 3: EXTRACT_D4Z4                             |
    |                                            |
    v                                            |
Step 4: BLAST_ANNOTATE                           |
(RU, anchors, BlnI chr, haplotype)               |
    |                                            |
    v                                            |
Step 5: ALLELE_ASSIGN  <-------------------------'
    |
    |  Phase 1: Separate alleles
    |      HP tags? -> group by HP
    |      No HP?   -> RU-anchor (if complete reads) or Jenks
    |
    |  Phase 2: Assign complete reads
    |      Already have HP? -> done
    |      No HP? -> assign to allele with closest median RU
    |
    |  Phase 3: Per-allele strategy
    |      Has complete reads? -> Step 6 (read-level)
    |      No complete reads?  -> Step 7 (assembly)
    |
    |               .-----------------.
    v               v                 |
Step 6:         Step 7:               |
READLEVEL       ASSEMBLY              |
(Racon          (Flye vs raw read     |
 consensus)      -> pick more RUs)    |
    |               |                 |
    '-------+-------'                 |
            |                         |
            |         .---------------'
            v         v
        Step 8: METHYLATION
        (modkit --phased from T2T BAM)
            |
            v
        Step 9: REPORT
        (repeats.tsv + methylation.tsv)
```

---

## Output

Four published files per sample:

| File | Description |
|------|-------------|
| `{SAMPLE_ID}.d4z4.bam` | D4Z4 candidate reads aligned to T2T-CHM13v2.0 |
| `{SAMPLE_ID}.d4z4.bam.bai` | BAM index |
| `{SAMPLE_ID}.repeats.tsv` | D4Z4 repeat unit counts per allele |
| `{SAMPLE_ID}.methylation.tsv` | CpG methylation at distal D4Z4/DUX4 region |

### repeats.tsv

One row per (Allele, Haplotype, Complete) group. RU_Count is the max within each group.

| Column | Description |
|--------|-------------|
| Chromosome | chr4 or chr10 |
| Allele | HAP1 / HAP2 / − (unphased) |
| Haplotype | 4qA / 10qA / − (no pLAM = 4qB or untyped) |
| RU_Count | Max RU count among reads in this group |
| Complete | complete / head-only / tail-only / incomplete |
| Count_Method | read-level / assembly / − (unphased) |
| Reference_Source | consensus / flye / raw-read / − (unphased) |
| Read_Count | Number of reads in this group |

### methylation.tsv

| Column | Description |
|--------|-------------|
| Haplotype | HP1, HP2, or Combined |
| Region | DR1 or DUX4-PAS |
| Methylation | Mean % modified across CpGs in region |
| N_CpG | Number of CpGs with data |
| Ref_Atlas (Mean (SD) %) | ONT methylation reference atlas (~45–49 healthy controls; Huang et al., 2024) |
| Lit_Caputo (%) | Diagnostic threshold (Caputo et al., 2022; Strafella et al., 2024) |

**Reference values (T2T coordinates):**

| Region | T2T coordinates | Ref_Atlas | Lit_Caputo threshold |
|--------|----------------|-----------|---------------------|
| DR1 | chr4:193,541,000–193,541,600 | 54.70% (1.97%) | CpG1 ≤ 53% |
| DUX4-PAS | chr4:193,543,430–193,543,634 | 51.81% (0.23%) | CpG6 ≤ 78% |

---

## Usage

```bash
bash nf-fshd/run.sh /path/to/sample.bam
```

The `run.sh` script auto-detects the environment, builds the Apptainer container if needed, and launches Nextflow with `-resume`.

**Input:** hg38-aligned BAM/CRAM with HP/PS tags (from any phasing tool) and MM/ML tags (from basecalling).

---

## Upstream prerequisites

These steps are performed outside this pipeline. The pipeline input is an hg38-aligned BAM/CRAM carrying both methylation tags (MM/ML) from basecalling and haplotype tags (HP/PS) from phasing. Any phasing tool that writes standard HP/PS tags is compatible.

| Step | Tool | Version | Output | Reference |
|------|------|---------|--------|-----------|
| Basecalling | Dorado (SUP + `--modified-bases 5mCG_5hmCG`) | ≥0.7.0 | Unaligned modBAM with MM/ML tags | Oxford Nanopore Technologies (2024b) |
| Genome-wide alignment | minimap2 (`-ax lr:hq`) | ≥2.28 | hg38-aligned BAM | Li (2018) |
| Variant calling | e.g., Clair3 | — | VCF (SNPs for phasing) | Zheng et al. (2022) |
| Phasing + haplotagging | e.g., LongPhase, WhatsHap | — | BAM with HP/PS tags | Lin et al. (2022); Martin et al. (2023) |

---

## Reference genome and query sequences

**Reference:** T2T-CHM13v2.0 (Nurk et al., 2022; `GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`)

**Key T2T-CHM13v2.0 coordinates:**

| Region | Coordinates | Notes |
|--------|------------|-------|
| Chr4 D4Z4 array | chr4:193,282,466–193,558,261 | ~276 kb; 33 RUs in reference |
| Chr4 DUX4 methylation target | chr4:193,540,172–193,543,634 | Distal D4Z4 + DUX4 promoter/coding |
| Chr4 p13E-11 anchor | chr4:192,667,301–192,902,247 | Upstream unique anchoring region |
| Chr10 D4Z4 array | chr10:134,615,120–134,741,148 | ~126 kb |

**Query sequences (`nf-fshd/query/`):**

| File | Purpose |
|------|---------|
| `d4z4.fa` | Single D4Z4 RU (~3.3 kb) for RU counting |
| `pLAM-polyA.fa` | hapA (4qA ATTAAA PAS) / hapB (10qA ATCAAA PAS) for haplotype + distal anchor; 4qB has no pLAM |
| `P13E11.fasta` | p13E-11 probe for proximal anchor |

---

## Software

### Container (`environment.sif`, built from `environment.def`)

| Tool | Version | Purpose | Reference |
|------|---------|---------|-----------|
| Python | 3.12 | Scripting runtime | |
| samtools | 1.23 | BAM/CRAM extraction, sort, index, merge, FASTQ conversion | Danecek et al. (2021) |
| minimap2 | 2.30 | T2T realignment; alignment-based unphased read assignment | Li (2018) |
| BLAST+ | 2.17.0 | RU counting, anchor detection, haplotyping, chromosome assignment | Camacho et al. (2009) |
| modkit | 0.6.1 | Phased CpG methylation pileup | Oxford Nanopore Technologies (2024a) |
| Flye | 2.9.6 | Per-allele local assembly | Kolmogorov et al. (2019) |
| Racon | 1.5.0 | Consensus polishing from spanning reads | Vaser et al. (2017) |
| seqkit | 2.13.0 | FASTA extraction and sorting | Shen et al. (2016) |
| pysam | 0.23.3 | BAM read iteration and HP tag extraction in Python | |
| pandas | 3.0.1 | Tabular data processing and report generation | |
| jenkspy | 0.4.1 | Jenks natural breaks for fallback allele separation | |

### Runtime requirements (host)

| Tool | Purpose | Reference |
|------|---------|-----------|
| Nextflow | Workflow engine (DSL2) | Di Tommaso et al. (2017) |
| Apptainer | Container runtime for `environment.sif` | Kurtzer et al. (2017) |

---

## Known limitations

1. **Read length constraint.** At N50 ~25–30 kb, complete reads are limited to alleles ≤~10 RU. Ultra-long prep extends this range but is not standard.

2. **Subtelomeric coverage.** Effective per-allele coverage at D4Z4 is 30–50% lower than genome-wide average. At 30× WGS, expect ~4–6 reads per allele.

3. **Flye tandem repeat collapse.** Flye may underestimate RU count for large arrays. The pipeline mitigates this by comparing Flye RU with max single-read RU and selecting the larger.

4. **Methylation coverage.** At 30× WGS, per-allele methylation coverage may be only 3–5 reads. ≥60× recommended for clinical-grade calls.

5. **CHM13 reference gap for 4qA.** T2T-CHM13v2.0 carries 4qB on chr4 (Nurk et al., 2022). 4qA-specific pLAM/PAS is absent from the reference. BLAST-based annotation detects pLAM in read sequence directly.

6. **Somatic mosaicism.** Mosaic samples may show multiple RU count clusters from the same allele. The pipeline does not explicitly detect mosaicism but unusual clustering patterns should be flagged for review.

7. **Cis duplications.** Two adjacent D4Z4 arrays on one allele (~1.5% of Europeans) may produce two distinct RU groups in one HP. The pipeline may misinterpret this as two alleles.

---

## References

### FSHD biology and clinical

Caputo, V., Megalizzi, D., Fabrizio, C., Termine, A., Colantoni, L., Bax, C., Gimenez, J., Monforte, M., Tasca, G., Ricci, E., Caltagirone, C., Giardina, E., Cascella, R., & Strafella, C. (2022). D4Z4 methylation levels combined with a machine learning pipeline highlight single CpG sites as discriminating biomarkers for FSHD patients. *Cells*, *11*(24), Article 4114. https://doi.org/10.3390/cells11244114

Lemmers, R. J. L. F., Tawil, R., Petek, L. M., Balog, J., Block, G. J., Santen, G. W. E., Amell, A. M., van der Vliet, P. J., Almomani, R., Straasheijm, K. R., Krom, Y. D., Klooster, R., Sun, Y., den Dunnen, J. T., Helmer, Q., Donlin-Smith, C. M., Padberg, G. W., van Engelen, B. G. M., de Greef, J. C., ... van der Maarel, S. M. (2012). Digenic inheritance of an SMCHD1 mutation and an FSHD-permissive D4Z4 allele causes facioscapulohumeral muscular dystrophy type 2. *Nature Genetics*, *44*(12), 1370–1374. https://doi.org/10.1038/ng.2454

Strafella, C., Megalizzi, D., Trastulli, G., Proietti Piorgo, E., Colantoni, L., Tasca, G., Monforte, M., Zampatti, S., Primiano, G., Sancricca, C., Bortolani, S., Torchia, E., Ravera, B., Torri, F., Gadaleta, G., Risi, B., Caria, F., Gerardi, F., Carraro, E., ... Giardina, E. (2024). Integrating D4Z4 methylation analysis into clinical practice: Improvement of FSHD molecular diagnosis through distinct thresholds for 4qA/4qA and 4qA/4qB patients. *Clinical Epigenetics*, *16*, Article 148. https://doi.org/10.1186/s13148-024-01747-2

### FSHD long-read pipelines

Huang, M., Zhang, Q., Jiao, J., Shi, J., Xu, Y., Zhang, C., Zhou, R., Liu, W., Liang, Y., Chen, H., Wang, Y., Xu, Z., & Hu, P. (2024). Comprehensive genetic analysis of facioscapulohumeral muscular dystrophy by Nanopore long-read whole-genome sequencing. *Journal of Translational Medicine*, *22*, Article 451. https://doi.org/10.1186/s12967-024-05259-8

Lowenstern, T., Madritsch, S., Horner, D., Brait, N., Lafci, N. G., Schachner, A., Bujalkova, M. G., Hengstschlager, M., Dremsek, P., & Laccone, F. (2026). DUCKS4: A comprehensive workflow for Nanopore sequencing analysis of facioscapulohumeral muscular dystrophy (FSHD). *Human Genomics*, *20*, Article 21. https://doi.org/10.1186/s40246-026-00921-2

Xiao, L. C., Semwal, A., St John, B., Zeglinski, K., Su, S., Lancaster, J., Xue, S., Reversade, B., Ritchie, M. E., Magdinier, F., Blewitt, M. E., & Gouil, Q. (2025). D4Z4End2End: Complete genetic and epigenetic architecture of D4Z4 macrosatellites in FSHD, BAMS and reference cohorts. *medRxiv*. https://doi.org/10.1101/2025.04.24.25326320

### Reference genome

Nurk, S., Koren, S., Rhie, A., Rautiainen, M., Bzikadze, A. V., Mikheenko, A., Vollger, M. R., Altemose, N., Uralsky, L., Gershman, A., Aganezov, S., Hoyt, S. J., Diekhans, M., Logsdon, G. A., Alonge, M., Antonarakis, S. E., Borchers, M., Bouffard, G. G., Brooks, S. Y., ... Phillippy, A. M. (2022). The complete sequence of a human genome. *Science*, *376*(6588), 44–53. https://doi.org/10.1126/science.abj6987

### Bioinformatics tools

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics*, *10*, Article 421. https://doi.org/10.1186/1471-2105-10-421

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, *10*(2), Article giab008. https://doi.org/10.1093/gigascience/giab008

Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, *35*(4), 316–319. https://doi.org/10.1038/nbt.3820

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, *37*(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8

Kurtzer, G. M., Sochat, V., & Bauer, M. W. (2017). Singularity: Scientific containers for mobility of compute. *PLOS ONE*, *12*(5), Article e0177459. https://doi.org/10.1371/journal.pone.0177459

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. *Bioinformatics*, *34*(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

Lin, J.-H., Chen, L.-C., Yu, S.-C., & Huang, Y.-T. (2022). LongPhase: An ultra-fast chromosome-scale phasing algorithm for small and large variants. *Bioinformatics*, *38*(7), 1816–1822. https://doi.org/10.1093/bioinformatics/btac058

Martin, M., Patterson, M., Garg, S., Fischer, S. O., Pisanti, N., Klau, G. W., Schonhuth, A., & Marschall, T. (2023). Read-based phasing and analysis of phased variants with WhatsHap. In *Methods in Molecular Biology* (Vol. 2590, pp. 127–138). Springer. https://doi.org/10.1007/978-1-0716-2819-5_8

Oxford Nanopore Technologies. (2024a). *Modkit* [Computer software]. https://github.com/nanoporetech/modkit

Oxford Nanopore Technologies. (2024b). *Dorado* [Computer software]. https://github.com/nanoporetech/dorado

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*, *11*(10), Article e0163962. https://doi.org/10.1371/journal.pone.0163962

Vaser, R., Sovic, I., Nagarajan, N., & Sikic, M. (2017). Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Research*, *27*(5), 737–746. https://doi.org/10.1101/gr.214270.116

Zheng, Z., Li, S., Su, J., Leung, A. W.-S., Lam, T.-W., & Luo, R. (2022). Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nature Computational Science*, *2*(12), 797–803. https://doi.org/10.1038/s43588-022-00387-x
