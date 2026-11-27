"""Biological consistency checks and the validation report."""

from collections.abc import Sequence

from fshdetector import annotate
from fshdetector import layout as layout_module

GENERIC_HAPLOTYPES = frozenset({"4q", "10q", "Unknown"})


def check(
    calls: Sequence[annotate.ReadCall],
) -> tuple[list[tuple[str, str]], dict[str, list[str]]]:
    """Checks the calls of one sample.

    Returns:
        A list of (level, message) pairs about haplotype counts, and the
        warning flags of each read that has any.
    """
    messages = []
    flags: dict[str, list[str]] = {}

    def flag(call: annotate.ReadCall, name: str) -> None:
        flags.setdefault(call.read_id, []).append(name)

    for chrom in layout_module.CHROMOSOMES:
        on_chrom = [c for c in calls if c.chromosome == chrom]
        if not on_chrom:
            continue
        typed = sorted({c.haplotype for c in on_chrom} - GENERIC_HAPLOTYPES)
        level = "WARN" if len(typed) > 2 else "PASS"
        messages.append(
            (level, f"{chrom} has {len(typed)} typed haplotype(s): {typed}")
        )

        for haplotype in typed:
            group = [c for c in on_chrom if c.haplotype == haplotype]
            complete = [c for c in group if c.complete]
            if not complete:
                continue
            shortest = min(c.read_length for c in complete)
            largest = max(c.ru_count for c in complete)
            for call in group:
                if call.complete:
                    continue
                if call.read_length >= shortest:
                    flag(call, "longer_than_complete")
                if call.ru_count > largest:
                    flag(call, "RU_exceeds_complete")

        alleles = [
            c
            for c in on_chrom
            if c.complete and c.haplotype not in GENERIC_HAPLOTYPES
        ]
        if alleles:
            largest = max(c.ru_count for c in alleles)
            for call in on_chrom:
                if (
                    call.haplotype in GENERIC_HAPLOTYPES
                    and call.ru_count > largest
                ):
                    flag(call, "RU_exceeds_all_alleles")

    for call in calls:
        if call.chimeric:
            flag(call, "chimeric")
        if call.mixed_strand:
            flag(call, "mixed_strand")
    return messages, flags


def report(sample: str, calls: Sequence[annotate.ReadCall]) -> str:
    """Formats the validation report as text."""
    messages, flags = check(calls)
    complete = sum(c.complete for c in calls)
    head = sum(c.has_head and not c.has_tail for c in calls)
    tail = sum(c.has_tail and not c.has_head for c in calls)
    lines = [
        f"=== Biological Validation: {sample} ===",
        "",
        "Summary:",
        f"  Total D4Z4 reads: {len(calls)}",
        f"  Complete (head+tail): {complete}",
        f"  Head only: {head}",
        f"  Tail only: {tail}",
        f"  Fragment (no anchors): {len(calls) - complete - head - tail}",
        "",
        "Haplotype counts per chromosome:",
    ]
    lines += [f"  [{level}] {message}" for level, message in messages]
    lines += ["", "Per-allele summary:"]

    groups = []
    for chrom in layout_module.CHROMOSOMES:
        on_chrom = [c for c in calls if c.chromosome == chrom]
        for haplotype in sorted(
            {c.haplotype for c in on_chrom} - GENERIC_HAPLOTYPES
        ):
            group = [c for c in on_chrom if c.haplotype == haplotype]
            groups.append((chrom, haplotype, group))

    for chrom, haplotype, group in groups:
        full = [c for c in group if c.complete]
        permissive = sum(c.permissive for c in group)
        chimeric = sum(c.chimeric for c in group)
        source = "complete" if full else "incomplete"
        lines += [
            "",
            f"  {chrom}/{haplotype}:",
            f"    Permissive reads: {permissive} of {len(group)}",
            f"    Total reads: {len(group)} ({len(full)} complete, "
            f"{len(group) - len(full)} incomplete)",
            f"    Max RU: {max(c.ru_count for c in full or group)} "
            f"(from {source} reads)",
        ]
        if chimeric:
            lines.append(f"    Chimeric reads: {chimeric}")

    lines += ["", "Complete reads (Has_Head AND Has_Tail):"]
    found = False
    for chrom, haplotype, group in groups:
        full = [c for c in group if c.complete]
        if not full:
            continue
        found = True
        permissive = sum(c.permissive for c in full)
        lines += [
            "",
            f"  {chrom}/{haplotype} (permissive {permissive} of "
            f"{len(full)}): {len(full)} complete read(s)",
        ]
        for call in sorted(full, key=lambda c: (-c.ru_count, c.read_id)):
            chimeric = "  CHIMERIC" if call.chimeric else ""
            lines += [
                f"    {call.read_id}",
                f"      RU={call.ru_count}  Read_Length={call.read_length}  "
                f"Chr_Confidence={call.chr_confidence}",
                f"      Proximal_Marker={call.proximal_marker or '-'}  "
                f"Distal_Marker={call.distal_marker or '-'}{chimeric}",
                f"      BlnI_Chr4={call.blni_chr4}  "
                f"BlnI_Chr10={call.blni_chr10}",
            ]
    if not found:
        lines.append("  No complete reads found.")
    lines.append("")

    flagged = [c for c in calls if c.read_id in flags]
    if flagged:
        lines += [
            f"Low confidence reads ({len(flagged)}):",
            "  Possible reasons: longer_than_complete, RU_exceeds_complete,",
            "                    RU_exceeds_all_alleles, chimeric,"
            " mixed_strand",
            "",
        ]
        for call in flagged:
            lines += [
                f"  {call.read_id}",
                f"    {call.chromosome}  {call.haplotype}  "
                f"RU={call.ru_count}  len={call.read_length}",
                f"    reason: {'; '.join(flags[call.read_id])}",
            ]
    else:
        lines.append("No low confidence reads.")
    return "\n".join(lines) + "\n"
