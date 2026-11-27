"""Turns the BLAST hits of one sample into read calls and barcode fragments."""

from collections.abc import Mapping, Sequence
import dataclasses

import numpy as np

from fshdetector import annotate
from fshdetector import assemble
from fshdetector import barcode
from fshdetector import blast
from fshdetector import layout as layout_module
from fshdetector import sequences as seqlib

# A site that reads of one allele often disagree on is not usable. The limit
# rises with the error rate of the sample.
MAX_SITE_ERROR = 0.05
MAX_RELATIVE_ERROR = 3.0  # Times the median site error.
MIN_COMPLETE_UNIT = 2800  # Aligned length of a complete repeat unit.


@dataclasses.dataclass
class SampleReads:
    """The analyzed reads of one sample.

    Attributes:
        calls: calls of every read with a repeat unit.
        units: repeat units of every read with a unit hit.
        usable: calls on chr4 or chr10 that can be placed on an array.
        sites: variant sites of the sample.
        fragments: read name to its assembly view (usable reads).
        placed_units: read name to the units in the rows of its fragment.
        misplaced: reads with an array-end marker on the wrong side.
        model: scoring model of the sample.
    """

    calls: list[annotate.ReadCall]
    units: dict[str, barcode.ReadUnits]
    usable: list[annotate.ReadCall]
    sites: barcode.Sites
    fragments: dict[str, assemble.Fragment]
    placed_units: dict[str, list]
    misplaced: set[str]
    model: barcode.Model

    def partial_by_chromosome(self) -> dict[str, list[assemble.Fragment]]:
        """Reads without both array ends, per chromosome."""
        return {
            chromosome: [
                self.fragments[c.read_id]
                for c in self.usable
                if c.chromosome == chromosome
                and not (
                    self.fragments[c.read_id].head
                    and self.fragments[c.read_id].tail
                )
            ]
            for chromosome in layout_module.CHROMOSOMES
        }


def make_fragment(
    call: annotate.ReadCall,
    read: barcode.ReadUnits,
    sites: barcode.Sites,
) -> tuple[assemble.Fragment, list[barcode.Unit | None]]:
    """Builds the assembly view of a read.

    The rows are the counted units of the read call. The partial units that
    flank an array never add to a contig, and at an anchored end the first (or
    last) row is the first (or last) unit of the array.
    """
    start, end = read.counted_range()
    units = read.units[start:end]
    trimmed = barcode.ReadUnits(
        read.read_id, read.reverse, read.read_length, units
    )
    kind = call.haplotype[-1] if call.haplotype[-1] in "AB" else ""
    fragment = assemble.Fragment(
        call.read_id,
        barcode.matrix(trimmed, sites),
        call.has_head,
        call.has_tail,
        kind,
        annotate.marker_chromosome(call.proximal_marker, call.distal_marker),
        # The first unit of an unanchored read and the last unit of any read
        # may be partial.
        any(
            u is None or u.length < MIN_COMPLETE_UNIT
            for u in units[0 if call.has_head else 1 : -1]
        ),
    )
    return fragment, units


def anchors_out_of_order(
    read: barcode.ReadUnits, anchors: Sequence[blast.Hit], proximal: bool
) -> bool:
    """True if every anchor hit lies on the wrong side of the unit array."""
    units = [u for u in read.units if u is not None]
    first = (units[0].start + units[0].end) / 2
    last = (units[-1].start + units[-1].end) / 2
    for hit in anchors:
        middle = (hit.sstart + hit.send) / 2
        if read.reverse:
            middle = read.read_length - middle
        if (middle < last) if proximal else (middle > first):
            return False
    return bool(anchors)


def _anchor_hits(
    hits: dict[str, list[blast.Hit]],
) -> tuple[dict[str, list], dict[str, list]]:
    proximal: dict[str, list[blast.Hit]] = {}
    distal: dict[str, list[blast.Hit]] = {}
    for hit in hits[blast.MARKERS]:
        if hit.identity < annotate.MIN_MARKER_IDENTITY:
            continue
        side = proximal if hit.query in annotate.PROXIMAL_MARKERS else distal
        side.setdefault(hit.read_id, []).append(hit)
    for hit in hits[blast.P13E11]:
        if (
            hit.identity >= annotate.MIN_MARKER_IDENTITY
            and hit.length >= annotate.MIN_P13E11_LENGTH
        ):
            proximal.setdefault(hit.read_id, []).append(hit)
    return proximal, distal


def _sample_model(
    calls: Sequence[annotate.ReadCall],
    fragments: dict[str, assemble.Fragment],
    sites: barcode.Sites,
    units: dict[str, barcode.ReadUnits],
) -> barcode.Model:
    """Scoring model with error rates learned from spanning reads.

    Spanning reads are grouped by chromosome, unit count, and A/B type. The
    groups do not depend on how well the barcodes agree, which would hide the
    errors that are to be measured.
    """
    rows = [barcode.matrix(read, sites) for read in units.values()]
    pool = (
        np.concatenate(rows)
        if rows
        else np.zeros((0, len(sites)), dtype=np.int8)
    )
    groups: dict[tuple, list[np.ndarray]] = {}
    for call in calls:
        fragment = fragments[call.read_id]
        if fragment.head and fragment.tail:
            key = (call.chromosome, call.ru_count, fragment.kind)
            groups.setdefault(key, []).append(fragment.matrix)
    error = barcode.estimate_error(
        [barcode.same_allele(group) for group in groups.values()]
    )
    return barcode.Model(pool, error)


def analyze(
    sequences: Mapping[str, str],
    aligned_chromosomes: Mapping[str, str],
    hits: Mapping[str, Sequence[blast.Hit]],
) -> SampleReads:
    """Calls every read and prepares the reads for allele calling.

    Args:
        sequences: read name to full sequence.
        aligned_chromosomes: contig of each mapped primary record.
        hits: BLAST hits per query file, from `blast.search`.
    """
    consensus = next(
        iter(seqlib.read_fasta(seqlib.data_path(blast.UNIT)).values())
    )
    lengths = {name: len(s) for name, s in sequences.items()}
    units = barcode.build_units(hits[blast.UNIT], lengths, consensus)
    calls = annotate.annotate(
        units,
        hits[blast.UNIT],
        hits[blast.MARKERS],
        hits[blast.P13E11],
        sequences,
        aligned_chromosomes,
    )
    # Fold-back reads cannot be placed on an array; sites come from the reads
    # that were assigned to chr4 or chr10.
    usable = [
        c
        for c in calls
        if not c.mixed_strand and c.chromosome in layout_module.CHROMOSOMES
    ]
    placed = {c.read_id: units[c.read_id] for c in usable}
    sites = barcode.find_sites(placed.values(), consensus)
    proximal, distal = _anchor_hits(hits)
    fragments = {}
    placed_units = {}
    misplaced: set[str] = set()
    for call in usable:
        read = units[call.read_id]
        fragment, kept = make_fragment(call, read, sites)
        for hits_of, near, attribute in (
            (proximal, True, "head"),
            (distal, False, "tail"),
        ):
            if getattr(fragment, attribute) and anchors_out_of_order(
                read, hits_of.get(call.read_id, []), near
            ):
                setattr(fragment, attribute, False)
                misplaced.add(call.read_id)
        fragments[call.read_id] = fragment
        placed_units[call.read_id] = kept
    model = _sample_model(usable, fragments, sites, placed)
    reliable = np.ones(len(sites), dtype=bool)
    if len(sites):
        limit = max(
            MAX_SITE_ERROR,
            MAX_RELATIVE_ERROR * float(np.median(model.error)),
        )
        reliable = model.error <= limit
    if not reliable.all():
        sites = sites.subset(reliable)
        for fragment in fragments.values():
            fragment.matrix = np.ascontiguousarray(fragment.matrix[:, reliable])
        model = _sample_model(usable, fragments, sites, placed)
    return SampleReads(
        calls,
        units,
        usable,
        sites,
        fragments,
        placed_units,
        misplaced,
        model,
    )
