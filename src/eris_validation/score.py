"""
Matches eris's flanks to the truth and gives every truth relation and every false first flank a status:

- truth relations: TP (found on its side), FN_wrong_side (found only on the other side), FN (not found);
- eris first flanks that match no truth relation on their side: FP_wrong_side (they match one on the other side),
  FP_unplaceable (the gene has no place on the reference), FP.

An eris gene matches a truth CDS when one of its placements on the reference lies on the same replicon and strand
and overlaps it by at least `min_overlap` of the shorter of the two. Only flanks of loci in the same isoform group
count. Recall is over truth relations and looks at every flank; precision is over first flanks only, because the
truth holds only the nearest CDS on each side.
"""
from collections.abc import Callable
from dataclasses import dataclass, field

from eris_validation.eris_output import Flank, disrupting
from eris_validation.place import Placement
from eris_validation.truth import SIDES, TruthRelation

Place = Callable[[str], list[Placement]]  # eris gene_id -> its placements on the reference


@dataclass(frozen=True)
class Thresholds:
    min_overlap: float = 0.5  # of the shorter of the eris gene and the truth CDS
    near: int = 150  # bp; pyfgs start codons and eris's IS ends differ from RefSeq by a few to tens of bp
    min_identity: float = 0.99  # nucleotide check, over the aligned part
    min_cds_coverage: float = 0.8  # nucleotide check, of the truth CDS; waived when it is partial or pseudo


@dataclass
class FirstFlank:
    flank: Flank
    status: str  # 'TP' or one of the FP statuses
    matched: list[TruthRelation]  # truth relations on its side that it matches


@dataclass
class Result:
    rows: list[dict] = field(default_factory=list)  # relations.tsv
    first: list[FirstFlank] = field(default_factory=list)
    truth_total: dict[str, int] = field(default_factory=dict)


def band(distance: int, overlapping: bool = False, near: int = Thresholds.near) -> str:
    if overlapping or distance < 0:
        return 'overlapping'
    return 'near' if distance <= near else 'far'


def first_flanks(flanks: list[Flank]) -> list[Flank]:
    """
    The nearest flank of each branch on each side of each locus. A flank on the IS's own contig (hop 0) is the only
    branch on its side; otherwise each contig at the fewest hops that carries a gene starts a branch, and flanks
    further along are left out, since eris does not report which branch a deeper contig continues.
    """
    loci: dict[str, list[Flank]] = {}
    for f in flanks:
        loci.setdefault(f.locus_id, []).append(f)
    firsts = []
    for locus_flanks in loci.values():
        for side in SIDES:
            on_side = sorted((f for f in locus_flanks if f.side == side), key=lambda f: (f.hops, f.distance_bp))
            if not on_side:
                continue
            nearest_hop, contigs = on_side[0].hops, set()
            for f in on_side:
                if f.hops == nearest_hop and f.contig not in contigs:
                    contigs.add(f.contig)
                    firsts.append(f)
                if nearest_hop == 0:
                    break
    return firsts


def matches(p: Placement, t: TruthRelation, th: Thresholds) -> bool:
    if p.replicon != t.replicon or p.strand != t.strand:
        return False
    overlap = min(p.end, t.end) - max(p.start, t.start)
    return overlap >= th.min_overlap * min(p.end - p.start, t.end - t.start)


def sequence_agrees(p: Placement, t: TruthRelation, th: Thresholds) -> bool:
    """Whether a matched gene has the truth CDS's sequence; boundaries of partial or pseudo CDS are not checked."""
    covered = (min(p.end, t.end) - max(p.start, t.start)) / (t.end - t.start)
    return p.identity >= th.min_identity and (covered >= th.min_cds_coverage or t.partial or t.pseudo)


def score(truths: list[TruthRelation], flanks: list[Flank], place: Place, reach_flanks: list[Flank] | None = None,
          th: Thresholds = Thresholds()) -> Result:
    """
    Scores eris's flanks against the truth. `reach_flanks` come from the same assembly run with a larger hops
    budget: a miss they find is put down to the budget rather than to a wrong or absent gene.
    """
    def matching(t: TruthRelation, candidates: list[Flank]) -> list[Flank]:
        return [f for f in candidates if f.group == t.group and any(matches(p, t, th) for p in place(f.gene_id))]

    result = Result(truth_total={side: sum(t.side == side for t in truths) for side in SIDES})
    for t in truths:
        found = matching(t, flanks)
        on_side = sorted((f for f in found if f.side == t.side), key=lambda f: (f.hops, f.distance_bp))
        if on_side:
            result.rows.append(_found_row(t, on_side[0], place, th))
            continue
        cause = ''
        if reach_flanks is not None:
            cause = 'budget' if any(f.side == t.side for f in matching(t, reach_flanks)) else 'gene'
        other = sorted(found, key=lambda f: (f.hops, f.distance_bp))
        row = _row('FN_wrong_side' if other else 'FN', truth=t, flank=other[0] if other else None, th=th)
        row['miss_cause'] = cause
        result.rows.append(row)

    for f in first_flanks(flanks):
        placements = place(f.gene_id)
        hits = [t for t in truths if t.group == f.group and any(matches(p, t, th) for p in placements)]
        matched = [t for t in hits if t.side == f.side]
        if matched:
            result.first.append(FirstFlank(f, 'TP', matched))
            continue
        status = 'FP_unplaceable' if not placements else 'FP_wrong_side' if hits else 'FP'
        result.first.append(FirstFlank(f, status, []))
        row = _row(status, truth=hits[0] if hits else None, flank=f, th=th)
        row['side'] = f.side
        if placements:
            p = placements[0]
            row['eris_placement'] = f"{p.replicon}:{p.start}-{p.end}({p.strand})"
        result.rows.append(row)
    return result


def _found_row(t: TruthRelation, f: Flank, place: Place, th: Thresholds) -> dict:
    row = _row('TP', truth=t, flank=f, th=th)
    placements = [p for p in place(f.gene_id) if matches(p, t, th)]
    p = placements[0]
    row.update(
        eris_placement=f"{p.replicon}:{p.start}-{p.end}({p.strand})",
        distance_band_ok=row['truth_band'] == row['eris_band'],
        effect_ok=None if t.effect is None else disrupting(f.effect) == (t.effect == 'disrupted'),
        distance_error=abs(f.distance_bp - max(0, t.distance_bp)),  # eris reports an overlap as distance 0
        nucleotide_ok=any(sequence_agrees(p, t, th) for p in placements))
    return row


def _row(status: str, truth: TruthRelation | None, flank: Flank | None, th: Thresholds) -> dict:
    t, f = truth, flank
    row = {
        'status': status,
        'group': t.group if t else f.group,
        'side': t.side if t else f.side,
        'is_id': t.is_id if t else '',
        'locus_tag': t.locus_tag if t else '',
        'truth_location': f"{t.replicon}:{t.start}-{t.end}({t.strand})" if t else '',
        'truth_distance': t.distance_bp if t else None,
        'truth_band': band(t.distance_bp, near=th.near) if t else '',
        'truth_effect': t.raw_effect if t else '',
        'eris_locus': f.locus_id if f else '',
        'eris_target': f.target if f else '',
        'eris_gene': f.gene_id if f else '',
        'eris_side': f.side if f else '',
        'eris_context': f.context if f else '',
        'eris_hops': f.hops if f else None,
        'eris_distance': f.distance_bp if f else None,
        'eris_band': band(f.distance_bp, f.overlapping, th.near) if f else '',
        'eris_effect': f.effect if f else '',
        'eris_placement': '',
        'distance_band_ok': None,
        'effect_ok': None,
        'distance_error': None,
        'nucleotide_ok': None,
        'miss_cause': '',
    }
    return row | (t.extra if t else {})
