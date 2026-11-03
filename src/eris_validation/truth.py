"""
Reads and checks a relation truth TSV (schema in docs/design.md). The truth is always curated outside this repo; a
file that breaks the schema stops the run rather than being scored wrong.
"""
import csv
from dataclasses import dataclass, field
from pathlib import Path

from eris_validation.isoforms import group_of

SIDES = ('upstream', 'downstream')
REQUIRED = ('is_id', 'replicon', 'side', 'locus_tag', 'cds_start', 'cds_end', 'cds_strand', 'distance_bp', 'effect')
SCORED_EFFECTS = ('disrupted', 'none')
RelationKey = tuple[str, str, str]  # (isoform group, side, locus tag)


class TruthError(ValueError):
    """The truth TSV breaks the schema."""


@dataclass(frozen=True)
class TruthRelation:
    is_id: str
    group: str
    side: str
    locus_tag: str
    replicon: str
    start: int  # 0-based, half-open
    end: int
    strand: str
    partial: bool
    pseudo: bool
    distance_bp: int  # negative when the CDS overlaps the IS
    effect: str | None  # 'disrupted' or 'none'; None when the effect is not scored
    raw_effect: str
    extra: dict[str, str] = field(default_factory=dict, compare=False, hash=False)  # passed through to relations.tsv

    @property
    def key(self) -> RelationKey:
        return self.group, self.side, self.locus_tag


def read_truth(path: Path, groups: dict[str, str], reference_ids: set[str] | None = None,
               replicons: set[str] | None = None) -> list[TruthRelation]:
    """
    The flank relations of the truth (rows whose side is not upstream or downstream are ignored), optionally only
    those on `replicons`. Two rows with the same relation key are kept as one, the nearer.
    """
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise TruthError(f"{path}: missing required column(s): {', '.join(missing)}")
        nearest: dict[RelationKey, TruthRelation] = {}
        for line, r in enumerate(reader, start=2):
            if r['side'] not in SIDES:
                if r['side'] != 'inside':
                    raise TruthError(f"{path} line {line}: side {r['side']!r} is not upstream, downstream or inside")
                continue
            t = _relation(r, groups, f"{path} line {line}")
            if reference_ids is not None and t.replicon not in reference_ids:
                raise TruthError(f"{path} line {line}: replicon {t.replicon!r} is not in the reference "
                                 f"({', '.join(sorted(reference_ids))})")
            if replicons is not None and t.replicon not in replicons:
                continue
            if t.key not in nearest or t.distance_bp < nearest[t.key].distance_bp:
                nearest[t.key] = t
    return list(nearest.values())


def _relation(r: dict[str, str], groups: dict[str, str], where: str) -> TruthRelation:
    try:
        start, end, distance = int(r['cds_start']), int(r['cds_end']), int(r['distance_bp'])
    except ValueError as e:
        raise TruthError(f"{where}: {e}") from None
    if end <= start:
        raise TruthError(f"{where}: cds_end {end} is not after cds_start {start}")
    if r['cds_strand'] not in ('+', '-'):
        raise TruthError(f"{where}: cds_strand {r['cds_strand']!r} is not + or -")
    effect = r['effect'] if r['effect'] in SCORED_EFFECTS and not r.get('inside_is') else None
    extra = {k: v for k, v in r.items() if k not in REQUIRED and k is not None}
    return TruthRelation(
        is_id=r['is_id'], group=group_of(r['is_id'], groups), side=r['side'], locus_tag=r['locus_tag'],
        replicon=r['replicon'], start=start, end=end, strand=r['cds_strand'], partial=r.get('partial') == 'True',
        pseudo=r.get('pseudo') == 'True', distance_bp=distance, effect=effect, raw_effect=r['effect'], extra=extra)
