"""
Places eris genes on the reference genome. eris genes are pyfgs predictions on the assembly, not RefSeq CDS, so each
one is aligned to the reference before it can be compared with the truth.
"""
from dataclasses import dataclass
from pathlib import Path

import mappy
from Bio import SeqIO
from Bio.Seq import reverse_complement


@dataclass(frozen=True)
class Placement:
    """Where an eris gene lies on the reference: `strand` is the strand its coding sequence reads on."""
    replicon: str
    start: int  # 0-based, half-open
    end: int
    strand: str
    identity: float  # over the aligned part


class GenePlacer:
    """
    Aligns genes to the reference with minimap2. A gene that aligns equally well to several places (e.g. a repeated
    gene) gets all of them; one too short to align is looked up by exact search.
    """

    def __init__(self, reference: Path):
        self._aligner = mappy.Aligner(str(reference), preset='asm5')
        if not self._aligner:
            raise ValueError(f"could not index {reference}")
        self._replicons = {r.id: str(r.seq).upper() for r in SeqIO.parse(reference, 'fasta')}

    @property
    def replicons(self) -> set[str]:
        return set(self._replicons)

    def place(self, sequence: str) -> list[Placement]:
        seq = sequence.upper()
        hits = list(self._aligner.map(seq))
        if hits:
            best = max(h.mlen for h in hits)
            return [Placement(h.ctg, h.r_st, h.r_en, '+' if h.strand == 1 else '-', h.mlen / h.blen)
                    for h in hits if h.mlen == best]
        return self._search(seq)

    def _search(self, seq: str) -> list[Placement]:
        placements = []
        for name, replicon in self._replicons.items():
            for strand, query in (('+', seq), ('-', reverse_complement(seq))):
                start = replicon.find(query)
                while start != -1:
                    placements.append(Placement(name, start, start + len(seq), strand, 1.0))
                    start = replicon.find(query, start + 1)
        return placements
