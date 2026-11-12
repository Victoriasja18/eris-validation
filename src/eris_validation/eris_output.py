"""
Reads the eris outputs the validation needs: the flanks in `_report.tsv`, the coding sequences of the reported genes
in `_cds.fna`, and every gene pyfgs called on the assembly. Older eris versions wrote every called gene to `_cds.fna`;
newer ones keep only the reported genes there, and every called gene is cut from the assembly by `_assembly.gff`.
Only the files are read, never eris itself, so any eris version's output can be scored.
"""
import csv
import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import reverse_complement

from eris_validation.isoforms import group_of

# A gene crossing an end of the IS is reported as a passenger, but it is a flank on the side it sticks out of
CONTEXT_SIDES = {'UPSTREAM': 'upstream', 'DOWNSTREAM': 'downstream',
                 'OVERLAPPING_START': 'upstream', 'OVERLAPPING_END': 'downstream'}
COORDINATES = re.compile(r'(\S+):(\d+)-(\d+)')  # '<contig>:<start>-<end>' in a pyfgs-style description


class ErisOutputError(ValueError):
    """The eris outputs disagree with each other."""


@dataclass(frozen=True)
class Flank:
    """One flank of one locus in the report."""
    locus_id: str
    target: str
    group: str
    side: str
    context: str
    gene_id: str
    contig: str
    distance_bp: int
    hops: int
    effect: str

    @property
    def overlapping(self) -> bool:
        return self.context.startswith('OVERLAPPING')


def split_gene_id(gene_id: str) -> tuple[str, int, int]:
    """'<contig>_<start>_<end>' (0-based, half-open) -> (contig, start, end)."""
    contig, start, end = gene_id.rsplit('_', 2)
    return contig, int(start), int(end)


def read_flanks(path: Path, groups: dict[str, str]) -> list[Flank]:
    """Every flank in the report, in report order. A locus is named by its first target."""
    flanks = []
    with open(path) as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            side = CONTEXT_SIDES.get(r['context'])
            if side is None or r['gene_id'] == 'None':
                continue
            target = r['target'].split(',')[0]
            flanks.append(Flank(
                locus_id=r['locus_id'], target=target, group=group_of(target, groups), side=side,
                context=r['context'], gene_id=r['gene_id'], contig=split_gene_id(r['gene_id'])[0],
                distance_bp=int(r['dist_bp']), hops=int(r['topo_hops']), effect=r['effect']))
    return flanks


def read_all_cds(path: Path) -> list[str]:
    """The coding sequence of every gene pyfgs called on the assembly."""
    return [str(record.seq) for record in SeqIO.parse(path, 'fasta')]


def read_called_genes(gff: Path, assembly: Path) -> list[str]:
    """
    The coding sequence of every gene pyfgs called on the assembly: each CDS in eris's `_assembly.gff` (1-based,
    inclusive), cut from the assembly (FASTA or GFA) and read on its own strand.
    """
    contigs = _read_assembly(assembly)
    sequences = []
    with open(gff) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 9 or fields[2] != 'CDS':
                continue
            contig, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            if contig not in contigs:
                raise ErisOutputError(f"called gene on contig {contig} in {gff} is not in {assembly}")
            seq = contigs[contig][start - 1:end]
            sequences.append(reverse_complement(seq) if strand == '-' else seq)
    return sequences


def _read_assembly(path: Path) -> dict[str, str]:
    """{contig: sequence} from a GFA (its S lines) or a FASTA."""
    with open(path) as handle:
        first = next((line for line in handle if line.strip()), '')
    if first.startswith('>'):
        return {record.id: str(record.seq) for record in SeqIO.parse(path, 'fasta')}
    contigs = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('S\t'):
                fields = line.rstrip('\n').split('\t')
                contigs[fields[1]] = fields[2]
    return contigs


def read_cds(path: Path, gene_ids) -> dict[str, str]:
    """
    {gene_id: coding sequence} for the report's genes. A record is found by its id when eris names it by gene_id,
    otherwise by the '<contig>:<start>-<end>' in its description, whose start may be 0- or 1-based depending on the
    eris version. The end is matched exactly, so a sequence whose length disagrees with its coordinates (pyfgs
    frameshift corrections) is still found.
    """
    by_id, by_end = {}, {}
    for record in SeqIO.parse(path, 'fasta'):
        seq = str(record.seq)
        by_id[record.id] = seq
        if m := COORDINATES.search(record.description.removeprefix(record.id)):
            by_end.setdefault((m[1], int(m[3])), []).append((int(m[2]), seq))
    found = {}
    for gene_id in gene_ids:
        if gene_id in by_id:
            found[gene_id] = by_id[gene_id]
            continue
        contig, start, end = split_gene_id(gene_id)
        hits = [seq for s, seq in by_end.get((contig, end), []) if s in (start, start + 1)]
        if len(hits) != 1:
            problem = "is missing from" if not hits else "matches several records in"
            raise ErisOutputError(f"reported gene {gene_id} {problem} {path}")
        found[gene_id] = hits[0]
    return found
