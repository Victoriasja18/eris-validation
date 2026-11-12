"""
erisval: score eris's IS–CDS relations against a curated truth.

    erisval score --report R --cds C --reference REF --truth T --groups G -o OUTDIR [--reach-report ... --reach-cds ...]
    erisval check-truth --truth T --reference REF --groups G
    erisval check-is-db --db IS_n.mmi
    erisval combine -o results/metrics.tsv run1/metrics.tsv run2/metrics.tsv ...
"""
import argparse
import sys
from functools import cache
from pathlib import Path

import pandas as pd

from eris_validation.eris_output import ErisOutputError, read_all_cds, read_called_genes, read_cds, read_flanks
from eris_validation.isoforms import read_groups
from eris_validation.place import GenePlacer
from eris_validation.score import Thresholds, score
from eris_validation.tables import details, metrics, write_tables
from eris_validation.truth import TruthError, read_truth

MIN_IS_DB_SEQUENCES = 100  # an ISfinder index holds thousands; a genome index a handful


def main(argv=None):
    parser = argparse.ArgumentParser(prog='erisval', description=__doc__.strip().splitlines()[0])
    commands = parser.add_subparsers(dest='command', required=True)

    s = commands.add_parser('score', help='score one eris run against the truth')
    s.add_argument('--report', type=Path, required=True, help='eris _report.tsv')
    s.add_argument('--cds', type=Path, required=True, help="eris _cds.fna of the same run")
    s.add_argument('--assembly', type=Path, help='the assembly eris ran on (FASTA or GFA), with --assembly-gff')
    s.add_argument('--assembly-gff', type=Path,
                   help="eris _assembly.gff: every gene pyfgs called, to tell not_called from not_reported. Without "
                        "it, every gene in --cds is taken as called, which holds only for eris versions whose "
                        "_cds.fna is not limited to the reported genes")
    s.add_argument('--reach-report', type=Path, help='_report.tsv of a larger --hops run, to tell budget misses')
    s.add_argument('--reach-cds', type=Path, help='_cds.fna of that run')
    s.add_argument('--reference', type=Path, required=True, help='reference genome FASTA')
    s.add_argument('--truth', type=Path, required=True, help='relation truth TSV')
    s.add_argument('--groups', type=Path, required=True, help="isoform groups TSV (eris tests/isoform_groups.tsv)")
    s.add_argument('--replicons', help='comma-separated; limit the truth to these replicons')
    s.add_argument('--label', action='append', default=[], metavar='KEY=VALUE',
                   help='column added to metrics.tsv and details.tsv, e.g. genome=EC958 (repeatable)')
    s.add_argument('-o', '--outdir', type=Path, required=True)
    for name, default in vars(Thresholds()).items():
        s.add_argument(f"--{name.replace('_', '-')}", type=type(default), default=default)

    c = commands.add_parser('check-truth', help='check a truth TSV against the schema and reference')
    c.add_argument('--truth', type=Path, required=True)
    c.add_argument('--reference', type=Path, required=True)
    c.add_argument('--groups', type=Path, required=True)

    d = commands.add_parser('check-is-db', help='check that an eris target index holds an IS database')
    d.add_argument('--db', type=Path, required=True)

    m = commands.add_parser('combine', help='concatenate metrics.tsv (or details.tsv) files')
    m.add_argument('tables', type=Path, nargs='+')
    m.add_argument('-o', '--out', type=Path, required=True)

    args = parser.parse_args(argv)
    try:
        {'score': run_score, 'check-truth': run_check_truth, 'check-is-db': run_check_is_db,
         'combine': run_combine}[args.command](args)
    except (TruthError, ErisOutputError) as e:
        sys.exit(f"erisval: {e}")


def run_score(args):
    if (args.reach_report is None) != (args.reach_cds is None):
        sys.exit("erisval: --reach-report and --reach-cds go together")
    if (args.assembly is None) != (args.assembly_gff is None):
        sys.exit("erisval: --assembly and --assembly-gff go together")
    labels = dict(label.split('=', 1) for label in args.label)
    groups = read_groups(args.groups)
    placer = GenePlacer(args.reference)
    replicons = set(args.replicons.split(',')) if args.replicons else None
    truths = read_truth(args.truth, groups, reference_ids=placer.replicons, replicons=replicons)

    flanks = read_flanks(args.report, groups)
    sequences = read_cds(args.cds, {f.gene_id for f in flanks})
    reach_flanks = None
    if args.reach_report:
        reach_flanks = read_flanks(args.reach_report, groups)
        sequences |= read_cds(args.reach_cds, {f.gene_id for f in reach_flanks} - sequences.keys())

    @cache
    def place(gene_id):
        return placer.place(sequences[gene_id])

    th = Thresholds(**{name: getattr(args, name) for name in vars(Thresholds())})
    # Every called gene, to tell not_called from not_reported
    called_genes = (read_called_genes(args.assembly_gff, args.assembly) if args.assembly_gff
                    else read_all_cds(args.cds))
    called = [p for seq in called_genes for p in placer.place(seq)]
    result = score(truths, flanks, place, reach_flanks, th, called)
    write_tables(result, args.outdir, labels)
    print(summary(metrics(result), details(result)), file=sys.stderr)


def summary(rows: list[dict], extra: list[dict]) -> str:
    value = {(r['dimension'], r['class'], r['metric']): r['value'] for r in extra}
    fmt = lambda v: '-' if v is None else f"{v:.2f}"
    lines = [f"{'':<12}{'P':>6}{'R':>6}{'F1':>6}"]
    for r in rows:
        if r['dimension'] != 'distance':
            lines.append(f"{r['class']:<12}" + "".join(f"{fmt(r[m]):>6}" for m in ('precision', 'recall', 'f1')))
    lines.append(f"distance band agreement {fmt(value['distance', 'all', 'agreement'])}, "
                 f"pooled recall {fmt(value['side', 'both', 'pooled_recall'])}")
    return "\n".join(lines)


def run_check_truth(args):
    placer_ids = {line[1:].split()[0] for line in open(args.reference) if line.startswith('>')}
    truths = read_truth(args.truth, read_groups(args.groups), reference_ids=placer_ids)
    sides = {side: sum(t.side == side for t in truths) for side in ('upstream', 'downstream')}
    scored = sum(t.effect is not None for t in truths)
    print(f"{args.truth}: {len(truths)} relations ({sides['upstream']} upstream, {sides['downstream']} downstream), "
          f"{scored} with a scored effect")


def run_check_is_db(args):
    import mappy
    aligner = mappy.Aligner(str(args.db))
    names = aligner.seq_names if aligner else []
    if len(names) < MIN_IS_DB_SEQUENCES:
        sys.exit(f"erisval: {args.db} holds {len(names)} sequence(s) ({', '.join(names[:3])}); an IS database "
                 f"holds thousands. Is it a genome index?")
    print(f"{args.db}: {len(names)} sequences")


def run_combine(args):
    tables = [pd.read_csv(path, sep='\t') for path in args.tables]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
