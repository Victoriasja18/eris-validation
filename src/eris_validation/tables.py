"""
Turns a scored result into the output tables: confusion matrices for side, distance and effect, and a long-format
metrics table (`dimension, class, metric, value`). A metric with nothing to compute it from is left empty (None).
"""
from pathlib import Path
from statistics import median

import pandas as pd

from eris_validation.eris_output import disrupting
from eris_validation.score import Result
from eris_validation.truth import SIDES

TRUTH_STATUSES = ('TP', 'FN', 'FN_wrong_side')
FALSE_FLANKS = ('FP', 'FP_unplaceable')  # flanks with no truth relation on either side
BANDS = ('overlapping', 'near', 'far')


def ratio(a: int, b: int) -> float | None:
    return a / b if b else None


def f1(p: float | None, r: float | None) -> float | None:
    if p is None or r is None:
        return None
    return 2 * p * r / (p + r) if p + r else 0.0


def truth_rows(result: Result) -> list[dict]:
    return [r for r in result.rows if r['status'] in TRUTH_STATUSES]


def confusion(result: Result, dimension: str) -> pd.DataFrame:
    """Truth class (rows) × eris class (columns) for 'side', 'distance' or 'effect'."""
    truths = truth_rows(result)
    false = [r for r in result.rows if r['status'] in FALSE_FLANKS]
    if dimension == 'side':
        pairs = [(r['side'], r['eris_side'] or 'not_reported') for r in truths]
        pairs += [('no_truth', r['eris_side']) for r in false]
        index, columns = [*SIDES, 'no_truth'], [*SIDES, 'not_reported']
    elif dimension == 'distance':
        pairs = [(r['truth_band'], r['eris_band']) for r in truths if r['status'] == 'TP']
        index = columns = list(BANDS)
    elif dimension == 'effect':
        pairs = [(_truth_effect(r), r['eris_effect'] if r['status'] == 'TP' else 'not_found') for r in truths]
        pairs += [('no_truth', r['eris_effect']) for r in false]
        index = ['disrupted', 'none', 'unscored', 'no_truth']
        columns = sorted({e for _, e in pairs} - {'not_found'}) + ['not_found']
    else:
        raise ValueError(f"unknown dimension {dimension!r}")
    table = pd.DataFrame(0, index=pd.Index(index, name='truth'), columns=pd.Index(columns, name='eris'))
    for truth, eris in pairs:
        table.loc[truth, eris] += 1
    return table


def _truth_effect(row: dict) -> str:
    return row['truth_effect'] if _scored(row) else 'unscored'


def _scored(row: dict) -> bool:
    return row['truth_effect'] in ('disrupted', 'none') and not row.get('inside_is')


def metrics(result: Result) -> list[dict]:
    rows = []

    def add(dimension, class_, metric, value):
        rows.append({'dimension': dimension, 'class': class_, 'metric': metric, 'value': value})

    # Side: recall over truth relations, precision over first flanks
    truths = truth_rows(result)
    for side in (*SIDES, 'both'):
        on_side = [r for r in truths if side in ('both', r['side'])]
        firsts = [f for f in result.first if side in ('both', f.flank.side)]
        found = sum(r['status'] == 'TP' for r in on_side)
        true_firsts = sum(f.status == 'TP' for f in firsts)
        p, r = ratio(true_firsts, len(firsts)), ratio(found, len(on_side))
        for metric, value in (('precision', p), ('recall', r), ('f1', f1(p, r)), ('truth', len(on_side)),
                              ('found', found), ('wrong_side', sum(r['status'] == 'FN_wrong_side' for r in on_side)),
                              ('first_flanks', len(firsts)), ('true_first_flanks', true_firsts),
                              ('unplaceable', sum(f.status == 'FP_unplaceable' for f in firsts))):
            add('side', side, metric, value)
    add('side', 'both', 'pooled_recall', ratio(sum(r['status'] != 'FN' for r in truths), len(truths)))
    for cause in ('gene', 'budget'):
        add('side', 'both', f'miss_{cause}', sum(r['miss_cause'] == cause for r in truths))

    # Distance: bands of the relations found on their side
    table = confusion(result, 'distance')
    for b in BANDS:
        p, r = ratio(table.loc[b, b], table[b].sum()), ratio(table.loc[b, b], table.loc[b].sum())
        add('distance', b, 'precision', p)
        add('distance', b, 'recall', r)
        add('distance', b, 'f1', f1(p, r))
    errors = [r['distance_error'] for r in truths if r['status'] == 'TP']
    add('distance', 'all', 'n', len(errors))
    add('distance', 'all', 'agreement', ratio(sum(table.loc[b, b] for b in BANDS), len(errors)))
    add('distance', 'all', 'median_abs_error', median(errors) if errors else None)
    add('distance', 'all', 'max_abs_error', max(errors, default=None))
    for within in (10, 50):
        add('distance', 'all', f'within_{within}bp', ratio(sum(e <= within for e in errors), len(errors)))

    # Effect: disrupted only. Precision over first flanks called DISRUPTED or TRUNCATED, leaving out those that
    # match a truth relation whose effect is not scored; recall over the truth's disrupted relations.
    called = [f for f in result.first if disrupting(f.flank.effect) and all(t.effect is not None for t in f.matched)]
    called_true = sum(any(t.effect == 'disrupted' for t in f.matched) for f in called)
    disrupted = [r for r in truths if _scored(r) and r['truth_effect'] == 'disrupted']
    disrupted_found = sum(r['status'] == 'TP' and disrupting(r['eris_effect']) for r in disrupted)
    p, r = ratio(called_true, len(called)), ratio(disrupted_found, len(disrupted))
    for metric, value in (('precision', p), ('recall', r), ('f1', f1(p, r)), ('truth', len(disrupted)),
                          ('called', len(called))):
        add('effect', 'disrupted', metric, value)
    return rows


def write_tables(result: Result, outdir: Path):
    """Writes the three confusion matrices, metrics.tsv and relations.tsv into `outdir`."""
    outdir.mkdir(parents=True, exist_ok=True)
    for dimension in ('side', 'distance', 'effect'):
        confusion(result, dimension).to_csv(outdir / f"{dimension}_confusion.tsv", sep='\t')
    pd.DataFrame(metrics(result)).to_csv(outdir / "metrics.tsv", sep='\t', index=False)
    pd.DataFrame(result.rows).to_csv(outdir / "relations.tsv", sep='\t', index=False)
