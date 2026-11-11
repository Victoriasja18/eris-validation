"""
Turns a scored result into the output tables: confusion matrices for side, distance and effect, and a long-format
summary table (counts and precision, recall, F1 per class), plus a long table of further statistics. A metric with
nothing to compute it from is left empty (None).
"""
from pathlib import Path
from statistics import median

import pandas as pd

from eris_validation.eris_output import disrupting
from eris_validation.score import Result
from eris_validation.truth import SIDES

TRUTH_STATUSES = ('TP', 'FN')
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
    # False flanks with no truth relation of their IS on either side; a wrong-side flank is already in its truth's row
    false = [r for r in result.rows if r['status'] == 'FP' and r['reason'] != 'wrong_side']
    if dimension == 'side':
        pairs = [(r['side'], r['eris_side'] if r['reason'] in ('', 'wrong_side') else 'not_reported') for r in truths]
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


def _row(dimension: str, class_: str, tp: int, fn: int, tp_calls: int, fp: int) -> dict:
    p, r = ratio(tp_calls, tp_calls + fp), ratio(tp, tp + fn)
    return {'dimension': dimension, 'class': class_, 'TP': tp, 'FN': fn, 'TP_calls': tp_calls, 'FP': fp,
            'precision': p, 'recall': r, 'f1': f1(p, r)}


def metrics(result: Result) -> list[dict]:
    """
    The quick summary: one row per scored class with its counts and precision, recall and F1. Recall is TP / (TP + FN)
    over the truth relations; precision is TP_calls / (TP_calls + FP) over eris's calls. For side and effect the calls
    are eris's first flanks while recall accepts the true CDS at any flank, so TP_calls is usually lower than TP: a
    wrong gene nearer the IS than the true CDS is an FP call, yet the true CDS still counts as found. For distance
    they are the same relations.
    """
    rows = []
    truths = truth_rows(result)
    for side in (*SIDES, 'both'):
        on_side = [r for r in truths if side in ('both', r['side'])]
        firsts = [f for f in result.first if side in ('both', f.flank.side)]
        tp, tp_calls = sum(r['status'] == 'TP' for r in on_side), sum(f.status == 'TP' for f in firsts)
        rows.append(_row('side', side, tp, len(on_side) - tp, tp_calls, len(firsts) - tp_calls))

    # Distance: bands of the relations found on their side
    table = confusion(result, 'distance')
    for b in BANDS:
        tp = int(table.loc[b, b])
        rows.append(_row('distance', b, tp, int(table.loc[b].sum()) - tp, tp, int(table[b].sum()) - tp))

    # Effect: disrupted only. Calls are first flanks called DISRUPTED or TRUNCATED, leaving out those that match a
    # truth relation whose effect is not scored; the truth is its disrupted relations.
    called = [f for f in result.first if disrupting(f.flank.effect) and all(t.effect is not None for t in f.matched)]
    called_true = sum(any(t.effect == 'disrupted' for t in f.matched) for f in called)
    disrupted = [r for r in truths if _scored(r) and r['truth_effect'] == 'disrupted']
    found = sum(r['status'] == 'TP' and disrupting(r['eris_effect']) for r in disrupted)
    rows.append(_row('effect', 'disrupted', found, len(disrupted) - found, called_true, len(called) - called_true))
    return rows


def details(result: Result) -> list[dict]:
    """Statistics beyond the summary, long format: `dimension, class, metric, value`."""
    rows = []

    def add(dimension, class_, metric, value):
        rows.append({'dimension': dimension, 'class': class_, 'metric': metric, 'value': value})

    truths = truth_rows(result)
    add('side', 'both', 'pooled_recall', ratio(sum(r['reason'] in ('', 'wrong_side') for r in truths), len(truths)))
    for reason in ('wrong_side', 'wrong_IS'):
        add('side', 'both', reason, sum(r['reason'] == reason for r in truths))
    for cause in ('gene', 'budget'):
        add('side', 'both', f'miss_{cause}', sum(r['miss_cause'] == cause for r in truths))
    add('side', 'both', 'unplaceable', sum(f.reason == 'unplaceable' for f in result.first))

    table = confusion(result, 'distance')
    errors = [r['distance_error'] for r in truths if r['status'] == 'TP']
    add('distance', 'all', 'n', len(errors))
    add('distance', 'all', 'agreement', ratio(sum(table.loc[b, b] for b in BANDS), len(errors)))
    add('distance', 'all', 'median_abs_error', median(errors) if errors else None)
    add('distance', 'all', 'max_abs_error', max(errors, default=None))
    for within in (10, 50):
        add('distance', 'all', f'within_{within}bp', ratio(sum(e <= within for e in errors), len(errors)))
    return rows


def write_tables(result: Result, outdir: Path, labels: dict[str, str] | None = None):
    """
    Writes the three confusion matrices, metrics.tsv, details.tsv and relations.tsv into `outdir`. `labels` (e.g.
    genome, assembly_type) become the first columns of metrics.tsv and details.tsv.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    for dimension in ('side', 'distance', 'effect'):
        confusion(result, dimension).to_csv(outdir / f"{dimension}_confusion.tsv", sep='\t')
    for name, rows in (('metrics', metrics(result)), ('details', details(result))):
        table = pd.DataFrame(rows).round(3)
        for i, (key, value) in enumerate((labels or {}).items()):
            table.insert(i, key, value)
        table.to_csv(outdir / f"{name}.tsv", sep='\t', index=False)
    pd.DataFrame(result.rows).to_csv(outdir / "relations.tsv", sep='\t', index=False)
