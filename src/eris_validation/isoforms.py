"""
IS names to isoform groups, read from eris's `tests/isoform_groups.tsv` (built from ISfinder's synonym and isoform
annotation). A name not in the table is its own group; names are never grouped by how they look.
"""
import csv
from pathlib import Path


def read_groups(path: Path) -> dict[str, str]:
    """{IS name: group} for every name in the table."""
    with open(path) as handle:
        rows = csv.DictReader((line for line in handle if not line.startswith('#')), delimiter='\t')
        return {r['name']: r['group'] for r in rows}


def group_of(name: str, groups: dict[str, str]) -> str:
    """
    The isoform group of a truth IS id or an eris target, e.g. 'IS15DI_2' or 'IS15DI_IS6_unknown_stitched' -> 'IS26'.
    Everything after the first underscore (copy suffix, family, stitching tag) is dropped first.
    """
    name = name.split('_')[0]
    return groups.get(name, name)
