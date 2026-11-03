# Score from eris's output files, outside eris

eris already has a relation scorer (`tests/relation_scoring.py`), but it reads eris's in-memory objects inside eris's
test suite, where it serves as a regression gate on EC958. We wrote our own scorer here that reads only `_report.tsv`
and `_cds.fna`, so any eris version or old output can be scored on any genome without importing eris. We keep Dimas's
matching rules and thresholds so results stay comparable, and take isoform groups as data from the pinned eris commit.

## Considered options

- Move eris's scorer into `src/eris/benchmark/` and call it from here: one implementation, but it ties scoring to eris
  internals and to changes in Dimas's repo (it changed three times on 2026-11-03).
- Run `tests/relation_scoring.py` from a pinned eris checkout: no copy, but depends on a test file.

## Consequences

Two scorers can drift. If eris's rules change, the thresholds in `config.yaml` need updating by hand.
