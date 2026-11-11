# eris-validation: design

Agreed 2026-11-03 by Victoria Sjalim with Claude (design review, five rounds). Revised the same day: the truth is
always supplied from outside; this repo only validates.

## Goal

Measure how well eris reports, for each IS, which CDS lies upstream and downstream and how far away it is. Report
precision, recall and F1 for side and distance, and a row-level table of every wrong prediction. Run it on any genome
with Snakemake. (Effect, i.e. disrupted, was scored until 2026-11-11 and is now left out as not a priority; the truth
may still carry it.)

Terms are defined in [CONTEXT.md](../CONTEXT.md). Why the scorer lives here and not in eris:
[ADR 0001](adr/0001-independent-scorer.md).

## Inputs

Per assembly, one row of `config/samples.tsv`:

| column | meaning |
|---|---|
| `genome` | e.g. `EC958` |
| `assembly_type` | `complete`, `simulated` or `real` |
| `assembly` | FASTA or GFA given to eris |
| `reference` | complete reference genome FASTA, used to place eris genes |
| `truth` | relation truth TSV, curated outside this repo (schema below) |
| `replicons` | optional; limit the truth to these replicons (simulated EC958 GFA is chromosome-only) |

`config/config.yaml` holds everything that is not per assembly:

- `eris`: git URL and commit to install (default `adimascf/eris` main). The commit is written into every metrics row.
- `is_db`: ISfinder minimap2 index (`data/IS_n.mmi`). The workflow stops if it holds only a handful of sequences
  (the EC958 genome index `data/IS.mmi` was once used by mistake and made a GFA run hang).
- eris settings: collapse mode, `min_coverage 0.30`, scored run `hops 3`, reach run `hops 20`.
- Scoring thresholds (Dimas's values, so results compare with eris's own benchmark): placement overlap ≥ 50% of the
  shorter gene, near band ≤ 150 bp, nucleotide check ≥ 99% identity and ≥ 80% CDS coverage.

Isoform groups are fetched as data from `tests/isoform_groups.tsv` at the pinned eris commit.

### Truth schema

One row per IS–CDS relation, tab-separated, with a header. Same schema as `ec958_IS_cds_relations.tsv`.

| column | required | meaning |
|---|---|---|
| `is_id` | yes | IS name, optionally with a copy suffix (`IS26_6`) |
| `replicon` | yes | sequence id in the reference FASTA |
| `side` | yes | `upstream` or `downstream`, relative to the IS's strand; other values (e.g. `inside`) are ignored |
| `locus_tag` | yes | the CDS's identifier |
| `cds_start`, `cds_end` | yes | 0-based, half-open |
| `cds_strand` | yes | `+` or `-` |
| `distance_bp` | yes | gap between the IS and the CDS; negative when they overlap |
| `effect` | yes | `disrupted` or `none` are scored; anything else (`uncertain`, `upregulated`) is kept but unscored |
| `partial`, `pseudo` | no | `True`/`False`; a partial or pseudo CDS skips the coverage part of the nucleotide check |
| `inside_is` | no | set when the CDS lies inside a neighbouring IS; its effect is then unscored |

Other columns (`gene`, `product`, `note`, IS coordinates) are carried into `relations.tsv` untouched. The workflow
checks the schema first and stops on a missing column, an unknown replicon or a bad side. Two rows with the same
relation key (e.g. two copies of one IS with the same flank) are kept as one, the nearer.

## Steps

1. **Install eris** at the pinned commit (conda env).
2. **Check truth**: schema and replicons as above; limit to `replicons` if given.
3. **Run eris** twice per assembly (hops 3 scored, hops 20 reach), producing `_report.tsv` and `_cds.fna`.
4. **Score** (`erisval score`): reads `_report.tsv`, `_cds.fna`, the reference and the truth. Does not import eris.
   - Join report genes to `_cds.fna`: by `gene_id` (`<contig>_<start>_<end>`, 0-based) when records carry it, else by
     the coordinates in older pyfgs-style records (`<contig>_FGS_<n>`, 1-based). A reported gene missing from the FASTA
     stops the run.
   - Place each gene on the reference with minimap2; very short genes by exact search. A gene with several equally good
     placements matches if any of them does.
   - Match: same replicon, same strand, overlap ≥ 50% of the shorter of the two.
   - IS names: copy suffix stripped, then mapped to their isoform group. No family fallback.
5. **Combine** every assembly's metrics into `results/metrics.tsv`.

## Scoring rules

- **Relation key:** (isoform group, side, locus tag).
- **Recall:** truth relations found on their side. Pooled recall also accepts the wrong side.
- **Precision:** over first flanks only (nearest flank per side per graph branch), because the truth holds only the
  nearest CDS.
- **Status and reason:** every truth relation is TP or FN, every false first flank is FP. FN and FP carry one
  `reason`, checked in this order, so each points at one kind of bug:

  | status | reason | meaning |
  |---|---|---|
  | FN | `wrong_side` | its own IS reports the CDS only on the other side |
  | FN | `wrong_IS` | only a differently named IS reports the CDS, on the same side |
  | FN | `not_found` | not reported (`miss_cause` says budget or gene) |
  | FP | `unplaceable` | the gene has no placement on the reference |
  | FP | `wrong_side` | a truth CDS of its IS, but on the other side |
  | FP | `wrong_IS` | a truth CDS on the same side, but of a differently named IS |
  | FP | `not_in_truth` | in no truth relation |

- **Wrong side** is strict: FN on the true side and FP on the reported side.
- **Wrong IS** only says the gene belongs to another IS's relation; whether eris misnamed the IS or found a real
  neighbouring IS is judged per case, with `distance_error` as a hint (a small error means the same IS position).
- **Unplaceable genes** count as FP.
- **Distance:** on relations found on their side, truth and eris distances binned into overlapping (< 0, or eris overlap
  context), near (0–150 bp), far (> 150 bp). P/R/F1 per band, band agreement, and absolute error (median, max, share
  within 10 and 50 bp). eris reports overlaps as 0.
- **Effect:** not scored (since 2026-11-11). Truth and eris effects are carried into `relations.tsv` only.
- (Formerly: only disrupted was scored. eris `DISRUPTED` or `TRUNCATED` counted as disrupted. Upregulated was shown in the
  matrix, not scored.)
- **Miss cause:** a miss found by the hops-20 run is `budget`; otherwise `not_reported` when a gene pyfgs called
  (any record in `_cds.fna`) places on the truth CDS, else `not_called`.
- **Reports only:** no pass/fail. eris's own EC958 ratchet stays the only gate.

## Outputs

Per assembly in `<outdir>/<genome>/<assembly_type>/` (`outdir` in `config/config.yaml`):

| file | content |
|---|---|
| `side_confusion.tsv` | truth side (upstream, downstream, no_truth) × eris side (upstream, downstream, not_reported) |
| `distance_confusion.tsv` | truth band × eris band |
| `metrics.tsv` | quick summary, one row per class: `dimension, class, TP, FN, TP_calls, FP, precision, recall, f1` |
| `details.tsv` | further statistics, long format: `dimension, class, metric, value` |
| `relations.tsv` | one row per truth relation and per false flank: status, reason, IS, side, locus tag, eris gene, both distances, both effects (unscored), `distance_band_ok`, miss cause, nucleotide check |

`<outdir>/metrics.tsv` and `<outdir>/details.tsv` combine every assembly, with `genome, assembly_type` columns.
Recall is TP / (TP + FN) over truth relations; precision is TP_calls / (TP_calls + FP) over eris's first flanks (side)
or the found relations (distance).

## Layout

```
eris-validation/
  workflow/Snakefile, rules/*.smk, envs/*.yaml
  config/config.yaml, samples.tsv
  src/eris_validation/   truth.py (read + check), place.py, score.py, tables.py, cli.py
  tests/                 synthetic genome + unit tests; EC958 end-to-end (skipped without data/)
  docs/design.md, docs/adr/
  CONTEXT.md, README.md
```

Python with pandas, mappy and Biopython; Snakemake ≥ 8 with conda. Genome data stays out of git. Built test-first:
a tiny synthetic genome (about three IS on both strands, flanks on each side, one overlapping, one far) with a
hand-written eris report containing known errors, whose P/R/F1 is worked out by hand.

GitHub: `Victoriasja18/eris-validation`.

## Separate work

- `eris` fork: write an extra `_report.gff` holding only the CDS in `_report.tsv` (the full `_assembly.gff` stays).
  Its own PR, then merged into `eris-combined`. Not needed by the scorer.
- PR #29 (`synonyms.py`) is superseded by Dimas's isoform groups (#40–42): close it, or cut it down to
  `--merge-synonyms` if eris still needs that for copy counts.
