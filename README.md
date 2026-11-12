# eris-validation

Validates how well [eris](https://github.com/Victoriasja18/eris) reports, for each insertion sequence (IS), which CDS lies
**upstream** and **downstream** and **how far** away. For every assembly it gives precision, recall and F1 for side
and distance, confusion matrices, and a row-level table of every wrong prediction. Effects (disrupted, truncated) are
not scored for now; they are only carried into `relations.tsv`.

The truth is always curated outside this repo. The scorer reads only eris's output files, so any eris version can be
validated. Design: [docs/design.md](docs/design.md). Terms: [CONTEXT.md](CONTEXT.md).

## Setup

With [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

or with conda:

```bash
conda env create -f environment.yaml && conda activate eris-validation
```

The workflow clones and installs eris itself (at the commit in `config/config.yaml`), so eris does not need to be
installed. Genome data is not in git; the example sample sheet expects it in `../data/`.

## Run

```bash
uv run snakemake --cores 4
```

Per assembly, eris runs twice: the scored run (`hops: 3`) and a reach run (`reach_hops: 20`) used only to tell budget
misses from gene misses. Results land in `<outdir>/<genome>/<assembly_type>/`, with every assembly combined in
`<outdir>/metrics.tsv` and `<outdir>/details.tsv`. Set `outdir` in `config/config.yaml`, or per run with
`--config outdir=results`.

To score an eris run you already have, without Snakemake:

```bash
uv run erisval score --report run_report.tsv --cds run_cds.fna --reference ref.fa --truth truth.tsv \
    --groups isoform_groups.tsv -o results/my_run
```

`--groups` is eris's `tests/isoform_groups.tsv`. Add `--reach-report`/`--reach-cds` for miss causes.

## Adding a genome

1. Prepare a truth TSV in the schema in [docs/design.md](docs/design.md#truth-schema) and check it:
   `uv run erisval check-truth --truth truth.tsv --reference ref.fa --groups isoform_groups.tsv`
2. Add a row per assembly to `config/samples.tsv`:

   | column | meaning |
   |---|---|
   | `genome` | e.g. `EC958` |
   | `assembly_type` | `complete`, `simulated` or `real` |
   | `assembly` | FASTA or GFA given to eris |
   | `reference` | complete reference genome FASTA |
   | `truth` | relation truth TSV |
   | `replicons` | optional, comma-separated: score only these replicons (e.g. a chromosome-only assembly) |

## Outputs

| file | content |
|---|---|
| `metrics.tsv` | the quick summary: `genome, assembly_type, dimension, class, TP, FN, TP_calls, FP, precision, recall, f1`, one row per class (side: upstream/downstream/both; distance: overlapping/near/far) |
| `details.tsv` | further statistics, `dimension, class, metric, value`: pooled recall, wrong-side, wrong-IS and miss counts, unplaceable genes, distance band agreement and errors |
| `relations.tsv` | one row per truth relation (`TP`, `FN`) and per false first flank (`FP`), with a `reason` for FN and FP (FN: `wrong_side`, `wrong_IS`, `not_found`; FP: `unplaceable`, `wrong_side`, `wrong_IS`, `not_in_truth`), both distances and effects (unscored), `distance_band_ok`, `distance_error`, `nucleotide_ok` and, for FN, `miss_cause` (`budget`, `not_reported`, `not_called`) |
| `side_confusion.tsv` | truth side × eris side |
| `distance_confusion.tsv` | truth distance band × eris band, for relations found on their side |

In `metrics.tsv`, recall = TP / (TP + FN) over the truth relations, and precision = TP_calls / (TP_calls + FP) over
eris's calls. For side the calls are eris's first flank per side per graph branch, while recall accepts
the true CDS at any flank, so TP_calls is usually lower than TP: when a wrong gene (often the IS's own transposase or
a short ORF) sits nearer the IS than the true CDS, the true CDS still counts as found, but eris's call is an FP. For
distance they are the same relations.

Every FN gets a miss cause, checked in order: `budget` (the `reach_hops` run reports it), `not_reported` (pyfgs called
a gene that is the truth CDS, but eris does not report it as a flank of that IS on that side), `not_called` (no gene
in `_cds.fna`, which holds every gene pyfgs called, is the truth CDS).
This workflow only reports; it never fails on low scores.

## Tests

```bash
uv run python -m pytest
```
