# eris-validation

Validates how well [eris](https://github.com/adimascf/eris) reports, for each insertion sequence (IS), which CDS lies
**upstream** and **downstream**, **how far** away, and whether the insertion **disrupts** it. For every assembly it
gives precision, recall and F1 for side, distance and effect, confusion matrices, and a row-level table of every wrong
prediction.

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
misses from gene misses. Results land in `results/<genome>/<assembly_type>/`, all of them combined in
`results/metrics.tsv`.

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
| `metrics.tsv` | `genome, assembly_type, eris_commit, hops, dimension, class, metric, value` |
| `relations.tsv` | one row per truth relation (`TP`, `FN`, `FN_wrong_side`) and per false first flank (`FP`, `FP_wrong_side`, `FP_unplaceable`), with both distances and effects, `distance_band_ok`, `effect_ok`, `distance_error`, `nucleotide_ok` and `miss_cause` |
| `side_confusion.tsv` | truth side × eris side |
| `distance_confusion.tsv` | truth distance band × eris band, for relations found on their side |
| `effect_confusion.tsv` | truth effect × eris effect |

Metrics: side P/R/F1 per side and both (precision over eris's first flank per side per graph branch), pooled recall;
distance P/R/F1 per band, band agreement and absolute error; disrupted P/R/F1 (eris `TRUNCATED` counts as disrupted).
This workflow only reports; it never fails on low scores.

## Tests

```bash
uv run python -m pytest
```
