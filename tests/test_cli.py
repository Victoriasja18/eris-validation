import random
from pathlib import Path

import pandas as pd
import pytest

from eris_validation.cli import main

DATA = Path(__file__).parents[2] / "data"


def test_score_writes_every_table_with_labels(tmp_path):
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(20_000))
    (tmp_path / "ref.fa").write_text(f">chr\n{seq}\n")
    (tmp_path / "groups.tsv").write_text("name\tgroup\nIS15DI\tIS26\n")
    (tmp_path / "truth.tsv").write_text(
        "is_id\treplicon\tside\tlocus_tag\tcds_start\tcds_end\tcds_strand\tdistance_bp\teffect\n"
        "IS26_1\tchr\tupstream\tT1\t1000\t1900\t+\t40\tnone\n"
        "IS26_1\tchr\tdownstream\tT2\t5000\t5600\t+\t60\tdisrupted\n")
    (tmp_path / "report.tsv").write_text(
        "locus_id\ttarget\tgene_id\tcontext\tdist_bp\ttopo_hops\teffect\n"
        "L1\tIS15DI_IS6_unknown\tchr_1000_1900\tUPSTREAM\t42\t0\tNONE\n")
    (tmp_path / "cds.fna").write_text(f">chr_1000_1900 strand=+\n{seq[1000:1900]}\n")
    out = tmp_path / "out"

    main(["score", "--report", str(tmp_path / "report.tsv"), "--cds", str(tmp_path / "cds.fna"),
          "--reference", str(tmp_path / "ref.fa"), "--truth", str(tmp_path / "truth.tsv"),
          "--groups", str(tmp_path / "groups.tsv"), "-o", str(out), "--label", "genome=toy"])

    m = pd.read_csv(out / "metrics.tsv", sep="\t")
    assert list(m.columns) == ["genome", "dimension", "class", "metric", "value"]
    both = m.set_index(["dimension", "class", "metric"])["value"]
    assert both["side", "both", "recall"] == 0.5 and both["side", "both", "precision"] == 1.0
    statuses = pd.read_csv(out / "relations.tsv", sep="\t").set_index("locus_tag")["status"].to_dict()
    assert statuses == {"T1": "TP", "T2": "FN"}
    assert all((out / f"{d}_confusion.tsv").exists() for d in ("side", "distance", "effect"))


def test_a_genome_index_is_not_taken_for_an_is_database(tmp_path):
    (tmp_path / "ref.fa").write_text(">chr\n" + "ACGT" * 5000 + "\n")
    with pytest.raises(SystemExit, match="genome index"):
        main(["check-is-db", "--db", str(tmp_path / "ref.fa")])


@pytest.mark.slow
@pytest.mark.skipif(not (DATA / "ec958_IS_cds_relations.tsv").exists(), reason="EC958 data not in ../data")
def test_ec958_truth_passes_the_schema_check(tmp_path, capsys):
    (tmp_path / "groups.tsv").write_text("name\tgroup\n")
    main(["check-truth", "--truth", str(DATA / "ec958_IS_cds_relations.tsv"),
          "--reference", str(DATA / "EC958_complete_genome.fasta"), "--groups", str(tmp_path / "groups.tsv")])
    assert "117 relations (59 upstream, 58 downstream)" in capsys.readouterr().out
