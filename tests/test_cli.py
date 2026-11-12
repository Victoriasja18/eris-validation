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
    assert list(m.columns) == ["genome", "dimension", "class", "TP", "FN", "TP_calls", "FP", "precision", "recall",
                               "f1"]
    both = m.set_index(["dimension", "class"]).loc[("side", "both")]
    assert (both["TP"], both["FN"], both["recall"], both["precision"]) == (1, 1, 0.5, 1.0)
    assert list(pd.read_csv(out / "details.tsv", sep="\t").columns) == ["genome", "dimension", "class", "metric",
                                                                         "value"]
    relations = pd.read_csv(out / "relations.tsv", sep="\t").set_index("locus_tag")
    assert relations["status"].to_dict() == {"T1": "TP", "T2": "FN"}
    assert relations.loc["T2", "miss_cause"] == "not_called"  # cds.fna holds no gene on T2
    assert all((out / f"{d}_confusion.tsv").exists() for d in ("side", "distance"))
    assert not (out / "effect_confusion.tsv").exists()


def test_called_genes_come_from_the_assembly_gff_when_given(tmp_path):
    # eris's _cds.fna holds only reported genes; its _assembly.gff lists every called gene, here the one on T2
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
    (tmp_path / "assembly.gff").write_text(
        "##gff-version 3\n"
        "chr\tpyfgs\tCDS\t1001\t1900\t1\t+\t0\tID=chr_FGS_1\n"
        "chr\tpyfgs\tCDS\t5001\t5600\t1\t+\t0\tID=chr_FGS_2\n")
    out = tmp_path / "out"

    main(["score", "--report", str(tmp_path / "report.tsv"), "--cds", str(tmp_path / "cds.fna"),
          "--assembly", str(tmp_path / "ref.fa"), "--assembly-gff", str(tmp_path / "assembly.gff"),
          "--reference", str(tmp_path / "ref.fa"), "--truth", str(tmp_path / "truth.tsv"),
          "--groups", str(tmp_path / "groups.tsv"), "-o", str(out)])

    relations = pd.read_csv(out / "relations.tsv", sep="\t").set_index("locus_tag")
    assert relations.loc["T2", "miss_cause"] == "not_reported"


def test_assembly_and_assembly_gff_go_together(tmp_path):
    with pytest.raises(SystemExit, match="go together"):
        main(["score", "--report", "r", "--cds", "c", "--assembly-gff", "g", "--reference", "ref",
              "--truth", "t", "--groups", "g", "-o", str(tmp_path)])


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
