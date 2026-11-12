"""
A hand-worked scenario. Truth: IS26 has T1 upstream (near, none) and T2 downstream (overlapping, disrupted); ISEc9 has
T3 upstream (far, none) and T4 downstream (near, disrupted). eris:

- L1 (IS26) finds T1 (g1) and T2 (g2); g6 is a farther upstream flank, so not a first flank.
- L2 (ISEc9) reports T3's gene on the downstream side (g3), and an upstream gene that cannot be placed (g5). T4 is
  missed, but the hops-20 run finds it.
- L3 (IS1) reports a disrupted downstream gene that is in no truth relation (g7).
"""
import pytest

from eris_validation.eris_output import Flank
from eris_validation.place import Placement
from eris_validation.score import band, first_flanks, score
from eris_validation.tables import confusion, details, metrics
from eris_validation.truth import TruthRelation


def truth(is_id, group, side, tag, start, end, strand, distance, effect):
    return TruthRelation(is_id, group, side, tag, "chr", start, end, strand, False, False, distance, effect,
                         effect, {"gene": tag.lower()})


TRUTHS = [truth("IS26_1", "IS26", "upstream", "T1", 100, 400, "+", 50, "none"),
          truth("IS26_1", "IS26", "downstream", "T2", 2000, 2600, "-", -10, "disrupted"),
          truth("ISEc9_1", "ISEc9", "upstream", "T3", 5000, 5900, "+", 300, "none"),
          truth("ISEc9_1", "ISEc9", "downstream", "T4", 7000, 7600, "+", 20, "disrupted")]


def flank(locus, group, side, gene, distance, hops=0, effect="NONE", context=None, contig="ctg"):
    return Flank(locus, f"{group}_x", group, side, context or side.upper(), gene, contig, distance, hops, effect)


FLANKS = [flank("L1", "IS26", "upstream", "g1", 48),
          flank("L1", "IS26", "upstream", "g6", 900),
          flank("L1", "IS26", "downstream", "g2", 0, effect="DISRUPTED", context="OVERLAPPING_END"),
          flank("L2", "ISEc9", "upstream", "g5", 30),
          flank("L2", "ISEc9", "downstream", "g3", 310),
          flank("L3", "IS1", "downstream", "g7", 5, effect="DISRUPTED")]
REACH = FLANKS + [flank("L2", "ISEc9", "downstream", "g8", 25, hops=4)]
PLACES = {"g1": [Placement("chr", 102, 400, "+", 0.995)], "g2": [Placement("chr", 2000, 2590, "-", 0.95)],
          "g3": [Placement("chr", 5000, 5900, "+", 1.0)], "g5": [],
          "g6": [Placement("chr", 9000, 9500, "+", 1.0)], "g7": [Placement("chr", 12000, 12600, "+", 1.0)],
          "g8": [Placement("chr", 7000, 7600, "+", 1.0)]}


# Every gene pyfgs called on the assembly (every CDS in _assembly.gff), placed on the reference
CALLED = [p for gene in ("g1", "g2", "g3", "g6", "g7", "g8") for p in PLACES[gene]]


@pytest.fixture
def result():
    return score(TRUTHS, FLANKS, PLACES.__getitem__, reach_flanks=REACH, called=CALLED)


def by(rows, key):
    return {r[key]: r for r in rows}


def test_bands():
    assert [band(-1), band(0), band(150), band(151), band(0, overlapping=True)] == [
        "overlapping", "near", "near", "far", "overlapping"]


def test_first_flanks_are_the_nearest_per_side_per_locus():
    assert [f.gene_id for f in first_flanks(FLANKS)] == ["g1", "g2", "g5", "g3", "g7"]


def test_a_branch_beyond_the_locus_contig_gives_one_first_flank_per_contig():
    flanks = [flank("L", "A", "upstream", "a", 900, hops=1, contig="c1"),
              flank("L", "A", "upstream", "b", 100, hops=1, contig="c2"),
              flank("L", "A", "upstream", "c", 50, hops=2, contig="c3"),
              flank("L", "A", "upstream", "d", 950, hops=1, contig="c1")]
    assert [f.gene_id for f in first_flanks(flanks)] == ["b", "a"]


def status(row):
    return row["status"], row["reason"]


def test_every_truth_relation_and_false_flank_gets_a_status_and_reason(result):
    rows = by(result.rows, "eris_gene")
    assert [status(r) for r in result.rows if r["locus_tag"] == "T4"] == [("FN", "not_found")]
    assert [(g, status(rows[g])) for g in ("g1", "g2", "g3", "g5", "g7")] == [
        ("g1", ("TP", "")), ("g2", ("TP", "")), ("g3", ("FP", "wrong_side")), ("g5", ("FP", "unplaceable")),
        ("g7", ("FP", "not_in_truth"))]
    tags = by([r for r in result.rows if r["status"] in ("TP", "FN")], "locus_tag")
    assert status(tags["T3"]) == ("FN", "wrong_side") and tags["T3"]["eris_side"] == "downstream"


# Wrong IS: ISSen10 has T5 downstream; eris names that IS IS3411 and reports T5's gene (g9) as its nearest downstream
# flank. A third IS (ISEc8) reports T5's gene on the other side (g10), which is not a wrong-IS case.
T5 = truth("ISSen10_1", "ISSen10", "downstream", "T5", 20000, 20900, "+", 30, "none")
WRONG_IS_PLACES = {"g9": [Placement("chr", 20000, 20900, "+", 1.0)], "g10": [Placement("chr", 20000, 20900, "+", 1.0)]}


def test_a_truth_cds_reported_under_another_is_is_wrong_is_on_both_sides():
    flanks = [flank("L4", "IS3411", "downstream", "g9", 28), flank("L5", "ISEc8", "upstream", "g10", 28)]
    result = score([T5], flanks, WRONG_IS_PLACES.__getitem__)
    [truth_row] = [r for r in result.rows if r["locus_tag"] == "T5" and r["status"] == "FN"]
    assert status(truth_row) == ("FN", "wrong_IS") and truth_row["eris_target"] == "IS3411_x"
    flank_rows = by([r for r in result.rows if r["status"] == "FP"], "eris_gene")
    assert status(flank_rows["g9"]) == ("FP", "wrong_IS") and flank_rows["g9"]["is_id"] == "ISSen10_1"
    assert truth_row["distance_error"] == flank_rows["g9"]["distance_error"] == 2  # to judge: same IS position?
    assert status(flank_rows["g10"]) == ("FP", "not_in_truth")
    assert [f.reason for f in result.first] == ["wrong_IS", "not_in_truth"]


def test_wrong_side_is_checked_before_wrong_is():
    t5_up = truth("IS3411_1", "IS3411", "upstream", "T5", 20000, 20900, "+", 30, "none")
    flanks = [flank("L4", "IS3411", "downstream", "g9", 28)]
    result = score([T5, t5_up], flanks, WRONG_IS_PLACES.__getitem__)
    reasons = {(r["status"], r["is_id"]): r["reason"] for r in result.rows}
    assert reasons == {("FN", "ISSen10_1"): "wrong_IS", ("FN", "IS3411_1"): "wrong_side", ("FP", "IS3411_1"): "wrong_side"}


def test_matched_rows_flag_distance_and_sequence_but_not_effect(result):
    tags = by([r for r in result.rows if r["status"] == "TP"], "locus_tag")
    assert {k: tags["T1"][k] for k in ("distance_band_ok", "nucleotide_ok", "distance_error", "gene")} \
        == {"distance_band_ok": True, "nucleotide_ok": True, "distance_error": 2, "gene": "t1"}
    assert "effect_ok" not in tags["T1"]
    assert (tags["T2"]["truth_band"], tags["T2"]["eris_band"], tags["T2"]["nucleotide_ok"]) == (
        "overlapping", "overlapping", False)


def test_miss_cause_comes_from_the_reach_run_then_the_called_genes(result):
    tags = by([r for r in result.rows if r["status"] == "FN"], "locus_tag")
    # T4: found with --hops 20; T3: its gene (g3) is called, but never reported as an upstream flank of ISEc9
    assert (tags["T4"]["miss_cause"], tags["T3"]["miss_cause"]) == ("budget", "not_reported")


def test_a_miss_is_not_called_when_no_called_gene_places_on_the_truth_cds():
    t6 = truth("IS1_1", "IS1", "upstream", "T6", 30000, 30600, "+", 10, "none")
    t7 = truth("IS1_1", "IS1", "downstream", "T7", 31000, 31600, "-", 10, "none")
    called = [Placement("chr", 30000, 30600, "-", 1.0),  # on T6, but the other strand: a different gene
              Placement("chr", 31050, 31600, "-", 1.0)]  # T7's gene, called but in no report
    result = score([t6, t7], [], PLACES.__getitem__, reach_flanks=[], called=called)
    assert {r["locus_tag"]: r["miss_cause"] for r in result.rows} == {"T6": "not_called", "T7": "not_reported"}


def row(rows, dimension, class_):
    [r] = [r for r in rows if (r["dimension"], r["class"]) == (dimension, class_)]
    return r


def value(rows, dimension, class_, metric):
    [v] = [r["value"] for r in rows if (r["dimension"], r["class"], r["metric"]) == (dimension, class_, metric)]
    return v


def counts(r):
    return r["TP"], r["FN"], r["TP_calls"], r["FP"]


def test_side_metrics(result):
    m = metrics(result)
    assert counts(row(m, "side", "upstream")) == (1, 1, 1, 1)  # T1 found, T3 missed; g1 right, g5 unplaceable
    assert counts(row(m, "side", "downstream")) == (1, 1, 1, 2)  # T2 found, T4 missed; g2 right, g3 and g7 false
    both = row(m, "side", "both")
    assert counts(both) == (2, 2, 2, 3)
    assert both["precision"] == pytest.approx(0.4) and both["recall"] == 0.5 and both["f1"] == pytest.approx(4 / 9)
    assert row(m, "side", "downstream")["f1"] == pytest.approx(0.4)


def test_distance_metrics_and_no_effect_metrics(result):
    m = metrics(result)
    assert counts(row(m, "distance", "near")) == (1, 0, 1, 0) and row(m, "distance", "near")["f1"] == 1.0
    assert row(m, "distance", "far")["recall"] is None
    assert {r["dimension"] for r in m} == {"side", "distance"}


def test_details(result):
    d = details(result)
    assert value(d, "side", "both", "pooled_recall") == 0.75 and value(d, "side", "both", "wrong_side") == 1
    assert value(d, "side", "both", "wrong_IS") == 0 and value(d, "side", "both", "unplaceable") == 1
    assert [value(d, "side", "both", f"miss_{c}") for c in ("budget", "not_reported", "not_called")] == [1, 1, 0]
    assert value(d, "distance", "all", "agreement") == 1.0 and value(d, "distance", "all", "median_abs_error") == 1


def test_side_confusion(result):
    side = confusion(result, "side")
    assert side.loc["upstream"].to_dict() == {"upstream": 1, "downstream": 1, "not_reported": 0}
    assert side.loc["downstream"].to_dict() == {"upstream": 0, "downstream": 1, "not_reported": 1}
    assert side.loc["no_truth"].to_dict() == {"upstream": 1, "downstream": 1, "not_reported": 0}


def test_effect_is_not_scored(result):
    with pytest.raises(ValueError):
        confusion(result, "effect")
