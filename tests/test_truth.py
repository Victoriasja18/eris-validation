import pytest

from eris_validation.isoforms import group_of, read_groups
from eris_validation.truth import TruthError, read_truth

HEADER = "is_id\treplicon\tside\tlocus_tag\tcds_start\tcds_end\tcds_strand\tdistance_bp\teffect\tpartial\tinside_is\tgene\n"
GROUPS = {"IS15DI": "IS26", "IS26": "IS26"}


def write(tmp_path, *rows, header=HEADER):
    path = tmp_path / "truth.tsv"
    path.write_text(header + "".join("\t".join(r) + "\n" for r in rows))
    return path


def row(is_id="IS26_6", replicon="chr", side="upstream", tag="T1", start="100", end="400", strand="+",
        distance="72", effect="none", partial="False", inside_is="", gene="blaX"):
    return [is_id, replicon, side, tag, start, end, strand, distance, effect, partial, inside_is, gene]


def test_isoform_group_strips_the_copy_suffix_and_target_tail():
    assert group_of("IS15DI_2", GROUPS) == "IS26"
    assert group_of("IS15DI_IS6_unknown_stitched", GROUPS) == "IS26"
    assert group_of("ISEc9_1", GROUPS) == "ISEc9"


def test_isoform_groups_file_skips_comments(tmp_path):
    path = tmp_path / "groups.tsv"
    path.write_text("# from ISfinder\nname\tgroup\nIS1A\tIS1\nIS1R\tIS1\n")
    assert read_groups(path) == {"IS1A": "IS1", "IS1R": "IS1"}


def test_reads_flank_rows_and_ignores_inside_rows(tmp_path):
    path = write(tmp_path, row(), row(side="inside", tag="T2"))
    [t] = read_truth(path, GROUPS)
    assert (t.group, t.side, t.locus_tag, t.start, t.end, t.strand) == ("IS26", "upstream", "T1", 100, 400, "+")
    assert t.distance_bp == 72 and t.effect == "none" and t.extra["gene"] == "blaX"


def test_effect_is_unscored_inside_a_neighbouring_is_or_when_not_disrupted_or_none(tmp_path):
    path = write(tmp_path, row(tag="T1", effect="disrupted", inside_is="ISEc9_1"), row(tag="T2", effect="uncertain"),
                 row(tag="T3", effect="disrupted"))
    assert [t.effect for t in read_truth(path, GROUPS)] == [None, None, "disrupted"]


def test_two_copies_with_the_same_flank_keep_the_nearer(tmp_path):
    path = write(tmp_path, row(is_id="IS26_2", distance="451"), row(is_id="IS15DI_3", distance="256"))
    [t] = read_truth(path, GROUPS)
    assert t.distance_bp == 256 and t.is_id == "IS15DI_3"


def test_replicons_limit_the_truth(tmp_path):
    path = write(tmp_path, row(replicon="chr"), row(replicon="plasmid", tag="T2"))
    assert [t.locus_tag for t in read_truth(path, GROUPS, replicons={"chr"})] == ["T1"]


def test_missing_required_column_stops(tmp_path):
    path = write(tmp_path, row()[:8], header="\t".join(HEADER.split("\t")[:8]) + "\n")
    with pytest.raises(TruthError, match="effect"):
        read_truth(path, GROUPS)


def test_unknown_replicon_stops(tmp_path):
    with pytest.raises(TruthError, match="NZ_chr"):
        read_truth(write(tmp_path, row(replicon="NZ_chr")), GROUPS, reference_ids={"chr"})


@pytest.mark.parametrize("field, value", [("side", "left"), ("strand", "."), ("start", "x"), ("end", "50")])
def test_bad_values_stop_with_the_line_number(tmp_path, field, value):
    with pytest.raises(TruthError, match="line 2"):
        read_truth(write(tmp_path, row(**{field: value})), GROUPS)
