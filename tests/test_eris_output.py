import pytest

from eris_validation.eris_output import ErisOutputError, disrupting, read_cds, read_flanks

REPORT = (
    "locus_id\ttarget\ttarget_identity\tcontig_path\tgene_id\tcontext\tdist_bp\ttopo_hops\torientation\teffect\n"
    "L1\tIS15DI_IS6_unknown,ISEc9_IS5_unknown\t0.99\tchr+\tchr_10_90\tOVERLAPPING_START\t0\t0\tsame strand\tDISRUPTED\n"
    "L1\tIS15DI_IS6_unknown,ISEc9_IS5_unknown\t0.99\tchr+\tchr_100_200\tINSIDE\t0\t0\tsame strand\tNONE\n"
    "L1\tIS15DI_IS6_unknown,ISEc9_IS5_unknown\t0.99\tchr+\t5_0_60\tDOWNSTREAM\t40\t1\tsame strand\tNONE\n"
    "L2\tIS1A_IS1_unknown\t1.0\tctg+\tNone\tNONE\t-1\t0\t-\tNONE\n"
)


def test_flanks_take_the_first_target_s_group_and_skip_inside_and_empty_rows(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPORT)
    flanks = read_flanks(path, {"IS15DI": "IS26"})
    assert [(f.locus_id, f.group, f.side, f.gene_id, f.contig, f.distance_bp, f.hops, f.overlapping)
            for f in flanks] == [("L1", "IS26", "upstream", "chr_10_90", "chr", 0, 0, True),
                                 ("L1", "IS26", "downstream", "5_0_60", "5", 40, 1, False)]


def test_truncated_or_disrupted_counts_as_disrupting():
    assert disrupting("TRUNCATED") and disrupting("DISRUPTED|UPREGULATED") and not disrupting("NONE")


@pytest.mark.parametrize("fasta", [
    ">chr_10_90 strand=+\nATGAAA\n",  # named by the report's gene_id
    ">chr_FGS_1 chr:11-90 strand=+\nATGAAA\n",  # pyfgs name, 1-based start
    ">IS26_UPSTREAM_chr_FGS_1 chr:10-90 strand=+\nATGAAA\n",  # renamed, 0-based start
])
def test_cds_found_under_every_naming_style(tmp_path, fasta):
    path = tmp_path / "c.fna"
    path.write_text(fasta + ">chr_FGS_2 chr:301-400 strand=-\nATGCCC\n")
    assert read_cds(path, ["chr_10_90"]) == {"chr_10_90": "ATGAAA"}


def test_a_gene_missing_from_the_fasta_stops(tmp_path):
    path = tmp_path / "c.fna"
    path.write_text(">chr_FGS_2 chr:301-400 strand=-\nATGCCC\n")
    with pytest.raises(ErisOutputError, match="chr_10_90"):
        read_cds(path, ["chr_10_90"])
