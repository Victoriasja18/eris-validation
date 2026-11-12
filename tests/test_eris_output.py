import pytest

from eris_validation.eris_output import ErisOutputError, read_all_cds, read_called_genes, read_cds, read_flanks

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


def test_all_called_genes_are_read_whatever_their_names(tmp_path):
    path = tmp_path / "cds.fna"
    path.write_text(">chr_FGS_1 chr:11-90 strand=+\nATGAAA\n>chr_FGS_2 chr:101-200 strand=-\nATGCCC\nTAA\n")
    assert read_all_cds(path) == ["ATGAAA", "ATGCCCTAA"]


GFF = ("##gff-version 3\n"
       "ctgA\tpyfgs_v0.0.3\tCDS\t3\t8\t1.3\t+\t0\tID=ctgA_FGS_1\n"
       "ctgB\tpyfgs_v0.0.3\tCDS\t2\t7\t1.1\t-\t0\tID=ctgB_FGS_1\n")


@pytest.mark.parametrize("assembly", [
    ">ctgA\nCCATGAAAGG\n>ctgB length=10\nGTTACATTTT\n",  # FASTA
    "H\tVN:Z:1.0\nS\tctgA\tCCATGAAAGG\tdp:i:30\nS\tctgB\tGTTACATTTT\nL\tctgA\t+\tctgB\t+\t0M\n",  # GFA
])
def test_called_genes_are_cut_from_the_assembly_by_the_assembly_gff(tmp_path, assembly):
    (tmp_path / "a.gff").write_text(GFF)
    (tmp_path / "asm").write_text(assembly)
    # GFF coordinates are 1-based and inclusive; a '-' gene is read on the reverse strand
    assert read_called_genes(tmp_path / "a.gff", tmp_path / "asm") == ["ATGAAA", "ATGTAA"]


def test_a_called_gene_on_a_contig_missing_from_the_assembly_stops(tmp_path):
    (tmp_path / "a.gff").write_text(GFF)
    (tmp_path / "asm.fa").write_text(">ctgA\nCCATGAAAGG\n")
    with pytest.raises(ErisOutputError, match="ctgB"):
        read_called_genes(tmp_path / "a.gff", tmp_path / "asm.fa")


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
