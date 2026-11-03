import random

from Bio.Seq import reverse_complement

from eris_validation.place import GenePlacer, Placement


def genome(tmp_path, length=20_000, seed=1):
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    path = tmp_path / "ref.fa"
    path.write_text(f">chr\n{seq}\n")
    return path, seq


def test_a_gene_is_placed_on_the_strand_it_reads_on(tmp_path):
    path, seq = genome(tmp_path)
    placer = GenePlacer(path)
    assert placer.place(seq[1000:1900]) == [Placement("chr", 1000, 1900, "+", 1.0)]
    assert placer.place(reverse_complement(seq[5000:5600])) == [Placement("chr", 5000, 5600, "-", 1.0)]


def test_a_gene_too_short_to_align_is_found_by_exact_search(tmp_path):
    path, seq = genome(tmp_path)
    assert GenePlacer(path).place(reverse_complement(seq[300:330])) == [Placement("chr", 300, 330, "-", 1.0)]


def test_a_repeated_gene_gets_every_equally_good_place(tmp_path):
    path, seq = genome(tmp_path)
    gene = seq[1000:1900]
    path.write_text(f">chr\n{seq[:10_000]}{gene}{seq[10_000:]}\n")
    assert {p.start for p in GenePlacer(path).place(gene)} == {1000, 10_000}


def test_a_gene_from_elsewhere_has_no_place(tmp_path):
    path, _ = genome(tmp_path)
    _, other = genome(tmp_path / "..", seed=2)
    assert GenePlacer(path).place(other[:900]) == []
