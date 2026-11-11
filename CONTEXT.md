# eris-validation

Scores eris's IS–CDS relations against a relation truth set curated outside this repo. Terms owned by eris itself (flank,
relation, context, effect, disrupted, truncated, upregulated) keep eris's meaning.

## Language

**Truth relation**:
One CDS on one side of one IS in the reference, with its distance and effect. Always supplied from outside, never built here.
_Avoid_: ground truth row, expected flank

**Relation key**:
What makes two relations the same: isoform group, side and locus tag. IS copies are not told apart.

**Isoform group**:
ISfinder names that denote the same element (synonyms and isoforms), scored as one name. Never inferred from look-alike names.
_Avoid_: synonym, family

**Side**:
Upstream or downstream, read in the direction of the IS's own strand.

**First flank**:
eris's nearest reported flank on one side of one graph branch; the only flanks precision is computed over.

**Reason**:
Why an FN or FP is wrong, one per row, checked in order: FN `wrong_side`, `wrong_IS`, `not_found`; FP `unplaceable`,
`wrong_side`, `wrong_IS`, `not_in_truth`. Status stays TP / FN / FP.
_Avoid_: FN_wrong_side, FP_wrong_side, FP_unplaceable (old status names)

**Wrong side**:
A truth CDS that its own IS reports on the opposite side; an FN on its true side and an FP on the reported one.

**Wrong IS**:
A truth CDS that eris reports on the right side but under a differently named IS (another isoform group). Whether eris
misnamed the IS or found a real neighbouring IS is judged per case.

**Unplaceable gene**:
An eris gene that cannot be aligned to the reference, so it can match no truth relation.

**Distance band**:
Overlapping, near (0–150 bp) or far; the classes distance is scored in.

**Miss cause**:
Why a truth relation was not found, checked in order: **budget** (found with a larger hops budget), **not_reported**
(pyfgs called the gene, but eris does not report it as a flank of that IS on that side), **not_called** (pyfgs called
no gene that is the truth CDS).
_Avoid_: gene (the old catch-all for not_reported and not_called)

## Relationships

- A **truth relation** is found when an eris flank of the same **isoform group**, on the same **side**, places onto its CDS.
- Precision counts **first flanks**; recall counts **truth relations**, found at any flank. A truth relation is still found when a wrong gene sits nearer the IS than its CDS, while that wrong gene is the first flank and counts against precision.

## Flagged ambiguities

- "Synonym" was used for both ISfinder synonyms and isoforms; resolved: both are one **isoform group**.
- "Scorer" meant Dimas's benchmark inside eris; here it is this repo's own, which reads eris's output files only.
