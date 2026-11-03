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

**Wrong side**:
A truth CDS that eris reports on the opposite side; an FN on its true side and an FP on the reported one.

**Unplaceable gene**:
An eris gene that cannot be aligned to the reference, so it can match no truth relation.

**Distance band**:
Overlapping, near (0–150 bp) or far; the classes distance is scored in.

**Miss cause**:
Why a truth relation was not found: **budget** (found with a larger hops budget) or **gene** (wrong or absent gene call).

## Relationships

- A **truth relation** is found when an eris flank of the same **isoform group**, on the same **side**, places onto its CDS.
- Precision counts **first flanks**; recall counts **truth relations**. One collapsed repeat can find a relation from several branches.

## Flagged ambiguities

- "Synonym" was used for both ISfinder synonyms and isoforms; resolved: both are one **isoform group**.
- "Scorer" meant Dimas's benchmark inside eris; here it is this repo's own, which reads eris's output files only.
