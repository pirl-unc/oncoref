# Annotation-scoped CTA gene identities

Canonical CTA sets are reference-space IDs, not all source IDs that can appear in
an annotation. For example, OncoRef's shipped Ensembl alias table maps PRAME
`ENSG00000275013` to canonical `ENSG00000185686`. Both have translations in human
Ensembl 93 and 112. Comparing annotated source IDs directly with
`cta_unfiltered_gene_ids()` leaves the alternate PRAME translations in a negative
self-reference and can mask all 501 distinct canonical PRAME 9-mers.

Ensembl's [haplotype and patch documentation](https://www.ensembl.org/info/genome/genebuild/haplotypes_patches.html)
explains these alternative assembly paths. OncoRef verifies **explicit ID
mappings** from its shipped table; a shared symbol or identical protein does not
create a new gene alias or justify removing a genuine non-CTA match.

## Use the same resolver for admission and exclusion

```python
from pyensembl import EnsemblRelease
from oncoref import (
    cta_annotation_gene_identities,
    cta_annotation_gene_ids,
    cta_gene_ids,
    resolve_gene_identity,
)

genome = EnsemblRelease(93, species="human")  # explicitly chosen, locally installed
identity = resolve_gene_identity("ENSG00000275013.2", genome=genome)
assert identity.verified
assert identity.canonical_gene_id == "ENSG00000185686"

# Admission uses canonical default membership after verified ID resolution.
admitted = identity.verified and identity.canonical_gene_id in cta_gene_ids()
# Or obtain the equivalent source-ID set for this annotation:
admission_source_ids = cta_annotation_gene_ids(genome)

# Negative references exclude the entire candidate universe, independently of
# patient expression, target-family policy or the default admitted panel.
exclude_source_ids = cta_annotation_gene_ids(genome, unfiltered=True)
assert "ENSG00000275013" in exclude_source_ids

# Persist these decisions alongside the consumer's own annotation-file hashes.
audit = [record.as_dict() for record in cta_annotation_gene_identities(genome)]
```

Keep each source's original gene, transcript, protein, species and occurrence IDs.
Use the resolved canonical gene ID only for joining OncoRef facts. Do not relabel
source occurrences or discard all matches to an identical sequence: a sequence
can also come from a genuine non-CTA source, which stays in the self-reference.

## Contract version 1

`GeneIdentity` is immutable. `input_gene_id` preserves the caller's original ID;
`source_gene_id` is the unversioned annotation lookup ID. Records include the
annotation species, assembly and release; source symbol, contig, coordinates and
strand; explicit alias rows and method; canonical reference releases; OncoRef
version; identity contract version; and SHA-256 hashes of sorted logical rows of
the shipped alias and canonical gene-space tables. These hashes are independent
of CSV compression and ordering. They identify the mapping reference, not the
consumer's complete GTF/FASTA contents, which consumers must hash separately.

Only `canonical` and `alias` outcomes are verified. `absent`, `unmapped`,
`ambiguous`, `unverified` and `annotation_conflict` have no resolved canonical ID.
`candidate_gene_ids` and `alias_rows` retain proposed targets and evidence even
when a decision is rejected. Rejected sources remain in negative references and
do not qualify for default CTA admission.

- A canonical stable ID must exist on the reference's recorded chromosome in the
  supplied annotation. Its symbol may have changed between releases.
- An `ensembl_alt_contig` row requires its recorded symbol to agree with the
  source annotation and unique primary-locus canonical reference. Its source
  contig must be non-primary. If the canonical target also exists in the supplied
  annotation, its symbol and contig must agree too.
- An `ensembl_id_history_name` row requires the same recorded symbol and primary
  chromosome, with the successor absent from that annotation. Coexisting IDs
  are retained conservatively because they need further evidence.
- Competing targets, conflicting evidence, unknown alias methods and aliases
  without symbols are unresolved. Unknown same-symbol loci are never inferred.
- An explicit human GRCh37/GRCh38 Ensembl annotation is required. A missing gene
  yields `absent`; other query failures propagate.

The underlying alias-generation algorithm is in
[`scripts/generate_ensembl_id_aliases.py`](https://github.com/pirl-unc/oncoref/blob/main/scripts/generate_ensembl_id_aliases.py).
The current canonical table is Ensembl 115. Alias rows were generated across
multiple installed historical releases and preserved curated rows; they are not
a claim that every alias is present in every release. Runtime annotation checks
establish which supplied sources can use a mapping. The legacy
`resolve_ensembl_id()` and `canonical_gene_id(source_version=...)` APIs remain
release-agnostic table-harmonization helpers; they do not implement this
annotation-verification contract.

## Peptide counts and reference versions

`cta_specific_9mer_counts()` now uses the same verified identities to exclude
annotated CTA source IDs from its background. It prefers a verified canonical
source's longest protein; if that source is absent, it uses the longest protein
of a verified annotated copy. All unverified and genuine non-CTA sources stay in
the background. The cache fingerprint includes identity decisions, reference
hashes and the contract version, so pre-fix counts cannot be reused.

This changes gene-ID handling, not CTA curation. The canonical default and
candidate sets and the downloadable expression/source bundle versions remain
unchanged. Downstream consumers must adopt these APIs to correct their own
annotation-specific references; upgrading OncoRef alone does not rewrite a
consumer's saved dataset or previously serialized admission evidence.
