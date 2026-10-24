# CTA, Unsafe and Neither

The [full CTA atlas](audits/cta-atlas-20261024/index.md) adds cancer occurrence,
source-linked biology, clinical outcomes, reconsideration tables and four PDF
reports. It flags 37 retained-gene literature conflicts, including urgent neural
protein review for ABCC11 and MAGEC2, without silently changing the frozen sets.

The [current report](audits/cta-neither-review-20261024/index.md) supplies three disjoint
sets across 20,097 canonical human protein-coding genes and the pinned Ensembl
115 protein inventory.

| Category | Target | Forbidden background |
|---|---|---|
| CTA | Yes, under the research screen | No |
| Unsafe | No | Yes |
| Neither | No | No |

CTA requires the previous curated candidate screen plus affirmative reproductive
evidence: reproductive RNA **≥5 nTPM**, reliable positive reproductive protein,
or cited reproductive plus tumor/recognition evidence. Reliable reproductive IHC
can satisfy admission below 5 nTPM. None of these bypass normal-expression gates.

The current screen leaves **399 CTA genes**. A hard ceiling of **>30 TPM or nTPM**
outside explicitly mapped Reproductive/Disposable contexts moves 253 prior Neither
genes and 37 prior CTAs into Unsafe. It applies regardless of IHC confidence or
reproductive enrichment. MAGEA3 and ROPN1 are benchmark losses; 19 of the prior 21
benchmark genes remain. CT55 also remains retained. The dated 23 October report
preserves the previous 436-gene decisions.

Unsafe requires positive normal IHC, above-ceiling normal RNA, a reviewed safety
hold, or cited normal protein production outside permitted compartments. ABCB5
and IL5 also now have the latter classification. A selectivity failure alone
does not make a gene Unsafe. Original bulk category/IHC ceilings are unchanged.

Neither preserves Low, Variable, Cellular, Antibody, Unassessed and other
unresolved subgroups. CTAG2 remains Variable. FOXR2 has an Antibody hold;
HORMAD1, TDRD1, BRDT and DAZL have Cellular review holds. ACRBP/OY-TES-1 and
TEX101 now enter Unsafe / High RNA. No remaining Neither gene has an audited
non-permitted/unsited TPM or nTPM observation above 30. Corneal CP10K remains a
separate rule. [Candidate follow-up](audits/cta-neither-review-20261024/follow-up.md)
records questions and evidence needed for reconsideration.
The [normal-cell guide](normal-cell-evidence.md) documents the separate native-unit
RNA review triggers. Missing data never become zero.

These are research classifications: CTA is not proven clinical safety, Unsafe
is a screening exclusion rather than proven toxicity, and Neither is not proof
of normal-protein absence. Original DICE columns cannot be paired across genes;
only CTAG2 has the previously reviewed vital-organ individual-donor RNA evidence.

## Python

```python
from oncoref import target_sets, normal_rna_review

evidence = target_sets.results()
cta = target_sets.gene_ids("CTA")
forbidden = target_sets.gene_ids("Unsafe")
neither = target_sets.gene_ids("Neither")
names = target_sets.gene_names("CTA")  # known symbols; IDs include unnamed loci
reviews = normal_rna_review.followup_reviews()
sequence_loci = normal_rna_review.sequence_gene_reviews()
all_forbidden_genes = normal_rna_review.forbidden_gene_ids()
```

All three ID sets are disjoint. Earlier `cta`, `cta_screen` and normal-expression
exports retain their dated definitions. Use `target_sets` and the current report
for this policy.

`target_sets` retains its canonical 20,097-gene scope. The sequence inventory
also contains 1,721 exact reference loci outside that table, including
immune-receptor segments. Their RNA audit is exposed through
`sequence_gene_reviews()`. `forbidden_gene_ids()` includes Unsafe decisions from
both scopes and matches the gene decisions used for the protein exports.

## Exact proteins and k-mers

The report rebuilds complete exact-sequence sets from the checksum-pinned
Ensembl 115 inventory, preserving every source transcript and encoding gene.
The priority for shared proteins is **Unsafe → Neither → CTA**. An unclassified
encoding gene is Neither: it can block a target assignment without creating a
forbidden assignment. Gene-level expression does not identify expressed isoforms.

The default forbidden union contains `unsafe.fasta.gz` and
`additional-forbidden.fasta.gz`. The latter contains incomplete/ambiguous
translations from Unsafe genes only. Neither, Variable and Unassessed are
excluded from this union. The previous Normal + Unassessed recipe is historical.

Full proteins are disjoint, but k-mers may overlap; subtraction is still required.
Never enumerate across stops, ambiguous residues or record boundaries. Missing
complete translations remain coverage gaps. The report gives downloads,
checksums, reference identity and the machine-readable exclusion recipe.

For example, MAGEA3 now fails at gene level because of high colonic cell RNA.
Its shared complete protein consequently enters Unsafe even though another
encoding locus is unresolved: Unsafe takes precedence over Neither.
XAGE1A/B and CTAG1A/B retain their shared CTA protein sequences.
