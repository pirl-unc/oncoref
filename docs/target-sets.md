# CTA, Unsafe and Neither

The [current report](audits/target-sets-20261023/index.md) supplies three disjoint
sets across 20,097 canonical human protein-coding genes and the pinned Ensembl
115 protein inventory.

| Category | Target | Forbidden background |
|---|---|---|
| CTA | Yes, under the research screen | No |
| Unsafe | No | Yes |
| Neither | No | No |

CTA requires the previous curated candidate screen plus reproductive RNA
**≥5 nTPM when there is no positive IHC anywhere measured**. Reliable positive
reproductive IHC can support a lower-RNA candidate under the existing rules.
No IHC and measured-but-negative IHC are separate evidence states; neither is
positive evidence. This removes 150 candidates with no IHC and IL5 with
all-negative measured IHC, leaving **589 CTA genes**.

Unsafe requires positive normal IHC, RNA above an applicable normal-tissue
ceiling, or an independent reviewed safety hold. A selectivity failure alone
does not make a gene Unsafe. The existing category/IHC ceilings are unchanged;
5 nTPM is not a universal normal-tissue allowance.

Neither preserves Low, Variable, Unassessed and other unresolved/nonselected
subgroups. CTAG2 is Neither/Variable. ZNF730 and IL5 are Neither/Low. Missing
measurements never become zero. Variable remains excluded from targets while
its donor evidence and reconsideration requirements stay visible.

These are research classifications: CTA is not proven clinical safety, Unsafe
is a screening exclusion rather than proven toxicity, and Neither is not proof
of normal-protein absence. Only CTAG2 has imported individual-donor RNA so far.

## Python

```python
from oncoref import target_sets

evidence = target_sets.results()
cta = target_sets.gene_ids("CTA")
forbidden = target_sets.gene_ids("Unsafe")
neither = target_sets.gene_ids("Neither")
names = target_sets.gene_names("CTA")  # known symbols; IDs include unnamed loci
```

All three ID sets are disjoint. Earlier `cta`, `cta_screen` and normal-expression
exports retain their dated definitions. Use `target_sets` and the current report
for this policy.

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

For example, MAGEA3 passes at gene level, but its identical complete protein also
occurs at unresolved alternate-contig gene ENSG00000292186. That protein remains
Neither pending identity review; its gene name is not used to invent a mapping.
XAGE1A/B and CTAG1A/B retain their shared CTA protein sequences.
