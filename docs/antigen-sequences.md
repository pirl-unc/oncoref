# Antigen sets and peptide subtraction

`antigen_sets` reads the current catalog. CTA is CA ∩ RA; CA includes cancer-only
proteins and CTA, RA includes reproductive-only proteins and CTA, and `CA|RA`
returns **candidate antigen proteins**. The union name describes a research
inventory: reproductive-only proteins need tumor evidence before tumor targeting.
`Unsafe` is the accepted-exclusion background. `Other` contains unresolved or
nonqualifying proteins, including technical gaps and unmeasured home tissues; it
does not assert absent normal expression. `Unsafe|Other` is the conservative
background for screening shared peptides against both categories.

```python
from oncoref import antigen_sets

ctas = antigen_sets.proteins("CTA")
cancer_antigens = antigen_sets.proteins("CA")
candidates = antigen_sets.proteins("CA|RA")
unsafe_gene_ids = antigen_sets.gene_ids("Unsafe")
other = antigen_sets.proteins("Other")
conservative_background = antigen_sets.sequences("Unsafe|Other")
```

For host-specific work, select once and use that object's views and background:

```python
from oncoref.host_context import HostContext

selection = antigen_sets.select("CA", host=HostContext(protect=("testis",)))
proteins = selection.proteins()
sequences = selection.sequences()
background = selection.background_index("Unsafe|Other")
```

The background adds every host-excluded or unresolved candidate to the selected
Unsafe/Other background, including biological CTAs. Missing translations stop index
construction. [Host assumptions, scope and consequences](host-selection.md).
Convenience `proteins`, `gene_ids` and `sequences` also accept matching profile
arguments; the selection object avoids repeating those choices.

Rows carry the current policy and evidence. The older `target_sets` and
`cta_selection` APIs remain explicitly versioned historical snapshots.

## Sequences, occurrences and retained peptides

Install `oncoref[genome]` and provision human Ensembl 115 locally. Translation
extraction preserves every annotated protein translation, including loss-of-function
and decay-associated transcripts and multiple occurrences of the same sequence.
`transcript_biotype` records the annotation for each occurrence. It inherits gene-level evidence, without
claiming isoform-specific expression. Missing translations appear in
`sequences.attrs["missing_gene_ids"]`; they are coverage gaps, not negatives.

```python
from oncoref.peptide_subtraction import PeptideIndex, MatchPolicy

unsafe = antigen_sets.sequences("Unsafe")
targets = antigen_sets.sequences("CA")  # or CTA, RA, CA|RA
background = PeptideIndex(unsafe)       # background lengths 8–15, exact matches
# To include unresolved and nonqualifying proteins:
# background = PeptideIndex(antigen_sets.sequences("Unsafe|Other"))
target_index = PeptideIndex(targets)

peptide = "ACDEFGHIK"                 # illustrative sequence
target_occurrences = target_index.matches(peptide)
unsafe_occurrences = background.matches(peptide)
retained = background.subtract(targets.iloc[0].sequence)
```

`matches` retains all source fields (gene, transcript, protein group, catalog category
and background exclusion reason), the original
matched peptide, distance, and zero-based half-open coordinates. `subtract`
returns each surviving query occurrence with its coordinates. Both retain the
reference, missing-gene list and matching parameters in DataFrame attributes.
They never join across proteins, stop codons or ambiguous residues. By default,
every 8–15mer is considered; `lengths=[9]` makes either operation length-specific.
The background's configured lengths determine which potential matches are tested.

## Optional similarity screens

```python
# Ordinary Levenshtein distance: substitutions, insertions and deletions.
near_self = PeptideIndex(unsafe, policy=MatchPolicy(max_distance=1))

# Reduced alphabet, still requiring zero distance after equivalence.
ilv = PeptideIndex(unsafe, policy=MatchPolicy(equivalence_groups=("ILV",)))

# Explicit, illustrative weights; unspecified substitutions cost 1.
weighted = PeptideIndex(unsafe, policy=MatchPolicy(
    max_distance=0.5,
    substitution_costs={( "I", "V"): 0.25, ("V", "I"): 0.25},
    insertion_cost=1,
    deletion_cost=1,
))
```

Thresholds are inclusive. Substitution costs are directed from query to background;
provide both directions for a symmetric rule. All edit costs must be positive and
finite; only explicit disjoint equivalence groups have zero cost. The search checks
different background lengths when indels permit them. It uses exact seeds followed
by weighted dynamic programming, and falls back to exhaustive windows when the
edit budget leaves no guaranteed seed. Large budgets or tiny costs can be expensive.
No proteome-wide 8–15mer occurrence table is materialized.

Use exact subtraction as the reproducible baseline and compare I/L or I/L/V as
optional sensitivity screens. Other groups can be supplied, but none is enabled
implicitly. Physicochemical grouping has found cross-reactive candidates in a
[primary TCR screening study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10371023/);
it is not a validated universal equivalence relation. TCR recognition also depends
on peptide position, HLA and receptor. Neither a retained sequence nor a low
distance establishes presentation, cross-reactivity or therapeutic safety.

## Reviewed findings

```python
from oncoref import gene_knowledge
magec2 = gene_knowledge.findings("MAGEC2", include_updates=True)
daz = gene_knowledge.findings("DAZ3", include_updates=True)
```

A later finding may explicitly `supersede` an earlier finding for the same gene
identities. Both remain in the provenance; only active findings impose admission
holds. MAGEC2's historical TC-1 cerebellar signal is adjudicated as unsupported
using HPA062230 and GTEx, without claiming the original signal has been
experimentally disproved. The independent MAGE family filter still applies.
DAZ findings preserve the gastric RNA evidence, uncertain individual-paralog
assignment, absent protein confirmation and Unsafe sequence containment rationale.
