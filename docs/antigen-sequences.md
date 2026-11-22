# Antigen sets and peptide subtraction

Sequence rows cover all annotated translations of selected member loci. Use `matches_representative_sequence` for exact amino-acid matches to a catalog representative. A custom selection that omits reference hashes yields a missing match value. Tumor RNA measurements are at locus or group resolution.


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
background = selection.background_index()  # defaults to Unsafe|Other
```

The background adds every host-excluded or unresolved candidate to the selected
Unsafe/Other background, including biological CTAs. Missing translations stop index
construction. [Host assumptions, scope and consequences](host-selection.md).
Convenience `proteins`, `gene_ids` and `sequences` also accept matching profile
arguments; the selection object avoids repeating those choices.

Rows carry the current policy and evidence. The older `target_sets` and
`cta_selection` APIs remain explicitly versioned historical snapshots.

## Sequences, occurrences and retained peptides

Human sequence lookups automatically use the verified offline `sequence_bundle` when
staged locally (`sequence_source="auto"` or `sequence_source="sequence_bundle"`),
requiring neither `oncoref[genome]` nor `pyensembl`, and fall back to local Ensembl
115 (`oncoref[genome]`, `sequence_source="pyensembl"`) otherwise. Translation
extraction preserves every annotated protein translation, including loss-of-function
and decay-associated transcripts and multiple occurrences of the same sequence.
`transcript_biotype` records the annotation for each occurrence. It inherits gene-level evidence, without
claiming isoform-specific expression. Missing translations appear in
`sequences.attrs["missing_gene_ids"]`; they are coverage gaps, not negatives.

```python
from oncoref import antigen_sets
from oncoref.peptide_subtraction import MatchPolicy, PeptideIndex, near_self_contract

# Direct bundle-backed index without pyensembl:
background = PeptideIndex.from_sequence_bundle(background="Unsafe|Other")
contract = near_self_contract()

# Or explicit forbidden-background table with provenance & uncertainty attrs:
background_records = antigen_sets.forbidden_background("Unsafe|Other")
targets = antigen_sets.sequences("CA")  # or CTA, RA, CA|RA
target_index = PeptideIndex(targets)

peptide = "ACDEFGHIK"                 # illustrative sequence
target_occurrences = target_index.matches(peptide)
background_occurrences = background.matches(peptide)
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
near_self = PeptideIndex(background_records, policy=MatchPolicy(max_distance=1))

# Reduced alphabet, still requiring zero distance after equivalence.
ilv = PeptideIndex(background_records, policy=MatchPolicy(equivalence_groups=("ILV",)))

# Explicit, illustrative weights; unspecified substitutions cost 1.
weighted = PeptideIndex(background_records, policy=MatchPolicy(
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

## Germline coding polymorphisms, isoform selectivity and ncORF/TE contracts

Reference proteomes represent a single haplotype. Common germline missense variants (`AF >= 0.001`) can create self-peptides in healthy patients that are absent from the reference FASTA. Augment any background table or `PeptideIndex` with `CodingVariant` records (validated fail-closed against 1-based reference coordinates):

```python
from oncoref.peptide_subtraction import CodingVariant

variants = [
    CodingVariant(
        target_id="ENST00000370225",
        position=15,
        ref_aa="A",
        alt_aa="V",
        af=0.042,
        source_id="gnomAD:v4.1",
    )
]
polymorphic_bg = background.with_coding_variants(variants, min_af=0.001)
```

For tumor-restricted splice isoforms on genes whose canonical isoform is `Unsafe`, gene-level `Unsafe` precedence holds unless healthy normal tissues are quantified at transcript/junction resolution **and** subtracting against the full `Unsafe` proteome leaves unique junction/exon peptides (`antigen_sets.evaluate_isoform_selectivity` / `antigen_sets.isoform_selectivity_contract`). Non-canonical ORFs and transposable-element (`HERV`/`LINE-1`/`LTR`) nominations follow the separate four-gate contract in `sequence_bundle.ncorf_te_screening_contract()`.

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
