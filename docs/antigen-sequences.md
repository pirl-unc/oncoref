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
from oncoref import antigen_sets, antigens

# Direct one-call subtraction (accepts a gene symbol or amino-acid sequence):
unique_9mers = antigens.subtract("MAGEA4", sex="female", max_mismatches=1, lengths=[9])

# Underlying set & sequence tables:
ctas = antigens.table("CTA")
cancer_antigens = antigens.table("CA")
candidates = antigens.table("CA|RA")
unsafe_gene_ids = antigens.gene_ids("Unsafe")
other = antigens.table("Other")
conservative_background = antigens.sequences("Unsafe|Other")
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

Sequence lookups automatically use the verified offline `sequence_bundle` (human) or packaged species proteome (`dog`, `mouse`) when staged locally (`sequence_source="auto"` or `sequence_source="sequence_bundle"`, probeable via `sequence_bundle.is_cached()`), requiring neither `oncoref[genome]` nor `pyensembl` and never triggering implicit downloads under `sequence_source="auto"`. Otherwise they fall back to the species' pinned Ensembl reference via `pyensembl` (`oncoref[genome]`, `sequence_source="pyensembl"`: Ensembl 115 for human/mouse, Ensembl 98 for dog). Translation
extraction preserves every annotated protein translation, including loss-of-function
and decay-associated transcripts and multiple occurrences of the same sequence.
`transcript_biotype` records the annotation for each occurrence. It inherits gene-level evidence, without
claiming isoform-specific expression. Missing translations appear in
`sequences.attrs["missing_gene_ids"]`; they are coverage gaps, not negatives.

```python
from oncoref import antigen_sets
from oncoref.peptide_subtraction import MatchPolicy, PeptideIndex, near_self_contract

# Direct offline index without pyensembl (supports species="human", "dog", or "mouse"):
background = PeptideIndex.from_sequence_bundle(background="Unsafe|Other")
dog_background = PeptideIndex.from_sequence_bundle(species="dog", background="Unsafe|Other")
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

## 1-mismatch Hamming near-self risk and tryptic-MS paralog auditing

Exact 9-mer uniqueness against `Unsafe` proteins (`unsafe_sequence_overlap_pct < 95%`) prevents targeting identical self-peptides, but paralogous cancer–germline families (`MAGE-A`, `PAGE`, `GAGE`) frequently differ from `Unsafe` somatic paralogs by only **one amino acid**—the exact mechanism behind fatal clinical TCR cross-reactivity (`MAGEA3` / `MAGEA12` `KMAELVHFL` vs. `KVAELVHFL`). Conversely, bottom-up healthy-tissue mass spectrometry digests proteins with trypsin and cannot distinguish isobaric isoleucine (`I`) from leucine (`L`), causing peptides shared with ubiquitous `Unsafe` paralogs to masquerade as healthy somatic translation of a CTA.

`oncoref.assay_boundary_adjudications` exposes both pure 1-substitution Hamming screening (`insertion_cost=2.0, deletion_cost=2.0`, forbidding indel window shifts at `max_distance=1.0`) and precomputed columns on `antigens.table()`, `target_selection.results()`, `antigen_sets.proteins()`, and `sequence_bundle.target_selection_dispositions()`:

- `exact_unique_9mers`: distinct candidate 9-mers absent at `0` mismatches from every `Unsafe` Ensembl 115 protein.
- `safe_2mm_9mers`: distinct candidate 9-mers at Hamming distance **≥ 2** from every `Unsafe` protein (`0` exact and `0` 1-mismatch `Unsafe` neighbors).
- `unsafe_near_self_1mm_pct`: percentage of `exact_unique_9mers` that sit within **1 amino-acid substitution** of an `Unsafe` protein (`81.4%` for `MAGEA6` [`86` exact → `16` `≥2-mm` safe against `MAGEA3`/`MAGEA12`]; `64.8%` for `MAGEA2/B` [`165` → `58`]; `42.9%` for `PAGE2` [`21` → `12`]).
- `unsafe_shared_tryptic_peptides`: canonical `I == L` tryptic peptides (`6–30 aa`) shared exactly with `Unsafe` Ensembl 115 proteins (`CTAGE1` shares `55 / 138` [`39.9%`] tryptic peptides with ubiquitous `Unsafe` `MIA2` [`49` peptides] and `CTAGE4/9`, explaining 5-organ healthy MS attributions without `CTAGE1` somatic transcription).

```python
from oncoref import antigen_sets, assay_boundary_adjudications
from oncoref.peptide_subtraction import PeptideIndex

# Pure 1-substitution Hamming index (no insertions or deletions) against Unsafe|Other:
bg_records = antigen_sets.forbidden_background("Unsafe|Other")
hamming_1mm_bg = PeptideIndex(
    bg_records,
    policy=assay_boundary_adjudications.hamming_1mm_policy(),
    lengths=(9,),
)

# Inspect precomputed 1-mm near-self and tryptic-MS paralog metrics on CTAs:
ctas = antigen_sets.proteins("CTA").set_index("protein_group_label")
ctas.loc[
    ["MAGEA6", "MAGEA2/B", "PAGE2", "CTAGE1", "CTAG2", "PLAC1", "DUX4"],
    [
        "exact_unique_9mers",
        "safe_2mm_9mers",
        "unsafe_near_self_1mm_pct",
        "unsafe_shared_tryptic_peptides",
    ],
]

# Detailed I==L tryptic-peptide MS paralog audit (CTAGE1/MIA2, TGM4/TGM2, SPATA31A3/SPATA31C):
ms_audit = assay_boundary_adjudications.tryptic_ms_paralog_audit()
```

## Downstream epitope-selection contracts: closed coding families & isoform clearances

Downstream epitope-selection pipelines (`tsarina`) subtract non-target human translations when building `cta_exclusive_peptides()`. Two structured column sets prevent false peptide subtraction and over-broad clinical safety claims:

1. **Closed coding-family sibling preservation:** Coding gene families whose non-identical paralogs split across `CTA` and `RA only` / `CA only` (`CT47A1..12` + `CT47B1` vs. `RA only` `CT47C1`; `TSPY1..10` vs. `RA only` `TSPY9`; `SPANXN1..5`; `CSAG1` vs. `CA only` `CSAG2/3`; `NY-ESO-1` vs. `CTAG2`) carry `closed_homolog_family_id` and JSON-encoded Ensembl gene IDs in `co_silenced_coding_sibling_ids` (plus `closed_self_paralog_ids` for `unprocessed_pseudogene` arrays like `DUX4` / `DUX4L*`). Every gene in `co_silenced_coding_sibling_ids` is verified somatically silent (`CTA`, `CA only`, or `RA only`)—never `Unsafe` (`TSPYL1..6` are excluded) and never `Other:Review`—so downstream selectors can exempt co-silenced germline siblings from self-subtraction.
2. **Isoform- and epitope-scoped clinical clearance:** `CTAG2` (`LAGE-1`) encodes `SLLMWITQC` on `ENST00000369585` (`LAGE-1a`, `180 aa`) but **not** on `ENST00000247306` (`LAGE-1b`, `210 aa`). Catalog and `sequence_bundle.target_selection_dispositions()` rows carry `clinical_clearance_scope`, `clinical_cleared_epitopes`, `clinical_hla_alleles`, `clinical_cleared_organs`, `clinical_cleared_transcripts`, and `clinical_unverified_isoforms` so epitope callers only apply clinical lung/bladder clearance to verified pMHC complexes.

```python
import json
from oncoref import sequence_bundle

disp = sequence_bundle.target_selection_dispositions().set_index("protein_group_label")
ctag2 = disp.loc["CTAG2"]
assert ctag2["clinical_clearance_scope"] == "epitope"
assert json.loads(ctag2["clinical_cleared_transcripts"]) == ["ENST00000369585"]
assert json.loads(ctag2["clinical_unverified_isoforms"]) == ["ENST00000247306"]
```

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
