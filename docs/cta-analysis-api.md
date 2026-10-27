# CTA analysis APIs and identity conventions

The active research report is the [harmonized CTA funnel](audits/cta-harmonized-20261027/index.md).
Earlier membership and threshold reports remain versioned comparisons; their
`results()` accessors return frozen datasets, not the latest research selection.

## Names and counting units

| Field | Meaning |
|---|---|
| `gene_id` / `canonical_gene_id` | One canonical human Ensembl gene locus. `canonical_gene_id` is used in source-match audits. |
| `source_gene_id`, `source_gene_name` | Original source identifiers, retained even if unresolved. |
| `historical_ensembl_gene_ids` | Registered historical replacements; a JSON list, not a guessed name match. |
| `alternate_ensembl_gene_ids` | Alternate-contig identifiers, kept distinct from historical replacements. |
| `gene_name_synonyms` | Registered gene names and synonyms, including ambiguous names where applicable. |
| `protein_equivalence_group_id` | Membership-based `PEG:` identifier for identical chosen representative proteins. It is not an ENSG or a sequence hash. |
| `protein_group_label` | Display label such as NY-ESO-1 or XAGE1A/B; never the gene identity. |
| `representative_gene_id` | A real member locus used for compatibility with existing measurement APIs. It does not represent an extra gene. |
| `member_gene_ids`, `member_identities` | Complete expected group membership and identity bundles, including explicit unresolved registry members. |
| `proteoform_key` | Existing compatibility display key from `oncoref.proteoforms`; new research tables also expose the explicit group ID. |

“Protein equivalence” here means identical selected reference sequences. It does
not establish equivalence of all splice isoforms, post-translational forms or
gene-specific regulation. Reports always label whether counts represent loci or groups.

## Processing boundaries

`gene_matching` handles identity, `protein_equivalence` handles additive measurements,
`cell_rna_policy` handles context-specific RNA thresholds, and `cta_funnel` handles
independent decision layers. Report scripts assemble these APIs and render artifacts.

```python
from oncoref.gene_matching import match_reference_genes, harmonize_reference_expression
from oncoref.protein_equivalence import collapse_reference_expression
from oncoref.cta_funnel import EvidenceLayer, classify_layers

matches = match_reference_genes(
    source_table, id_column="Gene", symbol_column="Gene name",
    source="source-name", source_version="pinned-release",
)

collapsed = collapse_reference_expression(
    sample_matrix, sample_cols=sample_ids, source="cohort-name",
    source_version="pinned-release", unit="TPM", scope="genome",
)
values = collapsed.expression
coverage = collapsed.coverage
identity_audit = collapsed.identity_audit

decisions = classify_layers(
    group_facts,
    [EvidenceLayer("coverage", missing_coverage, no_unsafe_evidence),
     EvidenceLayer("normal RNA", above_candidate_ceiling, above_unsafe_boundary)],
    entity_column="protein_equivalence_group_id",
)
```

Identity matching quarantines conflicting stable IDs/names and ambiguous symbols.
An explicit unknown ENSG is never replaced by a name-only guess. Duplicate aliases
of one locus are not additive. The strict expression API raises on conflicts; the
explicit `quarantine` mode excludes conflicting measurements and records why.

Sums require the same sample/context and the same native unit. Partial sums can
establish positive lower bounds; incomplete member coverage cannot establish
normal absence. Maxima, percentiles and ranked DICE slots are not additive inputs.
`sum_context_observations()` supports long-form matched contexts;
`unpaired_rank_bounds()` preserves ranked-sum uncertainty without pairing donors.

## Staining evidence

`ihc_observations.parse_staining_statement()` extracts the reported result,
intensity, localization and percentage bounds. It retains the source statement.
`external_ihc.assess_evidence()` applies explicit evidence-quality requirements;
`reconcile_hpa()` reconciles exact gene/tissue/cell/stage/specimen-state contexts.
Collection, categorical result and validation confidence are separate concepts.
An unreviewed negative is still a negative *report*, with explicit reasons it
cannot yet clear an uncertainty hold. Cross-reactivity is retained as a source
claim; control recommendations are not observations. Reagent resellers do not
create independent evidence, and whole-organ negatives do not clear rare cell types.

All report writers share deterministic CSV/gzip output and SHA-256 helpers from
`scripts/cta_report_common.py`. Source, input and output receipts accompany each
stage. Figures are generated from the current decision tables, with no hardcoded
candidate counts in the active report.

`input_hashes()` pins the shared writer as well as the stage inputs. Repository
paths are relative to the checkout root; explicitly external source files retain
their original paths. The full raw source-matrix hashes are also preserved.

To rebuild the current report from the pinned source inputs and registered
matrices, run these commands from the repository root:

```sh
python scripts/cta_tumor_union.py docs/audits/cta-expanded-funnel-20261027/tumor-union
python scripts/cta_expanded_normal.py docs/audits/cta-expanded-funnel-20261027/normal
python scripts/cta_expanded_funnel.py docs/audits/cta-expanded-funnel-20261027
python scripts/cta_protein_group_normal.py docs/audits/cta-harmonized-20261027/group-normal
python scripts/curate_external_ihc.py
python scripts/cta_harmonized_report.py docs/audits/cta-harmonized-20261027
```

The external curation command uses the checked-in categorical source extraction;
`scripts/collect_external_ihc.py` refreshes that extraction when an intentional
source update is required. Raw normal-expression downloads use the source
registry and local caches. Different upstream releases require a new audit.

Group decision exports contain only group-level quantities. They do not copy a
representative locus's old maxima or selection flags into group measurement
columns. `Expression category`, `Expression detail` and `Expression forbidden`
in the locus and concordance exports are explicitly historical annotations;
`category` and the named `hold:` / `unsafe:` columns describe the current proposal.

The current policy uses unusable or disputed antibodies to disable negative-IHC
threshold boosts, without an additional gene-wide antibody veto after the strict
RNA fallback. `Antibody concern` and `external_review_flag` remain evidence flags
even for selected groups. Reliable positive evidence still contributes exclusions.
Clinical target-context evidence is recorded separately; an expression-policy
`Unsafe` label does not establish treatment toxicity. The
[PRAME, CTAG2 and TBXT assessment](audits/cta-target-evidence-20261027/index.md)
exports the relevant measurements, thresholds and clinical-source denominators.
