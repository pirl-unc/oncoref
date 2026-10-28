# A common source and decision audit API

```python
from oncoref import cta_selection, human, canine

audit = cta_selection.audit(species="human")
audit.funnel          # cumulative categories and sequential candidate removals
audit.criteria        # actual rules, source IDs and evaluation status
audit.sources         # per-paper candidate-list counts and citations
audit.memberships     # source-to-entity mappings, member gene IDs and category
audit.overlaps        # all pairwise mapped and selected overlaps
audit.explain("CTAG2") # tissue, cell, value, unit, threshold and source URL

human.cta_selection_audit()   # the same audit through the species facade
dog = canine.cta_selection_audit()
dog.status                   # "unassessed"
dog.funnel.eligible          # nulls: a quantitative dog screen has not run

# Existing interfaces retain their meaning:
cta_selection.results(species="human", category="CTA")
cta_selection.sources(species="dog")   # all measurement sources and integration states
cta_selection.coverage(species="mouse")
```

`SelectionAudit` returns detached pandas tables from packaged, checksum-verified data. It works offline and does not rerun biological classification. `explain()` accepts an exact entity ID, member stable gene ID or case-insensitive display label. It raises for unknown or ambiguous entities in a populated audit; an unavailable mouse audit returns an empty table. Use the existing identity-resolution interfaces for alias and historical-ID reconciliation before querying an audit. This method does not introduce another synonym resolver.

| Shared element | Species-specific detail |
|---|---|
| `entity_id`, `gene_ids` | Human protein-equivalence groups; canine native gene/assay nominations, with the MAGE family ambiguity retained |
| `sources`, `memberships`, `overlaps` | Per-paper unions for humans; two native candidate studies for dogs; no integrated mouse nomination list |
| `funnel`, `criteria` | Completed human quantitative policy; explicit unassessed dog stages; unavailable mouse curation |
| `explanations` | Source/tissue/cell measurements and native units when available; coverage reasons when unassessed |
| `status`, `policy_id` | A human release policy ID; null policy for species without a quantitative decision policy |

The public source categories remain common across species: candidate nomination, reproductive/developmental RNA, tumor RNA, normal bulk RNA, normal cell/immune RNA, normal protein, clinical counterevidence and coverage. Source adapters retain native units, reference identity, ages, tissues, donors/pools and ambiguity. A common API does not imply equal source coverage or transferable thresholds.

`sources` inside the audit counts **published candidate-list membership**. The existing `cta_selection.sources()` inventories **all evidence resources**, including unavailable or metadata-only sources. Paper counts overlap and do not measure causal admission. Human counts use protein-equivalence groups; dog nomination counts must not be compared as a genome-wide denominator.

An explanation provides a representative supported context per RNA layer plus additional compartment observations for widely studied antigens. These observations link the complete frozen evidence exports. A pass means the layer did not exclude the entity under the current policy; it does not mean measured zero or complete coverage. Unavailable quantitative summaries remain null. Negative IHC in the same tissue is reported separately if the exact cell match remains unresolved. The cell-context layer can also exclude through matched positive IHC when supported RNA is missing; the observation unit distinguishes staining from RNA. Retained observations use `flag` or `observation`, separate from `hold` and `unsafe` decisions.

## Reproduction

```sh
.venv/bin/python scripts/build_cta_selection_audit.py
.venv/bin/python scripts/render_cta_source_walkthrough.py
```

The build reads the frozen 1.21.0 report, pinned native candidate curation and cached HPA v23 evidence. The packaged receipt records input/output SHA-256 hashes. The report renderer reads the shared API and exports the figures, citation tables and explicitly labeled bowel comparison. It does not alter `cta_selection.results()` or the released decisions.

[Human walkthrough](index.md) · [Dog funnel](dog.md) · [Policy discussion](policy-review.md).
