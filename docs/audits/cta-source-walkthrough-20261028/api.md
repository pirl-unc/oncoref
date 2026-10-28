# Shared CTA selection API

```python
from oncoref import cta_selection

results = cta_selection.results(species="human")
selected = cta_selection.results(species="human", category="CTA")
unassigned = results[results.category.isna()]
audit = cta_selection.audit(species="human")
audit.sources
audit.overlaps
audit.funnel
audit.criteria
audit.explain("PRAME")
```

The same calls accept `dog` and `mouse`. `category` is CTA, Unsafe, Neither, or null. Null is an assessment status, not a fourth biological class. `assessment_status="screened"` means a category was assigned, not that safety was established or every source was complete. `incomplete_evidence` means an unassigned group lacks required coverage, a tumor/developmental measurement or evidence resolution; `outside_category_definitions` means those tracked checks pass but no category definition fits. `unassessed` denotes native nominations without an applied species policy. [Status definitions and per-group gaps](index.md#what-the-assessment-statuses-mean) are in the report. Coverage counts include an `Unassigned` count so the denominator reconciles.

Human entities are protein-equivalence groups with all member ENSG IDs and curated synonyms/history. Dog entities are native nominations with their actual identity status. Species-specific units, references and coverage remain explicit. `sources()` inventories measurement roles; `audit.sources` counts published list annotations. Every source table retains its citation/URL.

`audit.explain()` accepts an exact entity ID, gene ID or case-insensitive display label. RNA explanations contain the actual tissue, cell type, value, native unit, tier, single ceiling and IHC-match scope. HPA protein and imported categorical observations remain on their own scales. Unknown queries fail rather than fuzzy-matching a different gene.

All results are packaged and checksum-verified. Frozen release audits remain separately accessible; changing the current policy does not rewrite their categories or figures.
