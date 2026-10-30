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

The same calls accept `dog` and `mouse`. `category` is CTA, Tumor-restricted, Unsafe, Neither, or null. `cta_selection.candidates()` returns both candidate subtypes. Null is an assessment status, not a fourth biological class. `assessment_status="screened"` means a category was assigned, not that safety was established or every source was complete. `incomplete_evidence` means required normal coverage or evidence resolution is missing; `normal_evidence_hold` means a normal-screen hold remains. Missing tumor or developmental evidence alone does not prevent Neither. `unassessed` denotes native nominations without an applied species policy. [Status definitions and per-protein gaps](index.md#what-the-assessment-statuses-mean) are in the report. Coverage counts include an `Unassigned` count so the denominator reconciles.

Human entities are proteins (identical selected reference sequences collapsed) with all member ENSG IDs and curated synonyms/history. Dog entities are native nominations with their actual identity status. Species-specific units, references and coverage remain explicit. `sources()` inventories measurement roles; `audit.sources` counts all published list annotations; `primary_source` selects the compact source view and `covered_by` records retained coverage of retired small sets. Every source table retains its citation/URL.

`audit.explain()` accepts an exact entity ID, gene ID or case-insensitive display label. RNA explanations contain the actual tissue, cell type, value, native unit, tier, single ceiling and IHC-match scope. HPA protein and imported categorical observations remain on their own scales. Unknown queries fail rather than fuzzy-matching a different gene.

Raw data share a version selector across the catalog:

```python
from oncoref import catalog, gtex
catalog.versions("gtex_tissue_rna")             # v11, v10, v8
catalog.ensure("gtex_tissue_rna")               # newest registered: v11
catalog.ensure("gtex_tissue_rna", version="v10")
catalog.provenance("gtex_tissue_rna", "v11")
gtex.tissue_rna("v11")                          # native matrix + provenance
```

[Source versions, historical caches and species-specific data](../../data-versions.md) explains the common contract. This saved analysis pins GTEx v11, HPA normal panels v23, Tyser and Petropoulos snapshots. Public HPA RNA summary loaders default to v25.1; that does not silently migrate this saved analysis. The raw GTEx reader preserves all 68 native contexts; this whole-organ adapter uses 52 and separately lists the 14 microdissected and two cultured/transformed contexts outside its scope.

All results are packaged and checksum-verified. Frozen release audits remain separately accessible; changing the current policy does not rewrite their categories or figures.
