# OncoRef 1.21.0

This release makes the rebuilt CTA/oncofetal research selection available through one species-scoped API. Human selection now contains **463 protein-equivalence groups / 481 loci**, with **17,388 Unsafe** and **2,012 Neither** groups across the complete reference universe.

```python
from oncoref import cta_selection
selected = cta_selection.results(species="human", category="CTA")
coverage = cta_selection.coverage(species="dog")
sources = cta_selection.sources(species="mouse")
```

The same schema, identity fields, category vocabulary and source categories apply to human, dog and mouse. Dog currently exposes 14 native nominations as Neither/unassessed. Mouse returns an empty typed table with unavailable coverage. Neither case substitutes human measurements or converts native counts/CPM/iBAQ to TPM.

The release preserves independent reasons from every filter, identity history and synonyms, matched protein-group RNA sums, optional DICE/cornea coverage, native units and uncertainty. The code review corrected 2,902 grouped IHC labels, required controls for every group member before negative-IHC credit, and replayed complete immune/corneal contexts. Eleven intermediate decisions changed; final categories remained unchanged from the preceding 463-group proposal.

PRAME, MAGE-A4, NY-ESO-1 (CTAG1A/B), XAGE1A/B, SSX2/B, CT83, TBXT, SMC1B and PAGE1 pass the uniform rules. PRAME retains its Low adrenal-protein and bowel-progenitor flags. TBXT retains external-IHC/adult-disc review context. CTAG2 remains Unsafe, including bronchial basal respiratory-cell RNA. Among 213 previously annotated CTA loci, the new funnel assigns **162 CTA / 16 Unsafe / 35 Neither**.

The [report](index.md), [literature review](literature-review.md), editable [16-slide presentation](slides/cta-selection-1.21.0.pptx), [rendered slide PDF](cta-selection-1.21.0.pdf), SVG/PNG figures and downloadable decision tables use the same audited results. The PDF is an image-based rendering with bookmarks and a report link on every page; the PPTX retains native tables and charts. [API documentation](../../cta-selection.md) explains the source inventory, provenance and migration.

These are research prioritization decisions, not clinically validated safety classifications. PRAME informed calibration; the cutoffs are not independent evidence of its safety. Adult data do not establish pediatric safety, and the integrated embryonic atlas covers one early embryo. The separate v1.20 exact-protein near-self background and historical compatibility APIs retain their own policy/reference identities; this release does not silently rewrite those backgrounds.

Literature retrieval is not complete primary-paper adjudication. The published audit records review depth and unassessed genes; the broader adjudication work in [#638](https://github.com/pirl-unc/oncoref/issues/638) remains open. The descriptive [normal-expression breadth analysis](../normal-breadth-20261026/index.md) is available separately and does not introduce a validated classification cutoff.

Validation: the full local suite passed 2,723 tests with no skips; three additional release-artifact tests passed after their addition. Lint, formatting, strict documentation build, package build, dependency check and installed-wheel checks for all three species passed. Isolated source-archive identity tests passed 52 cases; two real-genome tests deliberately skipped with its empty genome cache, and both ran in the full local suite. CI results are attached to [PR #642](https://github.com/pirl-unc/oncoref/pull/642).
