# Current normal-expression breadth

Descriptive policy `oncoref.normal-expression-breadth.2026-11-04.1` uses bulk tiers from
`oncoref:tissue-hazards:v1:sha256:3370cb388d11e33820e0332bb224eac576cb1ff9bd7e7536db1cae6dd212cbda`. **No admission cutoff or membership change is introduced.**
The report covers 20,099 distinct member loci in the qualified reference
universe, each once; these are gene profiles, not protein-group counts.

Skin, spleen and thymus now resolve to **Protected**. The historical October panel
remains available as `historical_panel()`; `anatomical_panel()` uses the registry.
The non-retinal eye residual uses the bulk eye rule. Cell exceptions never lower
a whole-organ observation. Each organ has equal weight within its tier; brain
regions share one organ's weight. Missing compartments retain their fixed weight.

For each tier, **B = sum(w × x / (x + L))**. Reproductive expression does not enter
this calculation and cannot improve it. The separate reproductive concentration
score `mean(N / (N + R))` answers a different question.

Current descriptive baseline ceilings are {'vital': 4.0, 'protected': 8.0, 'manageable': 16.0} nTPM, without
gene-specific IHC or clinical-credit adjustments. These anchors are research
choices, not injury thresholds. The [before/after table](before-after.csv.gz)
separates historical tiers/ceilings, current tiers with historical ceilings,
and current tiers/ceilings. [All scores](breadth.csv.gz) include brain-weight
sensitivity 0.5, 1 and 2. [Panels](panels.csv) expose every weight and assessment ID.

The [current](current-profiles.csv.gz) and [historical](historical-profiles.csv.gz)
profiles retain native source, version, unit, measurement and missing coverage.
Lower/upper bounds describe missing anatomical weight, not confidence intervals.
An entirely missing profile remains [0, 1]. HPA consensus is a published aggregate,
not an independent donor cohort. Duplicate observations cannot increase weight.

[Donor summaries](donor-summaries.csv) and [donor burdens](donor-breadth.csv) remain
separate. The available four-donor CTAG2-heart series supports a median but not
the prespecified ten-donor q90; missing q90 is not zero. Marginal bootstrap
intervals do not provide a joint tissue interval or population-risk estimate.

Rebuild with `python scripts/build_current_normal_breadth.py`. The
[receipt](manifest.json) pins source content, ontology, ceilings, code and outputs.
The [October report](../normal-breadth-20261026/index.md) stays frozen and reproducible.
