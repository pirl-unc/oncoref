# Continuous normal-expression breadth

Descriptive policy `oncoref.normal-expression-breadth.2026-10-26.1`; active eligibility remains `oncoref.target-sets.2026-10-25.2`.
All 2,540 mapped nominees and three prespecified ordinary-normal controls are
represented. Scores do not admit or exclude any gene. Clinical benchmark genes
were chosen before scoring; this is a descriptive comparison, not validation.

For each category, `B = sum(weight * RNA / (RNA + ceiling))`. Vital, Protected
and Manageable are separate. Reproductive abundance never enters this formula.
The [complete profiles](applicable_ihc-profiles.csv.gz) use the current archived
per-gene IHC-adjusted RNA ceilings. A [fixed-baseline comparison](fixed_baseline-profiles.csv.gz)
uses 2/6/10 nTPM for every gene, including ACTB, GAPDH and PPIA controls.
[All scores and weighting sensitivity](breadth.csv.gz) · [Benchmarks](benchmarks.csv).

## Fixed anatomy and missing coverage

The [panel](anatomical-panel.csv) is specified independently of each gene's
observed tissues. Each organ has equal weight within its category; eight brain
regions divide one brain weight. Retina and unmeasured non-retinal eye divide
one eye weight. Broad small intestine is used once, without adding overlapping
duodenum. Bronchus and embedded vascular/immune populations are not counted as
additional independent whole organs. Cellular blind spots still require the
separate cell-expression audit. This panel is not whole-body cell coverage.

Blood, peripheral nerve, inner ear, non-retinal eye and several protected
compartments are explicit missing coverage where this bulk source has no sample.
Observed weighted contribution is the lower bound; adding missing weight gives
the upper bound. These are coverage bounds, **not confidence intervals**. An
all-missing gene gets [0,1], never a negative finding. Six nominees absent from
the archived tissue screen retain their missing profiles in that scenario.
Sensitivity halves/doubles the prespecified brain organ weight, conserving
within-brain weights and renormalizing the fixed panel, never just observed rows.

## Donors, sources and units

HPA consensus is a published cross-source summary, not a donor median or an
independent cohort. Individual-donor counts and sampling intervals cannot be
inferred from it. [Available reviewed donor summaries](donor-summaries.csv)
currently cover CTAG2 heart only: four independent donors permit a median;
the prespecified ten-donor minimum leaves q90 unavailable. Separate
[donor burdens](donor-breadth.csv) and [profiles](donor-profiles.csv.gz) use the
fixed baseline ceilings, because these RNA donors are not paired to IHC. The
unavailable q90 retains entirely missing coverage, never a zero measurement.
A marginal 95%
percentile-bootstrap interval uses 2,000 donor resamples with seed 0; it is
reported separately from anatomical coverage and is unstable with four donors.
This tiny selected
sample is not an estimate of population risk, nor matched RNA/IHC specimens.
No whole-profile sampling interval is claimed. Native TPM/CP10K, cell atlases,
other species and HPA versions are excluded pending their own calibration.

## Interpretation and reproduction

No classification cutoff has been fitted. Accepted focal protein findings,
hard RNA exclusions and current holds remain in the eligibility policy; a low
burden cannot erase them. This score is a review-priority descriptor, not injury
probability or a body-wide expression fraction.

Run `.venv/bin/python scripts/build_normal_expression_breadth.py --out <empty-directory>`.
The manifest pins inputs, implementation, source receipt and all outputs.
