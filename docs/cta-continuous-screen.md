# Continuous CTA screening

The [23 October ordered-rule screen](audits/cta-screen-20261023/index.md) adds
reproductive positive controls, tighter selectivity and donor-based exclusions.
Its **Variable** category excludes CTAG2 while recording the evidence needed for
reconsideration. The 22 October report below remains a historical pilot.

The [22 October pilot report](audits/cta-continuous-20261022/index.md) applies a
continuous low-expression allowance to the published HPA v23 screen. It keeps
two outputs separate:

| Output | Meaning |
| --- | --- |
| Point screen | Do the recorded RNA values and existing curation/IHC gates pass? |
| Uncertainty result | Does a separately fitted donor model support inclusion, exclusion, or further review? |

The current pilot has **no fitted donor model**. Point-screen passes therefore
remain **Review** in the uncertainty assessment. Frozen public panels and
sequence exports are still governed by their named release policies.

## Continuous rule

For each screened tissue, use `N <= C`, `N <= ((1-s)/s)*R + delta`, and `R > N`.
Here `R` is the reproductive maximum, `N` is normal-tissue expression, and `C`
and `s` are the existing category/local-IHC thresholds. The new fixed allowance
`delta` defaults to **0.1 nTPM**. It is an inspectable pilot parameter, not a
measurement uncertainty estimate or a validated limit. Larger uncertainty must
not increase it. All positive normal IHC, including Low, still vetoes inclusion.

```python
from oncoref.cta_continuous_screen import ContinuousPolicy, continuous_tissue_screen

# audit is the output of adaptive_normal_rna_limits(), or the pinned report table.
policy = ContinuousPolicy(allowance_ntpm=0.1)
continuous = continuous_tissue_screen(audit, policy=policy)
continuous[["tissue", "point_result", "continuous_selectivity_margin_ntpm",
            "absolute_margin_ntpm", "reproductive_margin_ntpm"]]
```

Negative margins identify failures. The reproductive margin must be strictly
positive; the other margins can be zero. The result includes the unchanged
source measurements, local IHC confidence, original policy ID and new policy ID.
Use `policy.manifest()` to record the exact numerical configuration, including
custom parameters. Never identify an analysis solely by policy ID when parameters
have been overridden.

Inclusive comparisons allow one floating-point representable step at the
boundary; the corresponding reported margin is zero. This fixes cases such as
0.3 / (0.3 + 0.1) evaluating just below 75% in binary arithmetic. It is separate
from the 0.1 nTPM project allowance. Positive IHC is resolved at its measured cell
compartment: permitted hair-follicle staining does not become a whole-skin veto,
and it still cannot establish complete negative bulk-skin IHC for relaxation.

This tissue function does not itself apply every gene admission or curation
rule. The report also requires measured reproductive expression, routine RNA
coverage, usable anatomy, family/nomination eligibility and no external hold.
Unmapped tissues and missing values are unassessed, not permissive zeros.

## Joint uncertainty interface

`assess_joint_expression()` accepts one gene's audit, independently specified
normal/reproductive tissue panels, and optional `JointExpressionEvidence`.
Omitting the evidence returns Review unless a hard exclusion applies.

```python
from oncoref.cta_continuous_screen import assess_joint_expression

decision = assess_joint_expression(
    audit,
    normal_tissues=expected_normal_tissues,
    reproductive_tissues=expected_reproductive_tissues,
    evidence=None,
    hard_exclusions=(),       # e.g. an external safety hold or family exclusion
    review_reasons=(),        # e.g. missing coverage or unresolved nomination
    other_gates_passed=False, # confirm only after the complete gene-level review
)
```

Panels must come from the source's independent coverage inventory, not this
gene's observed subset. The report uses the frozen 50-tissue HPA inventory and
the anatomy registry to distinguish reproductive, screened and disposable tissues.
Its scope is that measured panel; it does not certify all organs or cell types.
Human nTPM numerical limits do not transfer to mouse or dog measurements.

The RNA audit cannot reveal IHC positives in tissues without a corresponding RNA
row. Before setting `other_gates_passed=True`, check the complete IHC observations,
reproductive admission, coverage/anatomy, nomination/family decisions and external
holds. Pass the resulting hard exclusions and review flags explicitly. Leaving
this confirmation at its default prevents inclusion even with strong RNA draws.

### Required model outputs

`JointExpressionEvidence` contains:

- **`means`:** a DataFrame with independent posterior draw IDs as rows and tissue
  names as columns. Entries are draws of the tissue consensus mean on the
  **HPA v23 nTPM** scale. All required reproductive and normal tissues are needed.
- **`donor_upper`:** draws of each normal tissue's **95th donor percentile**,
  aligned to exactly the same draw IDs and order. This estimates biological donor
  variation rather than merely the sampling precision of a mean.
- **`provenance`:** `gene_id`, `model_id`, a lowercase 64-character `input_sha256`,
  `normalization="HPA-v23-nTPM"`,
  `estimand="posterior_tissue_consensus_mean"`, `donor_quantile=0.95`,
  `model_checks_passed=True`, and `independent_draws=True`.

The provenance assertions are supplied by the model fitter; this evaluator
does not establish calibration, convergence or correct normalization merely by
checking their presence. The input hash should identify a manifest that anchors
counts, sample/donor metadata, reference annotation, normalization, model code,
parameters and diagnostics. Record the posterior generation seed and settings.

A count-aware model must account for low counts, uncertainty after observed zeros,
between-donor variation, cohort differences and ambiguous gene assignment. A
precision interval around a mean cannot replace a donor distribution. Available
HPA sample transcript TPM/estimated counts require gene aggregation and verified
normalization; using unrelated TPM directly as nTPM is rejected. The
[HPA source definitions](https://v23.proteinatlas.org/about/download) and
[expression-unit guide](cta-expression-units.md) describe the scale. The
[apeglm paper](https://pubmed.ncbi.nlm.nih.gov/30395178/) supplies a statistical
precedent for stabilizing noisy count-based effects, not a validation of these
screening parameters.

### Decision calculation

Within every draw, recompute `R` from all reproductive tissue means and evaluate
all normal tissues together. This preserves the supplied correlations and
reference-selection uncertainty. The empirical passing fraction estimates the
posterior probability of meeting the mean-expression rules. It is not computed
by multiplying marginal tissue probabilities.

The default 95% probability cutoff, 95th donor percentile and minimum 1,000
independent draws are provisional implementation choices, not clinical criteria.
Monte Carlo estimates remain approximate; the fitter must establish adequate
effective sample size and check stability near a decision boundary.

| Result | Requirement |
| --- | --- |
| Included | At least 95% joint support for the mean rules **and** donor percentiles below their absolute ceilings; no review flags or hard exclusions. |
| Excluded | At least 95% support for failing the joint mean rules, or an explicit positive-IHC/curation/safety exclusion. |
| Review | Intermediate support, donor-tail concerns, missing evidence, insufficient draws, failed model checks or degenerate posterior inputs. |

The donor-tail check can only turn a possible inclusion into Review; it is not
used alone to assert clinical harm. Degenerate constant draws, including an
all-zero donor bootstrap, cannot establish uncertainty-qualified inclusion.
Unusable model evidence cannot justify an RNA-based exclusion either.

IHC remains an independently applied categorical gate and source of the existing
RNA thresholds. HPA reliability labels are not assigned invented sensitivities,
false-negative probabilities or numeric likelihoods. All probabilities describe
expression-rule compliance conditional on the supplied model, not treatment safety.

## Reproduce the current pilot

```bash
python scripts/cta_continuous_report.py --output /tmp/cta-continuous-audit
```

This verifies the historical report hashes, recomputes all tissue margins and
gene decisions, compares allowances, and writes a searchable report, CSV tables,
figures, policy and run manifest. It does not download or fit donor data. The
clinical panel is a retention/safety review check, not independent calibration
data. Any future donor-model integration needs a separately pinned, reviewed fit
and held-out checks before publishing uncertainty-qualified inclusion counts.
