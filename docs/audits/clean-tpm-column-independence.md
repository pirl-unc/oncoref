# Clean TPM column independence and partial matrices

Review for [oncoref #541](https://github.com/pirl-unc/oncoref/issues/541),
2026-09-22. The issue's final comment supersedes its original gene-subset
diagnosis: it reports a biological-compartment change when companion columns
are added. Missing technical genes do not explain that particular observation.

## Reproduction limits

At baseline `bbc97d2`, the local 19,781-row pan-cancer reference with all 121 raw
columns sums to 1e6 per column (within floating-point precision) under Python
3.12.6, pandas 3.0.5, and NumPy 2.5.2. The adrenal column agrees with its
single-column normalization. Checks of the 116 directly sourced/raw companion
columns also pass with pandas 2.2.3 / numexpr 2.13.1 and pandas 3.0.5 / numexpr
2.14.2. **The reported multi-percent drift was not reproduced**, and this change
does not claim to identify its historical dependency-level cause.

A separate precision defect is reproducible: the previous code converts the
output to float64 but sums the original biological input, including float32
blocks. A wide mixed-dtype matrix changes the reduction layout relative to
normalizing one column. The new float32 regression fails on the baseline
(maximum observed per-gene difference from a float64 calculation: about 0.0011).
This small error is not evidence for the much larger discrepancy in #541.

## Resolution and contract

The shared compartment scaler converts to float64 **before** summation, uses
contiguous columns and NumPy scaling, and checks the filled compartment budget.
Each column therefore uses the same arithmetic regardless of its companions,
input block layout, or pandas' optional expression engine. The censored
reference composition and 16/9/75 budgets are unchanged.

If a compartment cannot carry its budget, `RuntimeWarning` identifies it and
the affected columns. Missing values stay missing and the other compartments
do not absorb the unfilled share. A gene subset with no ribosomal rows still
sums to 840,000 when its other compartments are available; it now announces why.

Partial matrices use the supplied measured gene universe. There is no invented
whole-transcriptome mass. To retain a reference's normalization for a panel,
normalize before subsetting. `pan_cancer_expression(genes=...)` already follows
that order; its filtered result is not required to sum to 1e6.

Regression coverage in `tests/test_normalization.py` compares 121-column float32
and float64/mixed frames to independent per-column calculations, including
missing measurements, reordered/subselected columns, input immutability,
scale invariance, and unavailable compartments. `tests/test_expression.py`
checks the public pan-cancer budgets and verifies that gene filtering preserves
the full-reference values. The original failing environment remains useful for
confirming the historical reproduction.

## Downstream confirmation

On 2026-09-23, both pirlygenes guards named in #541 passed against this checkout
with `pytest --runxfail`, so expected-failure markers could not hide a failure:

- `test_pan_cancer_expression_normalize_tpm_clean_fixed_fraction`
- `test_pan_cancer_expression_normalize_tpm_clean_pins_cols_to_million`

The run used pirlygenes' Python 3.12.6 virtualenv, pandas 3.0.5, and NumPy 2.5.2.
Its checkout was `cefb19ce20ede3aa3d6e3d563eef9c31136930ba`; local changes to the
test file added the issue commentary and non-strict xfail decorators, leaving
the two test bodies unchanged. This confirms the current downstream contract,
not the root cause of the historical multi-percent drift. The test file and
its local changes were left intact.
