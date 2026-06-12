# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Expression normalization: clean TPM + helpers (#35, Phase N)."""

import numpy as np
import pandas as pd
import pytest

from cancerdata import gene_families as gf
from cancerdata import normalization as norm


def _matrix():
    tech = list(gf.gene_family_ids("mitochondrial"))[:2]
    gt = pd.DataFrame(
        {
            "Ensembl_Gene_ID": [*tech, "ENSG00000111111", "ENSG00000222222", "ENSG00000333333"],
            "Symbol": ["MT1", "MT2", "BIO1", "BIO2", "BIO3"],
        }
    )
    vals = pd.DataFrame(
        {"s1": [100.0, 200, 300, 400, 500], "s2": [10.0, 0, 100, 200, 300]}, index=gt.index
    )
    return gt, vals


def test_clean_tpm_two_compartment_budget():
    gt, vals = _matrix()
    clean = norm.clean_tpm(vals, gt)
    rem = norm._censored_mask(gt, exclude_ribosomal_proteins=True).to_numpy()
    # technical -> 250k, biological -> 750k, exactly, per sample.
    assert np.allclose(clean.loc[rem].sum(), 250_000.0)
    assert np.allclose(clean.loc[~rem].sum(), 750_000.0)
    # within-biology ratios preserved (300:400:500)
    bio = clean.loc[~rem, "s1"].to_numpy()
    assert np.allclose(bio / bio.min(), [1.0, 4 / 3, 5 / 3])


def test_clean_tpm_no_technical_mass():
    gt, vals = _matrix()
    rem = norm._censored_mask(gt, exclude_ribosomal_proteins=True)
    vals.loc[rem] = 0.0  # a sample/library with no technical reads
    clean = norm.clean_tpm(vals, gt)
    assert np.allclose(clean.loc[rem.to_numpy()].sum(), 0.0)  # technical stays 0
    assert np.allclose(clean.loc[~rem.to_numpy()].sum(), 750_000.0)  # biology still 75%


def test_clean_tpm_validates():
    gt, vals = _matrix()
    with pytest.raises(ValueError, match="technical_fraction"):
        norm.clean_tpm(vals, gt, technical_fraction=1.5)
    with pytest.raises(ValueError, match="gene_table or removable"):
        norm.clean_tpm(vals)


def test_drop_technical_vs_filter_technical_rna():
    gt, _ = _matrix()
    df = gt.assign(s1=1.0)
    # drop_technical_genes removes the censored set (incl. ribosomal); the two
    # mito genes here are technical, so both go.
    assert len(norm.drop_technical_genes(df)) == 3
    assert len(norm.filter_technical_rna(df)) == 3  # mito is technical RNA too


def test_log_and_rank_helpers():
    gt, vals = _matrix()
    df = vals.assign(Ensembl_Gene_ID=gt["Ensembl_Gene_ID"], Symbol=gt["Symbol"])
    log1p = norm.log1p_transform(df)
    assert np.allclose(log1p["s1"], np.log1p(vals["s1"]))
    rank = norm.percentile_rank(df)
    assert rank["s1"].max() <= 100.0 and rank["s1"].min() >= 0.0


def test_normalize_to_housekeeping():
    # housekeeping genes present -> column scaled by their median (1.0 = baseline)
    hk = list(gf.housekeeping_gene_ids())[:2]
    gt = pd.DataFrame({"Ensembl_Gene_ID": [*hk, "ENSG00000999999"], "Symbol": ["H1", "H2", "X"]})
    df = gt.assign(s1=[10.0, 30.0, 100.0])  # hk median = 20
    out = norm.normalize_to_housekeeping(df)
    assert np.allclose(out["s1"], [0.5, 1.5, 5.0])
