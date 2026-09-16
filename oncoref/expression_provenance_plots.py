# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Provenance and coverage figures for oncoref's reference expression data.

The CTA-curation figures answer "how is the gene set defined"; these answer the
matching question for the expression side — what reference expression data does
oncoref actually ship, where did it come from, how deep is each cancer type, how
trustworthy is each source's quantification, and how much of the cancer space
(both by type count and by disease burden) it covers.

Everything is read from three packaged registries, so the figures render with no
downloads and no per-sample matrices in cache:

``cancer-reference-expression-availability``
    One row per cancer-code x candidate source, with the project, publication,
    pipeline, scale class and sample/gene counts. ``selected`` marks the source
    actually used for that code.
``cancer-type-registry``
    The full cancer-type ontology — the denominator for "what is not covered".
``cancer-incidence-mortality``
    US and world incidence/mortality shares per burden category, which turn
    per-type coverage into burden-weighted coverage.

Burden-weighted coverage counts a burden category as covered when at least one
of its cancer codes has a selected expression source; a category like
``other_and_unknown_primary`` has no codes of its own and so reads as uncovered
by construction rather than by omission.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import figure_style
from .figure_style import DROP, KEPT, WEAK
from .load_dataset import get_data

FILENAMES = {
    "source_samples": "expr-source-samples.png",
    "samples_per_code": "expr-samples-per-code.png",
    "source_quality": "expr-source-quality.png",
    "family_coverage": "expr-family-coverage.png",
    "burden_coverage": "expr-burden-coverage.png",
}

#: Scale classes whose values are directly comparable as linear TPM. Everything
#: else is a proxy (microarray intensity, CPM, pseudobulk nTPM) that is usable
#: within a cohort but not on a shared absolute axis.
_TPM_EXACT = "linear_rnaseq_tpm"

#: Display names for the source_type tags, which are pipeline-level slugs.
_SOURCE_TYPE_LABEL = {
    "treehouse-compendium": "Treehouse compendium",
    "treehouse-derived": "Treehouse (derived)",
    "geo-matrix": "GEO (series matrix)",
    "geo-rnaseq": "GEO (RNA-seq)",
    "geo-microarray": "GEO (microarray)",
    "geo-rnaseq-cpm-proxy": "GEO (CPM proxy)",
    "geo-rnaseq-rpkm": "GEO (RPKM)",
    "geo-scrna-pseudobulk": "GEO (scRNA pseudobulk)",
    "geo": "GEO (other)",
    "gdc": "GDC",
    "recount3": "recount3",
    "openpbta-rnaseq-rsem": "OpenPBTA",
    "sra-ncbi-counts": "SRA/NCBI counts",
    "github-release": "GitHub release",
    "unc-case-series": "UNC case series",
    "zenodo": "Zenodo",
}


def _bool_series(series):
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def availability(selected_only=True):
    """Reference-expression availability rows, by default only the selected source."""
    df = get_data("cancer-reference-expression-availability").copy()
    if selected_only and "selected" in df.columns:
        df = df[_bool_series(df["selected"])].reset_index(drop=True)
    return df


def _registry():
    return get_data("cancer-type-registry").copy()


def coverage_summary() -> dict:
    """Headline counts behind the coverage figures."""
    av = availability()
    reg = _registry()
    covered = set(av["cancer_code"].astype(str))
    n_codes = len(reg)
    n_covered = int(reg["code"].astype(str).isin(covered).sum())
    scale = av["source_scale_class"].astype(str)
    return {
        "n_sources": int(av["source_project"].nunique()),
        "n_codes_total": n_codes,
        "n_codes_covered": n_covered,
        "n_samples": int(av["n_reference_samples"].sum()),
        "n_tpm_exact_codes": int((scale == _TPM_EXACT).sum()),
        "median_samples_per_code": float(av["n_reference_samples"].median()),
    }


def _save(fig, path, plt):
    figure_style.save(fig, path)
    plt.close(fig)


def _fig_source_samples(path, plt):
    """Samples contributed by each upstream data source."""
    av = availability()
    grp = (
        av.groupby(av["source_type"].astype(str))["n_reference_samples"]
        .agg(["sum", "size"])
        .sort_values("sum", ascending=False)
    )
    labels = [_SOURCE_TYPE_LABEL.get(k, k) for k in grp.index]
    samples = grp["sum"].to_numpy()
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.34 * len(labels) + 1.3))
    ax.barh(y, samples, color=KEPT, height=0.68)
    for i, (n, k) in enumerate(zip(samples, grp["size"].to_numpy())):
        ax.text(n + samples.max() * 0.012, i, f"{n:,} ({k})", va="center", fontsize=8.5)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("samples")
    ax.set_xlim(0, samples.max() * 1.22)
    _save(fig, path, plt)


def _fig_samples_per_code(path, plt):
    """How deep each covered cancer type is, in samples."""
    av = availability()
    n = av["n_reference_samples"].astype(float).to_numpy()
    n = n[n > 0]
    bins = np.logspace(0, np.log10(max(n.max(), 10)) + 0.05, 26)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.hist(n, bins=bins, color=KEPT)
    ax.axvline(np.median(n), color=figure_style.THRESHOLD, ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("samples per cancer type")
    ax.set_ylabel("cancer types")
    _save(fig, path, plt)


def _fig_source_quality(path, plt):
    """Quantification quality: exact linear TPM vs proxy scales."""
    av = availability()
    scale = av["source_scale_class"].astype(str)
    grp = av.groupby(scale)["n_reference_samples"].agg(["sum", "size"]).sort_values("sum")
    labels = [s.replace("_", " ") for s in grp.index]
    samples = grp["sum"].to_numpy()
    colors = [KEPT if s == _TPM_EXACT else WEAK for s in grp.index]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.42 * len(labels) + 1.3))
    ax.barh(y, samples, color=colors, height=0.66)
    for i, (n, k) in enumerate(zip(samples, grp["size"].to_numpy())):
        ax.text(n + samples.max() * 0.012, i, f"{n:,} ({k})", va="center", fontsize=8.5)
    ax.set_yticks(y, labels)
    ax.set_xlabel("samples")
    ax.set_xlim(0, samples.max() * 1.22)
    _save(fig, path, plt)


def _fig_family_coverage(path, plt):
    """Cancer types with and without a reference expression source, by family."""
    av = availability()
    reg = _registry()
    covered = set(av["cancer_code"].astype(str))
    reg["has"] = reg["code"].astype(str).isin(covered)
    grp = (
        reg.groupby(reg["family"].astype(str))["has"]
        .agg(["sum", "size"])
        .sort_values("size", ascending=False)
    )
    grp = grp[grp["size"] >= 2]
    labels = [s.replace("-", " ") for s in grp.index]
    have = grp["sum"].to_numpy()
    missing = grp["size"].to_numpy() - have
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(labels) + 1.5))
    ax.barh(y, have, color=KEPT, label="covered", height=0.68)
    ax.barh(y, missing, left=have, color=DROP, label="no source", height=0.68)
    total = have + missing
    for i, (h, t) in enumerate(zip(have, total)):
        ax.text(t + total.max() * 0.012, i, f"{h}/{t}", va="center", fontsize=8.5)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("cancer types")
    ax.set_xlim(0, total.max() * 1.16)
    ax.legend(loc="lower right")
    _save(fig, path, plt)


def burden_coverage(region="us"):
    """``(rows, totals)`` for burden-weighted coverage in ``region`` (us|world)."""
    from .incidence import burden_category

    av = availability()
    reg = _registry()
    im = get_data("cancer-incidence-mortality").set_index("burden_category")
    covered = set(av["cancer_code"].astype(str))
    reg["bc"] = [burden_category(str(c)) for c in reg["code"]]
    reg["has"] = reg["code"].astype(str).isin(covered)
    grp = reg.dropna(subset=["bc"]).groupby("bc")["has"].agg(["sum", "size"])
    inc_col, mor_col = f"{region}_incidence_pct", f"{region}_mortality_pct"
    grp = grp.join(im[[inc_col, mor_col]]).dropna(subset=[inc_col])
    grp["covered"] = grp["sum"] > 0
    totals = {
        "incidence_total": float(im[inc_col].sum()),
        "mortality_total": float(im[mor_col].sum()),
        "incidence_covered": float(grp.loc[grp["covered"], inc_col].sum()),
        "mortality_covered": float(grp.loc[grp["covered"], mor_col].sum()),
    }
    return grp.sort_values(inc_col, ascending=False), totals


def _fig_burden_coverage(path, plt, region="us"):
    """Share of cancer burden whose type has a reference expression source."""
    grp, _totals = burden_coverage(region)
    inc_col = f"{region}_incidence_pct"
    grp = grp[grp[inc_col] > 0.3]
    labels = [s.replace("_", " ") for s in grp.index]
    vals = grp[inc_col].to_numpy()
    frac = (grp["sum"] / grp["size"]).to_numpy()
    colors = [KEPT if f == 1 else (WEAK if f > 0 else DROP) for f in frac]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.31 * len(labels) + 1.6))
    ax.barh(y, vals, color=colors, height=0.68)
    for i, (v, s, n) in enumerate(zip(vals, grp["sum"], grp["size"])):
        ax.text(v + vals.max() * 0.015, i, f"{s}/{n}", va="center", fontsize=8.5)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel(f"{region.upper()} incidence share (%)")
    ax.set_xlim(0, vals.max() * 1.16)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=KEPT),
        plt.Rectangle((0, 0), 1, 1, color=WEAK),
        plt.Rectangle((0, 0), 1, 1, color=DROP),
    ]
    ax.legend(handles, ["all types", "some types", "none"], loc="lower right")
    _save(fig, path, plt)


_BUILDERS = {
    "source_samples": lambda p, plt, region: _fig_source_samples(p, plt),
    "samples_per_code": lambda p, plt, region: _fig_samples_per_code(p, plt),
    "source_quality": lambda p, plt, region: _fig_source_quality(p, plt),
    "family_coverage": lambda p, plt, region: _fig_family_coverage(p, plt),
    "burden_coverage": lambda p, plt, region: _fig_burden_coverage(p, plt, region),
}


def render(out_dir="expression_provenance_out", region="us") -> dict:
    """Write the expression-provenance figures into ``out_dir``.

    Returns ``{"summary": {...}, "burden": {...}, "paths": {key: Path}}``.
    """
    figure_style.apply()
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, builder in _BUILDERS.items():
        path = out / FILENAMES[key]
        builder(path, plt, region)
        paths[key] = path
    _, totals = burden_coverage(region)
    return {"summary": coverage_summary(), "burden": totals, "paths": paths}
