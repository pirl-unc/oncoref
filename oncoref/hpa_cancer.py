"""HPA cancer IHC patient counts and unpaired, native-pTPM RNA cohorts.

Gene IDs stay in the source Ensembl space. These measurements are neither
clean TPM nor matched RNA/protein observations. No antibody reliability or
fraction of stained cells can be recovered from the aggregate IHC download.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import reference_data

_DATA = Path(__file__).parent / "data"
THRESHOLDS = (0.1, 1.0, 2.0, 5.0)
_IHC_COUNTS = ["high", "medium", "low", "not_detected"]


def _threshold_suffix(threshold: float) -> str:
    if threshold not in THRESHOLDS:
        raise ValueError(f"threshold must be one of {THRESHOLDS}")
    return f"{threshold:g}".replace(".", "_")


def _select_genes(frame: pd.DataFrame, gene_ids) -> pd.DataFrame:
    if gene_ids is None:
        return frame.copy()
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]
    return frame.loc[frame["gene_id"].isin(gene_ids)].copy()


def hpa_cancer_sources() -> dict:
    """Fresh provenance snapshot: release, raw hashes, scope and legacy caveat."""
    return json.loads((_DATA / "hpa-cancer-sources.json").read_text())


def hpa_cancer_crosswalk() -> pd.DataFrame:
    """All 20 IHC groups, required RNA types, and explicit mapping limitations.

    ``rna_codes`` is pipe-delimited. ``scope_mismatch`` and ``unmatched`` rows
    are retained for audit; the comparison API does not compute RNA values for
    these groups. Anatomical correspondence does not imply matched histologies.
    """
    return pd.read_csv(_DATA / "hpa-cancer-crosswalk.csv", keep_default_na=False)


def hpa_cancer_rna_cohorts() -> pd.DataFrame:
    """Source cohort labels, codes and nominal sample counts, separated by cohort."""
    return pd.read_csv(_DATA / "hpa-cancer-rna-cohorts.csv")


def hpa_cancer_assay_limitations() -> pd.DataFrame:
    """Separate interpretation cautions; never alters the observed counts.

    Gene-specific cautions were curated for normal tissue. The cancer download
    has no antibody identifiers with which to establish that the same reagent
    was used. These are prompts to check assay evidence, not cancer exclusions.
    """
    general = pd.DataFrame(
        [
            (
                "",
                "cohort",
                "RNA and IHC cohorts are unpaired; comparisons are descriptive, not paired calibration.",
            ),
            (
                "",
                "measurement",
                "IHC prevalence is the fraction of scored patients, not the fraction of stained cells.",
            ),
            (
                "",
                "antibody",
                "Aggregate cancer counts do not identify antibodies, reliability, or cross-reactivity; gene/paralog specificity requires separate assay evidence.",
            ),
        ],
        columns=["gene_id", "scope", "caution"],
    )
    curated = pd.read_csv(_DATA / "cta-ihc-unreliable.csv")
    curated = curated.rename(
        columns={
            "Ensembl_Gene_ID": "gene_id",
            "Symbol": "symbol",
            "reason": "caution",
        }
    )
    curated["scope"] = "normal_tissue_curation_not_cancer_revalidation"
    return pd.concat([general, curated], ignore_index=True, sort=False)


def _summarize_ihc(raw: pd.DataFrame) -> pd.DataFrame:
    frame = raw.rename(
        columns={
            "Gene": "gene_id",
            "Gene name": "hpa_gene_name",
            "Cancer": "cancer",
            "High": "high",
            "Medium": "medium",
            "Low": "low",
            "Not detected": "not_detected",
        }
    ).copy()
    keys = ["gene_id", "cancer"]
    if frame[keys].isna().any().any() or frame.duplicated(keys).any():
        raise ValueError("IHC requires unique nonmissing gene/cancer identities")
    for col in _IHC_COUNTS:
        values = pd.to_numeric(frame[col], errors="raise")
        known = values.dropna()
        if (~np.isfinite(known) | (known < 0) | (known % 1 != 0)).any():
            raise ValueError("IHC counts must be nonnegative finite integers or missing")
        frame[col] = values.astype("Int64")
    # A missing category is unknown, never an implicit zero in the denominator.
    frame["total"] = frame[_IHC_COUNTS].sum(axis=1, min_count=4).astype("Int64")
    frame["detected"] = frame[["high", "medium", "low"]].sum(axis=1, min_count=3).astype("Int64")
    frame["medium_high"] = frame[["high", "medium"]].sum(axis=1, min_count=2).astype("Int64")
    denominator = frame["total"].where(frame["total"] > 0)
    frame["prevalence_detected"] = frame["detected"] / denominator
    frame["prevalence_medium_high"] = frame["medium_high"] / denominator
    frame["measurement_status"] = "measured"
    frame.loc[frame["total"].eq(0), "measurement_status"] = "no_scored_patients"
    frame.loc[frame["total"].isna(), "measurement_status"] = "missing_categories"
    return frame


def hpa_cancer_ihc_prevalence(gene_ids=None) -> pd.DataFrame:
    """Genome-wide scored-patient counts and any / medium-high IHC fractions.

    Unknown categories and zero denominators have missing fractions. A measured
    zero positive count remains zero. Gene names are the HPA source labels.
    """
    raw = pd.read_csv(reference_data.ensure("hpa_cancer_ihc"), sep="\t")
    frame = _summarize_ihc(raw)
    unknown = set(frame["cancer"]) - set(hpa_cancer_crosswalk()["cancer"])
    if unknown:
        raise ValueError(f"unmapped HPA IHC groups: {sorted(unknown)}")
    return _select_genes(frame, gene_ids)


def _summarize_rna_rows(raw: pd.DataFrame, cohorts: pd.DataFrame) -> pd.DataFrame:
    """Aggregate complete genes from raw sample rows (used by the offline builder)."""
    keys = ["Gene", "Cancer", "Sample"]
    if raw[keys].isna().any().any() or raw.duplicated(keys).any():
        raise ValueError("RNA requires unique nonmissing gene/cancer/sample identities")
    values = pd.to_numeric(raw["pTPM"], errors="raise")
    known = values.dropna()
    if (~np.isfinite(known) | (known < 0)).any():
        raise ValueError("RNA pTPM must be nonnegative and finite or missing")
    frame = raw.assign(pTPM=values)
    group = frame.groupby(["Gene", "Cancer"], sort=False, observed=True)
    summary = group.agg(
        source_samples=("Sample", "size"),
        samples=("pTPM", "count"),
        sum_ptpm=("pTPM", "sum"),
        max_ptpm=("pTPM", "max"),
    )
    summary["missing_samples"] = summary["source_samples"] - summary["samples"]
    denominator = summary["samples"].where(summary["samples"] > 0)
    summary["mean_ptpm"] = summary["sum_ptpm"] / denominator
    for threshold in THRESHOLDS:
        suffix = _threshold_suffix(threshold)
        count = f"expressed_samples_ptpm_ge_{suffix}"
        summary[count] = (
            frame["pTPM"]
            .ge(threshold)
            .groupby([frame["Gene"], frame["Cancer"]], sort=False, observed=True)
            .sum()
        ).astype("int64")
        summary[f"prevalence_ptpm_ge_{suffix}"] = summary[count] / denominator
    summary = summary.reset_index().rename(columns={"Gene": "gene_id", "Cancer": "cancer"})
    metadata = cohorts.copy()
    metadata["cancer"] = metadata["cancer_type"] + " (" + metadata["cohort"] + ")"
    # HPA has trailing whitespace in a few source labels; normalize only that.
    summary["cancer"] = summary["cancer"].str.replace(r"\s+\(", " (", regex=True)
    summary = summary.merge(metadata, on="cancer", how="left", validate="many_to_one")
    if summary["cohort"].isna().any():
        raise ValueError("unrecognized HPA RNA cancer/cohort label")
    if (summary["source_samples"] > summary["nominal_samples"]).any():
        raise ValueError("RNA sample counts exceed the source cohort size")
    return summary


def hpa_cancer_rna_prevalence(gene_ids=None, *, cohort: str | None = None) -> pd.DataFrame:
    """Genome-wide RNA counts at 0.1, 1, 2 and 5 native pTPM, plus mean/max.

    ``samples`` counts finite measurements; ``missing_samples`` counts explicit
    missing measurements. ``source_samples`` can be smaller than the cohort's
    ``nominal_samples``. Absent gene/cohort rows are not filled with negatives.
    TCGA and validation are always separate rows, including when cohort=None.
    """
    if cohort not in (None, "TCGA", "validation"):
        raise ValueError("cohort must be 'TCGA', 'validation', or None")
    frame = pd.read_csv(reference_data.ensure("hpa_cancer_rna"), sep="\t")
    if cohort is not None:
        frame = frame.loc[frame["cohort"] == cohort]
    return _select_genes(frame, gene_ids)


def _compare(ihc, rna, crosswalk, *, cohort, threshold):
    suffix = _threshold_suffix(threshold)
    count = f"expressed_samples_ptpm_ge_{suffix}"
    result = ihc.merge(crosswalk, on="cancer", how="left", validate="many_to_one")
    if result["mapping_status"].isna().any():
        raise ValueError("unmapped IHC cancer group")
    mapping = crosswalk.loc[crosswalk["mapping_status"].isin(["single", "pooled"])]
    mapping = mapping.assign(cancer_code=mapping["rna_codes"].str.split("|")).explode("cancer_code")
    frame = rna.loc[rna["cohort"] == cohort].copy()
    if frame.duplicated(["gene_id", "cancer_code"]).any():
        raise ValueError("duplicate RNA gene/type within cohort")
    frame = frame.merge(mapping[["cancer", "cancer_code"]], on="cancer_code", suffixes=("_rna", ""))
    pooled = frame.groupby(["gene_id", "cancer"], as_index=False, observed=True).agg(
        rna_samples=("samples", "sum"),
        rna_positive_samples=(count, "sum"),
        rna_sum_ptpm=("sum_ptpm", "sum"),
        rna_available_types=("cancer_code", "nunique"),
    )
    result = result.merge(pooled, on=["gene_id", "cancer"], how="left", validate="one_to_one")
    result["rna_required_types"] = (
        result["rna_codes"].map(lambda v: len(v.split("|")) if v else 0).astype(int)
    )
    result["rna_available_types"] = result["rna_available_types"].fillna(0).astype(int)
    result["comparison_status"] = "comparable"
    result.loc[result["rna_samples"].eq(0), "comparison_status"] = "no_measured_rna"
    result.loc[
        result["rna_available_types"] < result["rna_required_types"], "comparison_status"
    ] = "incomplete_rna_types"
    result.loc[result["rna_available_types"].eq(0), "comparison_status"] = "missing_rna"
    # A present type with no finite measurements cannot supply its share of a pool.
    empty = frame.loc[frame["samples"].eq(0), ["gene_id", "cancer"]].drop_duplicates()
    if not empty.empty:
        empty_index = pd.MultiIndex.from_frame(empty)
        mask = pd.MultiIndex.from_frame(result[["gene_id", "cancer"]]).isin(empty_index)
        result.loc[mask, "comparison_status"] = "no_measured_rna"
    excluded = result["mapping_status"].isin(["scope_mismatch", "unmatched"])
    result.loc[excluded, "comparison_status"] = result.loc[excluded, "mapping_status"]
    valid_rna = result["comparison_status"].eq("comparable")
    numeric = ["rna_samples", "rna_positive_samples", "rna_sum_ptpm"]
    result.loc[~valid_rna, numeric] = np.nan
    result["rna_mean_ptpm"] = result["rna_sum_ptpm"] / result["rna_samples"]
    result["rna_prevalence"] = result["rna_positive_samples"] / result["rna_samples"]
    result.loc[valid_rna & result["total"].eq(0), "comparison_status"] = "no_ihc_patients"
    result.loc[valid_rna & result["total"].isna(), "comparison_status"] = "missing_ihc_categories"
    result["cohort"] = cohort
    result["threshold_ptpm"] = threshold
    result["cohorts_paired"] = False
    return result


def hpa_cancer_rna_ihc_comparison(gene_ids=None, *, cohort="TCGA", threshold=1.0) -> pd.DataFrame:
    """Descriptive unpaired comparisons, retaining unmapped and incomplete rows.

    Colorectal, lung and renal RNA use sample-weighted pooling only when every
    required cancer type has measured RNA for that gene. Glioma/GBM is a scope
    mismatch. All fractions concern patients/samples, never stained cells.
    """
    if cohort not in ("TCGA", "validation"):
        raise ValueError("select one cohort: 'TCGA' or 'validation'")
    _threshold_suffix(threshold)
    if gene_ids is not None and not isinstance(gene_ids, str):
        gene_ids = list(gene_ids)
    return _compare(
        hpa_cancer_ihc_prevalence(gene_ids),
        hpa_cancer_rna_prevalence(gene_ids, cohort=cohort),
        hpa_cancer_crosswalk(),
        cohort=cohort,
        threshold=threshold,
    )
