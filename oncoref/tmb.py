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

"""Tumor mutational burden (TMB) reference data, keyed by cancer-type code."""

from __future__ import annotations

from functools import lru_cache

import pandas as pd

from ._evidence_resolution import evidence_record, resolve_evidence
from .cancer_types import (
    _registry_parent_by_code,
    cancer_evidence_source_code,
    cancer_type_registry,
    resolve_cancer_type,
)
from .load_dataset import _register_derived_cache, get_data

_TMB_EVIDENCE_OVERRIDES = {
    # MSI-H/dMMR colorectal estimates are published at the CRC source scope; anatomical
    # children such as COAD_MSI/READ_MSI resolve through this row rather than duplicating
    # the same source estimate.
    "CRC_MSI": {
        "source_scope": "aggregate_source",
    },
    # The full-text WGS study explicitly separates pancreatic and midgut NET.
    "NET_MIDGUT": {
        "source_scope": "advanced_site_specific_cohort",
    },
    "NET_PANCREAS": {
        "source_scope": "advanced_site_specific_cohort",
    },
    # The curated stomach median is the pooled intestinal-type panel value (5.0 mut/Mb,
    # Chalmers 2017 Table 1). TCGA-STAD analysed its 215 tumours below 11.4 mut/Mb,
    # "none of which were MSI-positive", separately from 74 hypermutated tumours, so
    # that pooled median demonstrably does not cover MSI-H disease. Neither source
    # publishes an MSI-stratified gastric median, so this stays an explicit audited gap
    # rather than an invented estimate.
    "STAD_MSI": {
        "estimate_type": "unknown",
        "source_scope": "source_rejected_for_subtype_value",
        "missing_reason": "no_supported_subtype_median",
    },
    # These broad entities have subtype/site evidence but no defensible median
    # over the complete aggregate represented by the ontology code.
    "RCC": {
        "estimate_type": "unknown",
        "source_scope": "subtype_sources_not_aggregated",
        "missing_reason": "no_supported_aggregate_median",
    },
    "THYM_EPITHELIAL": {
        "estimate_type": "unknown",
        "source_scope": "subtype_sources_not_aggregated",
        "missing_reason": "source_reports_subtype_medians_only",
    },
    "NEN": {
        "estimate_type": "unknown",
        "source_scope": "source_rejected_for_aggregate_scope",
        "missing_reason": "advanced_subcohorts_do_not_establish_full_aggregate_median",
    },
    "NET": {
        "estimate_type": "unknown",
        "source_scope": "source_rejected_for_aggregate_scope",
        "missing_reason": "advanced_subcohorts_do_not_establish_full_aggregate_median",
    },
    "NEC": {
        "estimate_type": "unknown",
        "source_scope": "source_rejected_for_aggregate_scope",
        "missing_reason": "advanced_subcohorts_do_not_establish_full_aggregate_median",
    },
    "NEC_LUNG": {
        "estimate_type": "unknown",
        "source_scope": "subtype_sources_not_aggregated",
        "missing_reason": "no_supported_aggregate_median",
    },
    "MPN": {"source_scope": "bone_marrow_mpn_cohort"},
    "CML": {
        "estimate_type": "unknown",
        "source_scope": "no_direct_source",
        "missing_reason": "no_published_per_mb_median_curated",
    },
    "SARC_RMS_ERMS": {"estimate_type": "approximate_literature"},
    "SARC_RMS_ARMS": {"estimate_type": "approximate_literature"},
    "WILMS": {"estimate_type": "approximate_literature"},
    "BL": {"estimate_type": "approximate_literature"},
    "T_ALL": {"estimate_type": "approximate_literature"},
    "CRANIO": {"source_scope": "adamantinomatous_and_papillary_discovery_cohort"},
    "HCL": {
        "estimate_type": "approximate_capture_normalized",
        "source_scope": "classic_hcl_five_exome_cohort",
    },
    "LUAD_EGFR": {
        "estimate_type": "broader_cohort_proxy",
        "source_scope": "egfr_mutant_lung_proxy_for_luad",
    },
    "CTCL": {
        "estimate_type": "subtype_proxy",
        "source_scope": "mycosis_fungoides_proxy_for_ctcl",
    },
    "HL": {
        "estimate_type": "reported_summary",
        "source_scope": "newly_diagnosed_classical_hodgkin_sorted_hrs",
    },
    "BRCA_Normal": {
        "estimate_type": "sample_recomputed",
        "source_scope": "tcga_pancancer_normal_like_expression_subtype",
    },
    "BRCA_TNBC": {
        "estimate_type": "sample_recomputed",
        "source_scope": "tcga_primary_receptor_defined_tnbc_rna_samples",
    },
    "UCEC_CNL": {
        "estimate_type": "reported_summary",
        "source_scope": "tcga_copy_number_low_molecular_class",
    },
    "UCEC_CNH": {
        "estimate_type": "reported_summary",
        "source_scope": "tcga_copy_number_high_molecular_class",
    },
    # TCGA-SARC is a soft-tissue sarcoma cohort. It does not span the full oncoref
    # SARC member-union scope, which also includes bone sarcomas and RMS.
    "SARC": {"source_scope": "soft_tissue_sarcoma_subset"},
    "SARC_CIC": {"source_scope": "renal_cic_rearranged_sarcoma_proxy"},
    "SARC_CHON": {"source_scope": "bone_chondrosarcoma_cohort"},
    "STAD": {"source_scope": "intestinal_type_gastric_subset"},
    "CHOL": {"source_scope": "liver_cholangiocarcinoma_cohort"},
    "MESO": {"source_scope": "pleural_mesothelioma_cohort"},
    "ADCC": {"source_scope": "salivary_adenoid_cystic_cohort"},
    "UCEC_POLE": {"source_scope": "pathogenic_pole_including_multiple_classifiers"},
    "UCEC_MSI": {"source_scope": "msi_high_pole_wild_type_cohort"},
    "NBL": {"source_scope": "high_risk_neuroblastoma_cohort"},
    "NBL_MYCNamp": {
        "estimate_type": "sample_recomputed",
        "source_scope": "high_risk_mycn_amplified_cohort",
    },
    "NBL_MYCNnonamp": {"source_scope": "east_asian_mycn_nonamplified_all_risk_cohort"},
}


@lru_cache(maxsize=1)
def _aggregate_tmb_source_codes() -> frozenset[str]:
    """Registry codes whose evidence represents a mixed source-scope cohort."""
    registry = cancer_type_registry()
    mixture_cohort = registry["mixture_cohort"].fillna(False).astype(bool)
    ontology_kind = registry["ontology_kind"].fillna("").astype(str)
    source_scope = ontology_kind.str.endswith("source_scope")
    return frozenset(registry.loc[mixture_cohort & source_scope, "code"].astype(str))


_register_derived_cache(_aggregate_tmb_source_codes.cache_clear)


@lru_cache(maxsize=1)
def _checked_tmb_codes() -> frozenset[str]:
    audit = get_data("cancer-tmb-source-audit")
    return frozenset(audit.loc[audit["source_review_status"] == "source_checked", "cancer_code"])


_register_derived_cache(_checked_tmb_codes.cache_clear)


def cancer_tmb_df():
    """Return curated TMB estimates (mut/Mb) with source-review provenance.

    ``tmb_mut_mb`` selects the median, then mean, then explicitly typed estimate;
    ``tmb_statistic`` identifies that choice. The original statistic-specific
    columns remain separate. ``unspecified`` means a reported summary whose
    statistic was not established; ``approximate`` marks a derived approximation.
    Cohorts with no estimate have blank values
    (and a ``confidence`` of ``none``) so the gap is
    explicit rather than silently absent. Retained estimates span WES
    (Lawrence 2013), panels (Chalmers 2017), genome-wide WGS and disease-specific
    studies; see the ``source``/``notes`` columns — panel and WES TMB are not
    strictly comparable in the low-TMB range. ``source_review_status`` and
    ``source_locator`` distinguish rechecked numbers from legacy estimates awaiting
    source review. A citation alone does not imply a published statistic."""
    return _tmb_evidence_frame().copy()


@lru_cache(maxsize=1)
def _tmb_evidence_frame():
    """Cached annotated TMB frame. Internal callers treat it as read-only."""
    df = get_data("cancer-tmb").copy()
    estimates = df["estimate_tmb_mut_mb"].notna()
    if not df.loc[estimates, "estimate_statistic"].isin({"unspecified", "approximate"}).all():
        raise ValueError("Explicit TMB estimates require an unspecified or approximate statistic")
    df["tmb_mut_mb"] = (
        df["median_tmb_mut_mb"]
        .combine_first(df["mean_tmb_mut_mb"])
        .combine_first(df["estimate_tmb_mut_mb"])
    )
    df["tmb_statistic"] = [
        "median"
        if pd.notna(row.median_tmb_mut_mb)
        else "mean"
        if pd.notna(row.mean_tmb_mut_mb)
        else row.estimate_statistic
        if pd.notna(row.estimate_tmb_mut_mb)
        else None
        for row in df.itertuples()
    ]
    evidence = [
        tmb_evidence_fields(
            row["cancer_code"],
            row["tmb_mut_mb"],
            statistic=row["tmb_statistic"] if pd.notna(row["tmb_statistic"]) else "median",
        )
        for _, row in df.iterrows()
    ]
    for col in ("estimate_type", "source_scope", "missing_reason"):
        df[col] = [record[col] for record in evidence]
    return df.merge(
        get_data("cancer-tmb-source-audit"), on="cancer_code", how="left", validate="one_to_one"
    )


_register_derived_cache(_tmb_evidence_frame.cache_clear)


def tmb_evidence_fields(
    cancer_type: str,
    median_tmb_mut_mb: float | None,
    *,
    statistic: str = "median",
) -> dict[str, object]:
    """Classify the provenance of one explicit TMB estimate.

    ``cancer_type`` accepts a canonical registry code, display name, or alias.
    ``median_tmb_mut_mb`` is the row's numeric estimate, or ``None``/``NaN``
    when no defensible estimate is available. The parameter name is retained for
    compatibility; ``statistic`` can be ``median``, ``mean``, ``unspecified``
    (reported summary), or ``approximate`` (derived estimate). The mapping contains
    ``estimate_type``, ``source_scope``, and ``missing_reason``.

    This helper classifies the supplied row; it does not search parent or
    source-scope rows. Use :func:`resolve_tmb_source` when inherited evidence
    should be resolved. Numeric estimates for registry entities designated as
    mixed source-scope cohorts are reported as ``source_scope="aggregate_source"``.
    Reviewed per-code overrides take precedence over that registry default.
    """
    statistic_types = {
        "median": "published_median",
        "mean": "published_mean",
        "unspecified": "reported_summary",
        "approximate": "approximate_derived",
    }
    if statistic not in statistic_types:
        raise ValueError("statistic must be 'median', 'mean', 'unspecified', or 'approximate'")
    code = resolve_cancer_type(cancer_type)
    override = _TMB_EVIDENCE_OVERRIDES.get(code, {})
    if pd.isna(median_tmb_mut_mb):
        return {
            "estimate_type": override.get("estimate_type", "unknown"),
            "source_scope": override.get("source_scope", "no_direct_source"),
            "missing_reason": override.get("missing_reason", "no_published_per_mb_median_curated"),
        }
    default_scope = (
        "aggregate_source" if code in _aggregate_tmb_source_codes() else "cancer_code_direct"
    )
    estimate_type = override.get(
        "estimate_type",
        statistic_types[statistic] if code in _checked_tmb_codes() else "curated_estimate",
    )
    if estimate_type == "sample_recomputed" and statistic in {"median", "mean"}:
        estimate_type += f"_{statistic}"
    return {
        "estimate_type": estimate_type,
        "source_scope": override.get("source_scope", default_scope),
        "missing_reason": override.get("missing_reason", float("nan")),
    }


@lru_cache(maxsize=1)
def _tmb_value_map() -> dict[str, float]:
    """Cached direct numeric map. Callers must treat it as read-only."""
    vals = _tmb_evidence_frame().dropna(subset=["tmb_mut_mb"])
    return dict(zip(vals["cancer_code"].astype(str), vals["tmb_mut_mb"].astype(float)))


_register_derived_cache(_tmb_value_map.cache_clear)


@lru_cache(maxsize=1)
def _tmb_rows_by_code() -> dict[str, object]:
    """Cached direct-row index, including audited gaps. Treat rows as read-only."""
    return {str(row["cancer_code"]): row for _, row in _tmb_evidence_frame().iterrows()}


_register_derived_cache(_tmb_rows_by_code.cache_clear)


def _record_from_row(row, *, requested_code: str, resolved_code: str, inheritance_kind: str):
    record = evidence_record(
        row,
        requested_code=requested_code,
        resolved_code=resolved_code,
        inheritance_kind=inheritance_kind,
    )
    record["has_tmb_source"] = True
    return record


def _resolve_tmb_row(requested_code: str, *, inherit: bool):
    values = _tmb_value_map()
    rows = _tmb_rows_by_code()
    resolution = resolve_evidence(
        requested_code,
        direct_lookup=lambda code: rows.get(code) if code in values else None,
        direct_gap_lookup=rows.get,
        source_code_for=cancer_evidence_source_code,
        parent_by_code=_registry_parent_by_code,
        inherit=inherit,
    )
    return resolution.resolved_code, resolution.inheritance_kind, resolution.payload


def resolve_tmb_source(cancer_type, *, inherit=True) -> dict:
    """Resolve the source row used for a TMB lookup.

    Returns metadata without requiring callers to inspect the raw table:

    - ``requested_cancer_code``: canonical code requested by the caller.
    - ``resolved_cancer_code``: direct/source-scope/ancestor row used, if any.
    - ``inheritance_kind``: ``"direct"``, ``"direct_missing"``, ``"source_scope"``,
      ``"ancestor"``, or ``"missing"``.
    - source/provenance fields from the selected row when available.

    This makes aggregate evidence explicit. For example, ``COAD_MSI`` and
    ``READ_MSI`` resolve through the curated ``CRC_MSI`` TMB row, preserving that the
    source estimate is CRC-level MSI-H/dMMR evidence rather than a colon- or
    rectum-specific median.
    """
    requested_code = resolve_cancer_type(cancer_type)
    resolved_code, inheritance_kind, row = _resolve_tmb_row(requested_code, inherit=inherit)
    if row is None:
        return {
            "requested_cancer_code": requested_code,
            "resolved_cancer_code": None,
            "inheritance_kind": inheritance_kind,
            "is_inherited_evidence": False,
            "has_tmb_source": False,
        }
    return _record_from_row(
        row,
        requested_code=requested_code,
        resolved_code=resolved_code,
        inheritance_kind=inheritance_kind,
    )


def cancer_tmb(cancer_type=None, *, inherit=True):
    """TMB (mut/Mb), preferring median, then mean, then an explicitly typed estimate.

    With no cancer type, return the whole ``{code: tmb}`` map, omitting gaps.
    Use :func:`cancer_tmb_record` to inspect the statistic, assay and population.

    ``cancer_type`` is resolved through :func:`resolve_cancer_type`, so aliases
    and display names work. When ``inherit`` (default), a code with no curated
    value of its own inherits its nearest ancestor's TMB by walking the registry
    ``parent_code`` chain — so molecular / histology subtypes (``SCLC_ASCL1`` ->
    ``SCLC``, rare ``SARC_*`` -> ``SARC``) resolve
    without a curated row each. Explicit audited gaps block inheritance. Returns
    ``None`` if no eligible source has a value."""
    mapping = _tmb_value_map()
    if cancer_type is None:
        return dict(mapping)
    code = resolve_cancer_type(cancer_type)
    resolved_code, _, row = _resolve_tmb_row(code, inherit=inherit)
    return mapping.get(resolved_code) if row is not None else None


def cancer_tmb_record(cancer_type=None, *, inherit=True):
    """Metadata-bearing TMB lookup.

    Mirrors :func:`cancer_tmb`, but returns the resolved source row as a dict instead
    of only the numeric estimate. The record includes ``tmb_mut_mb``,
    ``tmb_statistic``, and the derived evidence columns from
    :func:`cancer_tmb_df` plus lookup metadata:

    - ``requested_cancer_code``
    - ``resolved_cancer_code``
    - ``inheritance_kind``
    - ``is_inherited_evidence``

    With ``cancer_type=None`` the returned bulk map contains direct source rows only.
    Use :func:`resolve_tmb_source` for explicit requested-code resolution metadata.
    """
    if cancer_type is None:
        return {
            code: record
            for code in sorted(cancer_tmb())
            if (record := cancer_tmb_record(code, inherit=False))
        }
    record = resolve_tmb_source(cancer_type, inherit=inherit)
    return record if record.get("has_tmb_source") else None


def cancer_frameshift_burden_df():
    """Return the curated ``cancer-frameshift-burden.csv`` reference: per-type
    frameshift-indel burden (``cancer_code``, ``indel_class``, ``indel_score``,
    ``basis``, ``pmid_doi``, ``confidence``, ``notes``). A complement to TMB —
    frameshift indels yield disproportionately many high-affinity neoantigens."""
    return get_data("cancer-frameshift-burden")


def cancer_frameshift_burden(cancer_type=None, *, inherit=True):
    """Frameshift-indel ``indel_score`` for a cancer type (alias-resolved), or the
    full ``{code: score}`` map when ``cancer_type`` is None. Subtypes inherit a
    parent's score if unmapped (mirrors :func:`cancer_tmb`)."""
    df = cancer_frameshift_burden_df()
    scores = {
        str(c): int(s)
        for c, s in zip(df["cancer_code"], df["indel_score"])
        if str(s).strip() and str(s).strip().lower() != "nan"
    }
    if cancer_type is None:
        return scores
    code = resolve_cancer_type(cancer_type)
    if code in scores:
        return scores[code]
    if inherit:
        parent_of = cancer_type_registry().set_index("code")["parent_code"].to_dict()
        cur = code
        seen = set()
        while cur in parent_of and isinstance(parent_of[cur], str) and parent_of[cur]:
            cur = parent_of[cur]
            if cur in seen:
                break
            seen.add(cur)
            if cur in scores:
                return scores[cur]
    return None
