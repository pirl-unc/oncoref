# Licensed under the Apache License, Version 2.0.
"""Cross-tissue carryover and stoichiometric multiplex index-hopping QC (#764).

Detects sample-level multi-gene source-organ spillover (Illumina patterned-flowcell
index hopping and library-prep carryover; Nieuwenhuis et al., Genome Res 2020,
doi:10.1101/gr.255380.119) in multiplexed normal bulk RNA panels without ever
masking genuine single-gene somatic expression or signals corroborated by
orthogonal single-cell / non-multiplexed assays.
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files

import numpy as np
import pandas as pd

from .load_dataset import get_data
from .target_selection import summarize_normal_tissue_tpm

MULTIPLEX_BLEED_POLICY_ID = "oncoref.multiplex-bleed-qc.2026-11-23.1"

# Legacy 2-source advisory panels retained for backward compatibility.
CARRYOVER_PANELS = {
    "P": ("Pancreas", ("PRSS1", "CLPS", "PNLIP", "CELA3A")),
    "E": ("Esophagus_Mucosa", ("KRT4", "KRT13")),
}

# Multi-sentinel ultra-abundant source-organ panels for stoichiometric bleed QC (#764).
# Every source organ supplies >=5 high-abundance lineage markers so single-gene
# ectopic expression can never satisfy the co-detection requirement (min_sentinels >= 3).
STOICHIOMETRIC_SENTINEL_PANELS: dict[str, tuple[str, tuple[str, ...]]] = {
    "testis": (
        "Testis",
        ("PRM1", "PRM2", "TNP1", "LELP1", "SMCP", "CMTM2", "CRISP2"),
    ),
    "prostate": (
        "Prostate",
        ("KLK3", "MSMB", "KLK2", "ACP3", "KLK4", "CPLX3", "NKX3-1"),
    ),
    "pancreas": (
        "Pancreas",
        ("PRSS2", "PRSS1", "REG1A", "CPA1", "CLPS", "CELA3A", "PNLIP"),
    ),
    "liver": (
        "Liver",
        ("ALB", "HP", "ORM1", "APOA1", "APOC3", "SERPINA1", "APOA2"),
    ),
    "stomach": (
        "Stomach",
        ("PGC", "LIPF", "PGA5", "PGA3", "GKN1", "CBLIF", "TFF2"),
    ),
}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def carryover_reference(mapped: pd.DataFrame) -> pd.DataFrame:
    """Native site medians from the pinned GTEx gene crosswalk, keyed by locus."""
    panels = []
    for code, (native, symbols) in CARRYOVER_PANELS.items():
        rows = mapped.loc[mapped.Description.isin(symbols), ["gene_id", "Description", native]]
        if sorted(rows.Description) != sorted(symbols) or not rows[native].gt(0).all():
            raise ValueError(f"Missing or ambiguous carryover marker for {native}")
        panels.append(
            rows.set_index("gene_id")[[native]].rename(columns={native: "scale"}).assign(code=code)
        )
    return pd.concat(panels, verify_integrity=True)


def published_carryover(values: pd.DataFrame, site: str, reference: pd.DataFrame) -> pd.DataFrame:
    """Screen a single native-TPM site against its peers. NA means unassessed."""
    result = pd.DataFrame(index=values.columns)
    for code, (native, _) in CARRYOVER_PANELS.items():
        field = "carryover_" + code
        result[field] = pd.Series(pd.NA, index=result.index, dtype="boolean")
        if site.lower().startswith(native.split("_")[0].lower()):
            continue
        scale = reference.loc[reference.code.eq(code), "scale"]
        panel = values.reindex(scale.index).dropna(axis=1)
        if panel.shape[1] < 3:
            continue
        excess = panel.sub(panel.median(axis=1), axis=0).div(scale, axis=0)
        high = panel.gt(panel.quantile(0.95, axis=1), axis=0)
        result.loc[panel.columns, field] = (high & excess.ge(0.001)).all(axis=0)
    return result


def stoichiometric_sentinel_reference(
    mapped: pd.DataFrame,
    *,
    panels: dict[str, tuple[str, tuple[str, ...]]] = STOICHIOMETRIC_SENTINEL_PANELS,
    min_source_median_tpm: float = 100.0,
) -> pd.DataFrame:
    """Build a locus-indexed sentinel reference table with native source-organ medians."""
    frames = []
    for source_key, (native_col, symbols) in panels.items():
        if native_col not in mapped.columns:
            raise ValueError(f"Source column {native_col!r} missing from GTEx median table")
        rows = mapped.loc[
            mapped["Description"].isin(symbols),
            ["gene_id", "Description", native_col],
        ].drop_duplicates("gene_id")
        rows = rows[pd.to_numeric(rows[native_col], errors="coerce").ge(min_source_median_tpm)]
        if len(rows) < 3:
            raise ValueError(
                f"Sentinel panel {source_key!r} ({native_col}) requires >=3 high-abundance loci"
            )
        frames.append(
            rows.rename(columns={"Description": "symbol", native_col: "source_median_tpm"}).assign(
                source_key=source_key,
                source_organ=native_col,
            )
        )
    return pd.concat(frames, ignore_index=True)


def detect_stoichiometric_multiplex_bleed(
    sample_tpm: pd.DataFrame,
    *,
    site: str,
    sentinel_reference: pd.DataFrame,
    min_sentinels: int = 3,
    min_marker_tpm: float = 0.25,
    min_bleed_fraction: float = 1e-5,
    max_bleed_fraction: float = 5e-3,
    max_log10_iqr: float = 0.5,
) -> pd.DataFrame:
    """Detect biopsies carrying a proportional multi-sentinel bleed vector from another organ.

    Parameters
    ----------
    sample_tpm:
        Gene-by-sample native TPM matrix (`index` = canonical `gene_id`, `columns` = `sample_id`).
    site:
        Native tissue/site label of ``sample_tpm`` (self-organ comparisons are skipped).
    sentinel_reference:
        Output of :func:`stoichiometric_sentinel_reference`.
    min_sentinels:
        Minimum distinct sentinel genes from the same source organ co-detected at
        ``>= min_marker_tpm`` in the same biopsy (default 3; blocks single-gene ectopic signals).
    min_bleed_fraction, max_bleed_fraction:
        Accepted multiplex index-hopping / carryover band on ``TPM / source_median_TPM``.
    max_log10_iqr:
        Maximum interquartile range of ``log10(TPM / source_median_TPM)`` across detected
        sentinels in the biopsy, enforcing proportional bulk-transcriptome stoichiometry.
    """
    if min_sentinels < 3:
        raise ValueError("min_sentinels must be >= 3 to prevent single-gene masking")
    cols = [
        "sample_id",
        "site",
        "source_key",
        "source_organ",
        "sentinels_detected",
        "sentinel_symbols",
        "median_bleed_fraction",
        "log10_iqr",
    ]
    records: list[dict] = []
    site_prefix = site.lower().split("_")[0]
    for source_key, grp in sentinel_reference.groupby("source_key", sort=True):
        source_organ = str(grp["source_organ"].iloc[0])
        if site_prefix == source_organ.lower().split("_")[0]:
            continue
        # Prostate and testis are anatomical urogenital neighbors only via ductal proximity,
        # but prostate biopsies in GTEx also show flowcell testis bleed; keep evaluated unless
        # orthogonal assays show true expression.
        present = grp[grp["gene_id"].isin(sample_tpm.index)].set_index("gene_id")
        if len(present) < min_sentinels:
            continue
        sub = sample_tpm.loc[present.index]
        medians = present["source_median_tpm"].astype(float)
        syms = present["symbol"].astype(str)
        active_mask = sub.ge(min_marker_tpm)
        counts = active_mask.sum(axis=0)
        candidate_samples = counts[counts.ge(min_sentinels)].index
        for sid in candidate_samples:
            det_genes = active_mask.index[active_mask[sid]]
            ratios = (sub.loc[det_genes, sid].astype(float) / medians.loc[det_genes]).to_numpy()
            med_ratio = float(np.median(ratios))
            if not (min_bleed_fraction <= med_ratio <= max_bleed_fraction):
                continue
            log_r = np.log10(np.maximum(ratios, 1e-12))
            q75, q25 = np.percentile(log_r, [75.0, 25.0])
            iqr = float(q75 - q25)
            if iqr > max_log10_iqr:
                continue
            records.append(
                {
                    "sample_id": str(sid),
                    "site": site,
                    "source_key": str(source_key),
                    "source_organ": source_organ,
                    "sentinels_detected": len(det_genes),
                    "sentinel_symbols": ";".join(sorted(syms.loc[det_genes].tolist())),
                    "median_bleed_fraction": med_ratio,
                    "log10_iqr": iqr,
                }
            )
    return pd.DataFrame(records, columns=cols)


def compute_bleed_cleaned_donor_summary(
    values: pd.Series | pd.DataFrame,
    donors: pd.Series | pd.Index,
    *,
    flagged_samples: set[str],
    orthogonal_supported: bool = False,
    method: str = "p95",
    quantile: float = 0.95,
    n_std: float = 2.0,
) -> float:
    """Compute within-tissue donor summary after masking verified multiplex-bleed samples.

    If ``orthogonal_supported`` is True (e.g. target-organ single-cell RNA >= 0.5 nTPM,
    independent non-GTEx bulk >= 0.5 nTPM, or sorted immune RNA >= 0.5 TPM), **no** samples
    are masked even if ``flagged_samples`` is non-empty.
    """
    s = values.copy()
    if isinstance(s, pd.DataFrame):
        if s.shape[0] != 1:
            raise ValueError("Expected a 1-D sample Series or single-row DataFrame")
        s = s.iloc[0]
    s = pd.to_numeric(s, errors="coerce")
    donor_idx = pd.Series(list(donors), index=s.index, dtype="object")
    if not orthogonal_supported and flagged_samples:
        keep_mask = ~s.index.isin(flagged_samples)
        # Never mask an entire tissue to empty; require >=10 remaining donors.
        rem_donors = donor_idx[keep_mask & s.notna()].nunique()
        if rem_donors >= 10:
            s = s.where(keep_mask)
    donor_medians = s.groupby(donor_idx, sort=False).median().dropna()
    if donor_medians.empty:
        return float("nan")
    return float(
        summarize_normal_tissue_tpm(
            donor_medians,
            method=method,
            quantile=quantile,
            n_std=n_std,
        )
    )


def _verify_csv_receipt(stem: str) -> None:
    data_root = files("oncoref").joinpath("data")
    csv_bytes = data_root.joinpath(f"{stem}.csv").read_bytes()
    receipt = json.loads(data_root.joinpath(f"{stem}.json").read_text(encoding="utf-8"))
    if _sha256_bytes(csv_bytes) != receipt["sha256"]:
        raise ValueError(f"{stem}.csv does not match its SHA-256 receipt")


def multiplex_bleed_tissue_summary() -> pd.DataFrame:
    """Return SHA-256-verified per-tissue GTEx v11 stoichiometric bleed incidence (#764)."""
    _verify_csv_receipt("normal-multiplex-bleed-tissue-summary")
    return get_data("normal-multiplex-bleed-tissue-summary").copy()


def multiplex_bleed_gene_deltas(*, only_applied: bool = False) -> pd.DataFrame:
    """Return SHA-256-verified per-protein normal-baseline deltas under #764 bleed QC."""
    _verify_csv_receipt("normal-multiplex-bleed-gene-deltas")
    df = get_data("normal-multiplex-bleed-gene-deltas").copy()
    if only_applied:
        df = df[df["bleed_qc_applied"].astype(bool)].reset_index(drop=True)
    return df


def multiplex_bleed_receipt() -> dict:
    """Return the provenance and summary receipt for #764 multiplex bleed QC."""
    _verify_csv_receipt("normal-multiplex-bleed-gene-deltas")
    return json.loads(
        files("oncoref")
        .joinpath("data/normal-multiplex-bleed-gene-deltas.json")
        .read_text(encoding="utf-8")
    )


__all__ = [
    "CARRYOVER_PANELS",
    "MULTIPLEX_BLEED_POLICY_ID",
    "STOICHIOMETRIC_SENTINEL_PANELS",
    "carryover_reference",
    "compute_bleed_cleaned_donor_summary",
    "detect_stoichiometric_multiplex_bleed",
    "multiplex_bleed_gene_deltas",
    "multiplex_bleed_receipt",
    "multiplex_bleed_tissue_summary",
    "published_carryover",
    "stoichiometric_sentinel_reference",
]
