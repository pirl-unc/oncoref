# Licensed under the Apache License, Version 2.0.
"""Version- and provider-agnostic GTEx RNA loaders across ``(qc, cleanup, normalize)`` (#778)."""

from __future__ import annotations

import functools
import gzip
import warnings
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from .. import bulk_qc, normal_tissues, normalization, reference_data
from . import deconvolution as gtex_deconvolution
from .spec import (
    _HPA_TISSUE_TO_GTEX_SITE,
    _UNSET,
    CELL_LINE_COLUMNS,
    GTEX_V11_SUBCOMPARTMENT_COLUMNS,
    HPA_BRAIN_GTEX_SOURCE_ID,
    HPA_TISSUE_GTEX_SOURCE_ID,
    HPA_UNIQUE_GTEX_TISSUES,
    SHARD_SOURCE_ID,
    SOURCE_ID,
    GtexSpec,
    _effective_derived_cache_dir,
    normalize_column_to_site,
    primary_shard_sites,
    resolve_spec,
    spec_cache_key,
)


def versions(provider: str = "direct", *, granularity: str = "medians") -> tuple[str, ...]:
    """Registered concrete GTEx versions for ``provider`` (``'direct'`` or ``'hpa'``)."""
    prov = provider.strip().lower()
    gran = granularity.strip().lower()
    if prov == "direct":
        if gran == "medians":
            return reference_data.versions(SOURCE_ID)
        if gran == "samples":
            pins = reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"]
            rels = sorted({pk.split("-", 1)[0] for pk in pins}, reverse=True)
            return tuple(rels)
        raise ValueError(f"Unknown granularity {granularity!r}; expected 'medians' or 'samples'")
    if prov == "hpa":
        return reference_data.versions(HPA_TISSUE_GTEX_SOURCE_ID)
    raise ValueError(f"Unknown GTEx provider {provider!r}; expected 'direct' or 'hpa'")


def _is_boundary_columns(src_col: str, tgt_col: str) -> bool:
    s_site = normalize_column_to_site(src_col)
    t_site = normalize_column_to_site(tgt_col)
    if gtex_deconvolution.is_anatomical_boundary(s_site, t_site):
        return True
    return bool(
        str(tgt_col).lower().endswith("lymphoid_aggregate")
        and gtex_deconvolution.organ_lineage_family(s_site) == "lymphoid_hematopoietic"
    )


def apply_summary_spec(
    df: pd.DataFrame,
    *,
    spec: GtexSpec,
    provider: str = "direct",
    value_column: str = "TPM",
) -> pd.DataFrame:
    """Pass any GTEx summary matrix (direct ``v8``/``v10``/``v11`` or HPA ``v23``) through ``spec``."""
    if spec.sex is not None:
        raise ValueError(
            "Donor sex filtering (sex='XX' or 'XY') requires per-sample GTEx shards "
            "(provider='direct', version='v11')"
        )
    out = df.copy()
    id_cols = [c for c in ("Name", "Description") if c in out.columns]
    val_cols = [c for c in out.columns if c not in id_cols]

    if spec.drop_cell_lines:
        val_cols = [c for c in val_cols if c not in CELL_LINE_COLUMNS]
    if not spec.include_subcompartments:
        val_cols = [c for c in val_cols if c not in GTEX_V11_SUBCOMPARTMENT_COLUMNS]

    out = out[id_cols + val_cols].copy()
    if not val_cols:
        return out

    mat = out[val_cols].to_numpy(dtype=np.float64, copy=True)
    orig_nan = np.isnan(mat)

    clean_summary_splash = spec.reject_splash or spec.nnls_deconvolution
    if spec.subtract_index_hopping or clean_summary_splash:
        anchors = gtex_deconvolution.gtex_deconvolution_anchors()
        unver = out["Name"].astype(str).str.split(".").str[0]
        first_pos: dict[str, int] = {}
        for pos, gid in enumerate(unver):
            if gid not in first_pos:
                first_pos[gid] = pos

        mat_clean = np.fmax(np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0), 0.0)
        col_sites = [normalize_column_to_site(c) for c in val_cols]
        site_to_col_idxs: dict[str, list[int]] = {}
        for j, st in enumerate(col_sites):
            site_to_col_idxs.setdefault(st, []).append(j)

        nnls_w_mat: pd.DataFrame | None = None
        if spec.nnls_deconvolution:
            anchor_present = [g for g in anchors["gene_id"].unique() if g in first_pos]
            if len(anchor_present) >= 6:
                a_pos = [first_pos[g] for g in anchor_present]
                sub_df = pd.DataFrame(mat_clean[a_pos, :], index=anchor_present, columns=val_cols)
                basis_dict = {
                    st: sub_df[[val_cols[i] for i in idxs]].max(axis=1)
                    for st, idxs in site_to_col_idxs.items()
                }
                nnls_w_mat = bulk_qc.solve_nonnegative_archetype_mixture(
                    pd.DataFrame(basis_dict, index=anchor_present),
                    sub_df,
                    l1_penalty=5e-4,
                )

        best_by_col_fam: dict[tuple[int, str], tuple[float, float, np.ndarray]] = {}
        for src_site, grp in anchors.groupby("source_site", sort=True):
            if src_site not in site_to_col_idxs:
                continue
            gids = [g for g in grp["gene_id"] if g in first_pos]
            if len(gids) < 3:
                continue
            src_col_idxs = site_to_col_idxs[src_site]
            src_profile = mat_clean[:, src_col_idxs].max(axis=1)
            a_rows = np.array([first_pos[g] for g in gids], dtype=int)
            src_anchor_vals = src_profile[a_rows]
            valid_a = src_anchor_vals >= 20.0
            if int(valid_a.sum()) < 3:
                continue
            a_rows = a_rows[valid_a]
            src_anchor_vals = src_anchor_vals[valid_a]
            excl_sets = [
                {x for x in str(ex).split(";") if x}
                for ex, keep in zip(grp["excluded_target_sites"], valid_a)
                if keep
            ]
            src_fam = gtex_deconvolution.organ_lineage_family(str(src_site))

            for j, tgt_col in enumerate(val_cols):
                tgt_site = col_sites[j]
                if _is_boundary_columns(src_site, tgt_col):
                    continue
                keep_t = np.array([tgt_site not in ex for ex in excl_sets], dtype=bool)
                if int(keep_t.sum()) < 3:
                    continue
                use_rows = a_rows[keep_t]
                use_src = src_anchor_vals[keep_t]
                obs_a = mat[use_rows, j]
                det = np.isfinite(obs_a) & (obs_a >= 0.05)
                if int(det.sum()) < 3:
                    continue
                ratios = obs_a[det] / use_src[det]
                med_r = float(np.median(ratios))
                if med_r < 2e-5 or med_r > 0.25:
                    continue
                log_r = np.log10(np.maximum(ratios, 1e-12))
                q75, q25 = np.percentile(log_r, [75.0, 25.0])
                if float(q75 - q25) > 0.65:
                    continue
                if (
                    spec.nnls_deconvolution
                    and nnls_w_mat is not None
                    and med_r >= 5e-3
                    and src_site in nnls_w_mat.index
                    and float(nnls_w_mat.at[src_site, tgt_col]) < 0.15 * med_r
                ):
                    continue

                is_hop = med_r <= 5e-3
                if (is_hop and not spec.subtract_index_hopping) or (
                    not is_hop and not clean_summary_splash
                ):
                    continue

                coeff = float(np.percentile(ratios, 25.0)) if is_hop else med_r
                prev = best_by_col_fam.get((j, src_fam))
                if prev is None or med_r > prev[0]:
                    best_by_col_fam[(j, src_fam)] = (med_r, coeff, src_profile)

        expected_bleed = np.zeros_like(mat_clean)
        for (j, _sfam), (_med_r, coeff, src_profile) in best_by_col_fam.items():
            expected_bleed[:, j] += coeff * src_profile

        substantial = (
            ~orig_nan
            & (mat_clean > 0.0)
            & (expected_bleed > 0.0)
            & (expected_bleed >= 0.20 * mat_clean)
        )
        mat = np.where(substantial, np.fmax(0.0, mat_clean - expected_bleed), mat)
        mat[orig_nan] = np.nan
        out[val_cols] = mat

    if provider.strip().lower() == "hpa" and spec.hpa_unique_tissues_only:
        keep_unique = [c for c in val_cols if c.lower() in HPA_UNIQUE_GTEX_TISSUES]
        out = out[id_cols + keep_unique].copy()
        val_cols = keep_unique
        orig_nan = (
            np.isnan(out[val_cols].to_numpy(dtype=np.float64))
            if val_cols
            else np.zeros((0, 0), dtype=bool)
        )

    if spec.clean_tpm and val_cols:
        if value_column.upper() != "TPM":
            raise ValueError("normalize='clean_tpm' requires linear 'TPM', not 'nTPM'")
        gene_tbl = pd.DataFrame({"Ensembl_Gene_ID": out["Name"].astype(str).str.split(".").str[0]})
        filled = out[val_cols].fillna(0.0)
        norm_arr = normalization.clean_tpm(
            filled, gene_tbl, on_empty_compartment="record"
        ).to_numpy(dtype=np.float64, copy=True)
        norm_arr[orig_nan] = np.nan
        out[val_cols] = norm_arr

    return out


def _load_hpa_gtex_wide(version: str, *, value_column: str = "TPM") -> tuple[pd.DataFrame, dict]:
    if value_column not in {"TPM", "nTPM", "pTPM"}:
        raise ValueError(f"value_column must be 'TPM', 'nTPM', or 'pTPM'; got {value_column!r}")
    v_t = reference_data.resolve_version(HPA_TISSUE_GTEX_SOURCE_ID, version)
    v_b = reference_data.resolve_version(HPA_BRAIN_GTEX_SOURCE_ID, version)
    p_t = reference_data.ensure(HPA_TISSUE_GTEX_SOURCE_ID, v_t)
    p_b = reference_data.ensure(HPA_BRAIN_GTEX_SOURCE_ID, v_b)
    df_t = pd.read_csv(p_t, sep="\t")
    df_b = pd.read_csv(p_b, sep="\t").rename(columns={"Brain region": "Tissue"})

    t_existing = set(df_t["Tissue"].unique())
    df_b_extra = df_b[~df_b["Tissue"].isin(t_existing)]
    combined = pd.concat([df_t, df_b_extra], ignore_index=True)

    wide = (
        combined.pivot_table(
            index=["Gene", "Gene name"],
            columns="Tissue",
            values=value_column,
            aggfunc="first",
        )
        .reset_index()
        .rename(columns={"Gene": "Name", "Gene name": "Description"})
    )
    wide.columns.name = None
    prov_t = reference_data.provenance(HPA_TISSUE_GTEX_SOURCE_ID, v_t)
    prov_b = reference_data.provenance(HPA_BRAIN_GTEX_SOURCE_ID, v_b)
    up_sha = {
        f"{HPA_TISSUE_GTEX_SOURCE_ID}:{v_t}": str(prov_t.get("sha256") or ""),
        f"{HPA_BRAIN_GTEX_SOURCE_ID}:{v_b}": str(prov_b.get("sha256") or ""),
    }
    return wide, up_sha


def tissue_rna(
    version: str | None = None,
    *,
    provider: str = "direct",
    qc: GtexSpec | str | bool | Iterable[str] | None = None,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    value_column: str = "TPM",
    include_subcompartments: object = _UNSET,
    hpa_unique_tissues_only: object = _UNSET,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Return a GTEx gene x tissue summary matrix across any version, provider, ``qc``, and ``cleanup``."""
    prov = provider.strip().lower()
    spec = resolve_spec(
        qc=qc,
        cleanup=cleanup,
        normalize=normalize,
        include_subcompartments=include_subcompartments,
        hpa_unique_tissues_only=hpa_unique_tissues_only,
    )
    if spec.sex is not None:
        raise ValueError(
            "Donor sex filtering (sex='XX' or 'XY') requires per-sample GTEx shards "
            "(provider='direct', version='v11')"
        )

    if prov == "direct":
        if value_column.upper() != "TPM":
            raise ValueError(
                f"Direct GTEx releases provide 'TPM' only; value_column={value_column!r} requires provider='hpa'"
            )
        concrete = reference_data.resolve_version(SOURCE_ID, version)
        if spec.is_raw_default:
            path = reference_data.ensure(SOURCE_ID, concrete)
            result = pd.read_csv(path, sep="\t", skiprows=2)
            result.attrs.update(reference_data.provenance(SOURCE_ID, concrete))
            result.attrs["unit"] = "TPM"
            result.attrs["statistic"] = "donor median"
            result.attrs["provider"] = "direct"
            result.attrs["spec"] = spec.to_dict()
            return result

        prov_rec = reference_data.provenance(SOURCE_ID, concrete)
        up_sha = {f"{SOURCE_ID}:{concrete}": str(prov_rec.get("sha256") or "")}
        key, manifest_payload = spec_cache_key(
            provider="direct",
            version=concrete,
            granularity="medians",
            stat="median",
            value_column="TPM",
            spec=spec,
            upstream_sha256=up_sha,
        )

        def _build_direct() -> pd.DataFrame:
            p = reference_data.ensure(SOURCE_ID, concrete)
            raw_df = pd.read_csv(p, sep="\t", skiprows=2)
            return apply_summary_spec(raw_df, spec=spec, provider="direct", value_column="TPM")

        result = bulk_qc.read_or_compute_cached_parquet(
            _effective_derived_cache_dir(),
            key,
            manifest_payload,
            _build_direct,
            use_cache=use_cache,
            spec_attr_name="gtex_spec",
        )
        result.attrs.update(prov_rec)
        result.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
        result.attrs["statistic"] = "donor median"
        result.attrs["provider"] = "direct"
        result.attrs["spec"] = spec.to_dict()
        return result

    if prov == "hpa":
        concrete = reference_data.resolve_version(HPA_TISSUE_GTEX_SOURCE_ID, version)
        prov_t = reference_data.provenance(HPA_TISSUE_GTEX_SOURCE_ID, concrete)
        prov_b = reference_data.provenance(HPA_BRAIN_GTEX_SOURCE_ID, concrete)
        up_sha = {
            f"{HPA_TISSUE_GTEX_SOURCE_ID}:{concrete}": str(prov_t.get("sha256") or ""),
            f"{HPA_BRAIN_GTEX_SOURCE_ID}:{concrete}": str(prov_b.get("sha256") or ""),
        }
        key, manifest_payload = spec_cache_key(
            provider="hpa",
            version=concrete,
            granularity="means",
            stat="mean",
            value_column=value_column,
            spec=spec,
            upstream_sha256=up_sha,
        )

        def _build_hpa() -> pd.DataFrame:
            wide, _ = _load_hpa_gtex_wide(concrete, value_column=value_column)
            if spec.is_raw_default:
                return wide
            return apply_summary_spec(wide, spec=spec, provider="hpa", value_column=value_column)

        result = bulk_qc.read_or_compute_cached_parquet(
            _effective_derived_cache_dir(),
            key,
            manifest_payload,
            _build_hpa,
            use_cache=use_cache,
            spec_attr_name="gtex_spec",
        )
        result.attrs["version"] = concrete
        result.attrs["unit"] = "clean_TPM" if spec.clean_tpm else value_column
        result.attrs["statistic"] = "tissue arithmetic mean (GTEx v8 via HPA)"
        result.attrs["provider"] = "hpa"
        result.attrs["spec"] = spec.to_dict()
        return result

    raise ValueError(f"Unknown GTEx provider {provider!r}; expected 'direct' or 'hpa'")


def donor_sex() -> pd.DataFrame:
    """Return splash-robust chromosomal sex (``XX`` vs ``XY``) for all 946 GTEx v11 donors."""
    return gtex_deconvolution.gtex_donor_sex()


def coverage_grid(*, category: str | None = None, missing_from_gtex: bool = False) -> pd.DataFrame:
    """Return the comprehensive human normal-tissue x source coverage grid (#778)."""
    return normal_tissues.normal_tissue_source_grid(
        category=category, missing_from_gtex=missing_from_gtex
    )


def tissues(
    version: str | None = None,
    *,
    provider: str = "direct",
    include_subcompartments: bool = True,
    include_cell_lines: bool = False,
) -> pd.DataFrame:
    """List available GTEx tissues/shards with sample/donor counts and lineage families."""
    prov = provider.strip().lower()
    if prov == "hpa":
        concrete = reference_data.resolve_version(HPA_TISSUE_GTEX_SOURCE_ID, version)
        rows = [
            {
                "tissue": hpa_t,
                "canonical_site": site,
                "lineage_family": gtex_deconvolution.organ_lineage_family(site),
                "kind": (
                    "hpa_unique_eyegex" if hpa_t in HPA_UNIQUE_GTEX_TISSUES else "hpa_gtex_v8_mean"
                ),
                "provider": "hpa",
                "version": concrete,
            }
            for hpa_t, site in sorted(_HPA_TISSUE_TO_GTEX_SITE.items())
        ]
        return pd.DataFrame(rows)

    if prov != "direct":
        raise ValueError(f"Unknown GTEx provider {provider!r}; expected 'direct' or 'hpa'")
    concrete = reference_data.resolve_version(SOURCE_ID, version)
    tsum = gtex_deconvolution.gtex_deconvolution_tissue_summary().set_index("site")
    rows = []
    for site in primary_shard_sites():
        rec = tsum.loc[site] if site in tsum.index else None
        rows.append(
            {
                "tissue": site,
                "canonical_site": site,
                "lineage_family": gtex_deconvolution.organ_lineage_family(site),
                "kind": "primary_shard",
                "has_individual_samples": concrete == "v11",
                "total_samples": (
                    int(rec["total_samples"]) if rec is not None and concrete == "v11" else None
                ),
                "total_donors": (
                    int(rec["total_donors"]) if rec is not None and concrete == "v11" else None
                ),
                "provider": "direct",
                "version": concrete,
            }
        )
    if concrete == "v11" and include_subcompartments:
        for col in sorted(GTEX_V11_SUBCOMPARTMENT_COLUMNS):
            st = normalize_column_to_site(col)
            rows.append(
                {
                    "tissue": col,
                    "canonical_site": st,
                    "lineage_family": gtex_deconvolution.organ_lineage_family(st),
                    "kind": "v11_subcompartment_median",
                    "has_individual_samples": False,
                    "total_samples": None,
                    "total_donors": None,
                    "provider": "direct",
                    "version": concrete,
                }
            )
    if include_cell_lines:
        for col in ("Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes"):
            st = normalize_column_to_site(col)
            rows.append(
                {
                    "tissue": col,
                    "canonical_site": st,
                    "lineage_family": gtex_deconvolution.organ_lineage_family(st),
                    "kind": "transformed_or_cultured_cells",
                    "has_individual_samples": False,
                    "total_samples": None,
                    "total_donors": None,
                    "provider": "direct",
                    "version": concrete,
                }
            )
    return pd.DataFrame(rows)


def _resolve_shard_pin(tissue: str, version: str = "v11") -> tuple[str, str]:
    site = normalize_column_to_site(tissue)
    pin_key = f"{version}-{site}"
    pins = reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"]
    if pin_key not in pins:
        avail = ", ".join(primary_shard_sites())
        raise ValueError(f"Unknown GTEx {version} tissue shard {tissue!r}; available: {avail}")
    return site, pin_key


@functools.lru_cache(maxsize=2)
def _cached_v11_primary_medians(path_str: str, mtime_ns: int, size_bytes: int) -> pd.DataFrame:
    del mtime_ns, size_bytes
    raw = pd.read_csv(path_str, sep="\t", skiprows=2)
    raw["_gid"] = raw["Name"].astype(str).str.split(".").str[0]
    raw = raw.drop_duplicates("_gid").set_index("_gid")
    site_cols: dict[str, pd.Series] = {}
    for c in raw.columns:
        if (
            c in ("Name", "Description")
            or c in GTEX_V11_SUBCOMPARTMENT_COLUMNS
            or c in CELL_LINE_COLUMNS
        ):
            continue
        site_cols[normalize_column_to_site(c)] = pd.to_numeric(raw[c], errors="coerce").fillna(0.0)
    return pd.DataFrame(site_cols, index=raw.index)


def _v11_primary_medians() -> pd.DataFrame:
    p = reference_data.ensure(SOURCE_ID, "v11")
    st = p.stat()
    return _cached_v11_primary_medians(str(p.resolve()), st.st_mtime_ns, st.st_size)


def _stream_gtex_shard(
    path: Path,
    *,
    genes: set[str] | None,
) -> pd.DataFrame:
    want_all = genes is None
    want_set = {str(g).strip() for g in genes} if genes is not None else set()
    want_bytes = {g.encode("ascii") for g in want_set}

    names: list[str] = []
    descs: list[str] = []
    vecs: list[np.ndarray] = []

    with gzip.open(path, "rb") as fh:
        fh.readline()
        fh.readline()
        header = fh.readline().decode("utf-8").rstrip("\n").split("\t")
        samps = header[2:]
        for line in fh:
            rn_b, _, rest = line.partition(b"\t")
            unver_b = rn_b.split(b".", 1)[0]
            desc_b, _, vals_b = rest.partition(b"\t")
            if not want_all and not (
                rn_b in want_bytes or unver_b in want_bytes or desc_b in want_bytes
            ):
                continue
            names.append(rn_b.decode("ascii"))
            descs.append(desc_b.decode("ascii"))
            vecs.append(np.asarray(vals_b.split(), dtype=np.float64))

    mat = np.vstack(vecs) if vecs else np.zeros((0, len(samps)), dtype=np.float64)
    df = pd.DataFrame(mat, columns=samps)
    df.insert(0, "Description", descs)
    df.insert(0, "Name", names)
    return df


def sample_tpm(
    tissue: str,
    *,
    genes: list[str] | set[str] | tuple[str, ...] | None = None,
    version: str = "v11",
    qc: GtexSpec | str | bool | Iterable[str] | None = None,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    sex: object = _UNSET,
) -> pd.DataFrame:
    """Load per-sample TPM for a GTEx primary shard in strict ``qc -> cleanup -> normalize`` order."""
    spec = resolve_spec(qc=qc, cleanup=cleanup, normalize=normalize, sex=sex)
    site, pin_key = _resolve_shard_pin(tissue, version=version)
    shard_path = reference_data.ensure(SHARD_SOURCE_ID, pin_key)
    want_genes = {str(g).strip() for g in genes} if genes is not None else None
    need_live_anchors = spec.subtract_index_hopping or spec.reject_splash
    anchors_df = (
        gtex_deconvolution.gtex_deconvolution_anchors() if need_live_anchors else pd.DataFrame()
    )

    if spec.clean_tpm or want_genes is None:
        stream_genes = None
    elif need_live_anchors:
        stream_genes = want_genes | set(anchors_df["gene_id"])
    else:
        stream_genes = want_genes

    df_full = _stream_gtex_shard(shard_path, genes=stream_genes)
    samp_cols = [c for c in df_full.columns if c not in ("Name", "Description")]

    dsex_df = donor_sex().set_index("donor_id")["chromosomal_sex"].to_dict()
    swaps_set = (
        set(gtex_deconvolution.gtex_quarantined_swaps()["sample_id"])
        if spec.quarantine_swaps
        else set()
    )

    keep_samps: list[str] = []
    for s in samp_cols:
        if s in swaps_set:
            continue
        if spec.sex is not None:
            did = gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id
            if dsex_df.get(did) != spec.sex:
                continue
        keep_samps.append(s)

    anchor_obs_by_gid: dict[str, np.ndarray] = {}
    if need_live_anchors and len(df_full) > 0 and keep_samps:
        unver_all = df_full["Name"].astype(str).str.split(".").str[0].tolist()
        mat_keep_all = df_full[keep_samps].to_numpy(dtype=np.float64)
        a_set = set(anchors_df["gene_id"])
        for r_i, gid in enumerate(unver_all):
            if gid in a_set and gid not in anchor_obs_by_gid:
                anchor_obs_by_gid[gid] = mat_keep_all[r_i]

    if not spec.clean_tpm and want_genes is not None and stream_genes is not want_genes:
        unver_s = df_full["Name"].astype(str).str.split(".").str[0]
        keep_r = (
            df_full["Name"].isin(want_genes)
            | unver_s.isin(want_genes)
            | df_full["Description"].isin(want_genes)
        )
        df = df_full.loc[keep_r, ["Name", "Description", *keep_samps]].reset_index(drop=True)
    else:
        df = df_full[["Name", "Description", *keep_samps]].copy()

    purity_by_samp = dict.fromkeys(keep_samps, 1.0)
    if (
        (spec.reject_splash or spec.subtract_index_hopping or spec.nnls_deconvolution)
        and keep_samps
        and len(df) > 0
    ):
        ev_all = gtex_deconvolution.gtex_sample_splash_events(include_boundary=False)
        ev_site = ev_all[ev_all["target_site"].eq(site) & ev_all["sample_id"].isin(keep_samps)]
        med_primary: pd.DataFrame | None = None

        unver_rows = df["Name"].astype(str).str.split(".").str[0].tolist()
        vals = df[keep_samps].to_numpy(dtype=np.float64, copy=True)
        samp_idx = {s: j for j, s in enumerate(keep_samps)}
        src_med_vecs: dict[str, np.ndarray] = {}

        def _get_src_vec(ssrc: str) -> np.ndarray | None:
            nonlocal med_primary
            if ssrc in src_med_vecs:
                return src_med_vecs[ssrc]
            if med_primary is None:
                med_primary = _v11_primary_medians()
            if ssrc not in med_primary.columns:
                return None
            v = med_primary[ssrc].reindex(unver_rows).fillna(0.0).to_numpy(dtype=np.float64)
            src_med_vecs[ssrc] = v
            return v

        splashed_pairs: set[tuple[int, str]] = set()
        for sid, grp in ev_site.groupby("sample_id"):
            j = samp_idx[str(sid)]
            by_fam: dict[str, dict] = {}
            for rec in grp.to_dict("records"):
                if rec["regime"] == "tube_swap":
                    continue
                sfam = gtex_deconvolution.organ_lineage_family(str(rec["source_site"]))
                if sfam not in by_fam or float(rec["alpha_hat"]) > float(by_fam[sfam]["alpha_hat"]):
                    by_fam[sfam] = rec

            tot_alpha = sum(float(e["alpha_hat"]) for e in by_fam.values())
            purity_by_samp[str(sid)] = max(0.05, 1.0 - min(0.95, tot_alpha))

            for sfam, e in by_fam.items():
                ssrc = str(e["source_site"])
                a_hat = float(e["alpha_hat"])
                splashed_pairs.add((j, sfam))
                if (
                    spec.reject_splash
                    and a_hat >= 0.005
                    and e["regime"] in {"plate_pipetting_splash", "same_donor_carryover"}
                ):
                    src_vec = _get_src_vec(ssrc)
                    if src_vec is None:
                        continue
                    bleed_est = a_hat * src_vec
                    splash_mask = (
                        np.isfinite(vals[:, j])
                        & (vals[:, j] > 0.0)
                        & (bleed_est >= 0.15)
                        & (bleed_est >= 0.20 * vals[:, j])
                    )
                    vals[splash_mask, j] = np.nan

        if need_live_anchors and anchor_obs_by_gid:
            target_fam = gtex_deconvolution.organ_lineage_family(site)
            best_hop: dict[tuple[int, str], tuple[float, float, np.ndarray, bool]] = {}
            for src_site, grp in anchors_df.groupby("source_site", sort=True):
                src_fam = gtex_deconvolution.organ_lineage_family(str(src_site))
                if src_site == site or src_fam == target_fam:
                    continue
                if gtex_deconvolution.is_anatomical_boundary(str(src_site), site):
                    continue
                keep_a = [
                    g in anchor_obs_by_gid and site not in {x for x in str(ex).split(";") if x}
                    for g, ex in zip(grp["gene_id"], grp["excluded_target_sites"])
                ]
                sub_grp = grp.loc[keep_a]
                if len(sub_grp) < 3:
                    continue
                src_vec = _get_src_vec(str(src_site))
                if src_vec is None:
                    continue
                a_mat = np.vstack([anchor_obs_by_gid[g] for g in sub_grp["gene_id"]])
                a_meds = np.maximum(sub_grp["source_median_tpm"].to_numpy(dtype=np.float64), 1.0)
                det = np.isfinite(a_mat) & (a_mat >= 0.08)
                cand_js = np.where(det.sum(axis=0) >= 3)[0]
                for j_raw in cand_js:
                    j = int(j_raw)
                    if (j, src_fam) in splashed_pairs:
                        continue
                    dm = det[:, j]
                    r_vec = a_mat[dm, j] / a_meds[dm]
                    a_hat = float(np.median(r_vec))
                    if not (2e-5 <= a_hat < 5e-3):
                        continue
                    lr = np.log10(np.maximum(r_vec, 1e-12))
                    q75, q25 = np.percentile(lr, [75.0, 25.0])
                    if float(q75 - q25) > 0.65:
                        continue
                    a_lo = float(min(a_hat, np.percentile(r_vec, 25.0)))
                    samp_dsex = dsex_df.get(
                        gtex_deconvolution.parse_gtex_sample_barcode(keep_samps[j]).donor_id,
                        "unassessed",
                    )
                    cross_sex = (
                        samp_dsex == "XX"
                        and str(src_site) in gtex_deconvolution.MALE_SPECIFIC_SITES
                    ) or (
                        samp_dsex == "XY"
                        and str(src_site) in gtex_deconvolution.FEMALE_SPECIFIC_SITES
                    )
                    prev = best_hop.get((j, src_fam))
                    if prev is None or a_hat > prev[0]:
                        best_hop[(j, src_fam)] = (a_hat, a_lo, src_vec, cross_sex)

            for (j, _sfam), (a_hat, a_lo, src_vec, cross_sex) in best_hop.items():
                if spec.reject_splash and cross_sex:
                    bleed_est = a_hat * src_vec
                    splash_mask = (
                        np.isfinite(vals[:, j])
                        & (vals[:, j] > 0.0)
                        & (bleed_est >= 0.10)
                        & (bleed_est >= 0.20 * vals[:, j])
                    )
                    vals[splash_mask, j] = np.nan
                    continue
                if spec.subtract_index_hopping:
                    bleed_est = a_lo * src_vec
                    sub_mask = (
                        np.isfinite(vals[:, j])
                        & (vals[:, j] > 0.0)
                        & (bleed_est > 0.0)
                        & (bleed_est >= 0.25 * vals[:, j])
                    )
                    vals[sub_mask, j] = np.fmax(0.0, vals[sub_mask, j] - bleed_est[sub_mask])

        df[keep_samps] = vals

    if spec.clean_tpm and keep_samps and len(df) > 0:
        cur_vals = df[keep_samps].to_numpy(dtype=np.float64)
        nan_mask = np.isnan(cur_vals)
        gene_tbl = pd.DataFrame({"Ensembl_Gene_ID": df["Name"].astype(str).str.split(".").str[0]})
        clean_arr = normalization.clean_tpm(
            df[keep_samps].fillna(0.0), gene_tbl, on_empty_compartment="record"
        ).to_numpy(dtype=np.float64, copy=True)
        clean_arr[nan_mask] = np.nan
        df[keep_samps] = clean_arr

    if spec.clean_tpm and want_genes is not None:
        unver_s = df["Name"].astype(str).str.split(".").str[0]
        keep_r = (
            df["Name"].isin(want_genes)
            | unver_s.isin(want_genes)
            | df["Description"].isin(want_genes)
        )
        df = df.loc[keep_r].reset_index(drop=True)

    df.attrs["tissue"] = site
    df.attrs["version"] = version
    df.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
    df.attrs["spec"] = spec.to_dict()
    df.attrs["purity_weights"] = purity_by_samp
    df.attrs["donor_sex"] = {
        s: dsex_df.get(gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id, "unassessed")
        for s in keep_samps
    }
    return df


def _match_summary_columns(mat_columns: Iterable[str], requested: list[str]) -> list[str]:
    val_cols = [c for c in mat_columns if c not in ("Name", "Description")]
    by_exact = {c: c for c in val_cols}
    by_lower = {c.lower(): c for c in val_cols}
    by_site: dict[str, list[str]] = {}
    for c in val_cols:
        by_site.setdefault(normalize_column_to_site(c), []).append(c)

    out: list[str] = []
    seen: set[str] = set()
    for req in requested:
        r = str(req).strip()
        if r in by_exact:
            matched = [by_exact[r]]
        elif r.lower() in by_lower:
            matched = [by_lower[r.lower()]]
        else:
            st = normalize_column_to_site(r)
            if st in by_site:
                matched = by_site[st]
            else:
                raise ValueError(f"Unknown GTEx tissue {req!r}; available: {val_cols}")
        for m in matched:
            if m not in seen:
                seen.add(m)
                out.append(m)
    return out


def tissue_summary(
    tissues_select: str | list[str] | tuple[str, ...] | None = None,
    *,
    genes: list[str] | set[str] | tuple[str, ...] | None = None,
    version: str | None = None,
    provider: str = "direct",
    stat: object = _UNSET,
    qc: object = _UNSET,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    sex: object = _UNSET,
    value_column: str = "TPM",
    include_subcompartments: object = _UNSET,
    hpa_unique_tissues_only: object = _UNSET,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Summarize GTEx donor expression across one or all tissues (#778)."""
    prov = provider.strip().lower()
    if version is None or version == "latest":
        version = "v23" if prov == "hpa" else "v11"

    if isinstance(qc, GtexSpec):
        spec = resolve_spec(
            qc=qc,
            cleanup=cleanup,
            normalize=normalize,
            sex=sex,
            include_subcompartments=include_subcompartments,
            hpa_unique_tissues_only=hpa_unique_tissues_only,
        )
    else:
        spec = resolve_spec(
            qc="all" if qc is _UNSET else qc,  # type: ignore[arg-type]
            cleanup="all" if cleanup is _UNSET else cleanup,
            normalize=normalize,
            sex=sex,
            include_subcompartments=include_subcompartments,
            hpa_unique_tissues_only=hpa_unique_tissues_only,
        )

    if stat is _UNSET:
        stat_norm = "mean" if prov == "hpa" else ("median" if version != "v11" else "p95")
    else:
        stat_norm = str(stat).strip().lower()
    valid_stats = {"mean", "median", "p95", "p99", "max", "host_exposure"}
    if stat_norm not in valid_stats:
        raise ValueError(f"Unknown stat {stat!r}; choose from {sorted(valid_stats)}")

    if prov == "hpa":
        if stat_norm != "mean":
            raise ValueError(
                f"HPA GTEx redistributes tissue arithmetic means only; got stat={stat!r} "
                "(use stat='mean' or provider='direct', version='v11' for donor quantiles)"
            )
    elif prov == "direct":
        if value_column.upper() != "TPM":
            raise ValueError(
                f"Direct GTEx releases provide 'TPM' only; value_column={value_column!r} requires provider='hpa'"
            )
        if stat_norm == "mean":
            raise ValueError(
                "Direct GTEx summaries use donor quantiles ('median', 'p95', 'p99', 'max', 'host_exposure'), not 'mean'"
            )

    if prov == "hpa" or (prov == "direct" and version != "v11"):
        if prov == "direct" and stat_norm != "median":
            raise ValueError(
                f"GTEx {version} pins donor medians only; per-sample {stat!r} requires version='v11'"
            )
        mat = tissue_rna(
            version=version,
            provider=prov,
            qc=spec,
            value_column=value_column,
            use_cache=use_cache,
        )
        if genes is not None:
            wset = {str(g).strip() for g in genes}
            unver = mat["Name"].astype(str).str.split(".").str[0]
            mat = mat[
                mat["Name"].isin(wset) | unver.isin(wset) | mat["Description"].isin(wset)
            ].reset_index(drop=True)
        if tissues_select is not None:
            req = [tissues_select] if isinstance(tissues_select, str) else list(tissues_select)
            keep_c = _match_summary_columns(mat.columns, req) if req else []
            mat = mat[["Name", "Description", *keep_c]].copy()
        return mat

    if tissues_select is None:
        target_sites = list(primary_shard_sites())
    elif isinstance(tissues_select, str):
        target_sites = [_resolve_shard_pin(tissues_select, version=version)[0]]
    else:
        target_sites = [_resolve_shard_pin(t, version=version)[0] for t in tissues_select]

    if not target_sites:
        empty = pd.DataFrame(columns=["Name", "Description"])
        empty.attrs["provider"] = "direct"
        empty.attrs["version"] = version
        empty.attrs["statistic"] = stat_norm
        empty.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
        empty.attrs["spec"] = spec.to_dict()
        return empty

    can_cache_full = use_cache and genes is None and tissues_select is None

    def _compute_summary() -> pd.DataFrame:
        base_meta: pd.DataFrame | None = None
        site_cols: dict[str, np.ndarray] = {}
        q_map = {"median": 0.50, "p95": 0.95, "p99": 0.99}

        for site in target_sites:
            sdf = sample_tpm(site, genes=genes, version=version, qc=spec)
            if base_meta is None:
                base_meta = sdf[["Name", "Description"]].copy()
            s_cols = [c for c in sdf.columns if c not in ("Name", "Description")]
            if not s_cols or len(sdf) == 0:
                site_cols[site] = np.full(len(sdf), np.nan, dtype=np.float64)
                continue

            vals = sdf[s_cols].to_numpy(dtype=np.float64)
            pw_map = sdf.attrs.get("purity_weights", {})
            p_wts = np.array([float(pw_map.get(s, 1.0)) for s in s_cols], dtype=np.float64)

            donors = [gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id for s in s_cols]
            unique_donors = list(dict.fromkeys(donors))
            if len(unique_donors) != len(s_cols):
                d_idx_map: dict[str, list[int]] = {}
                for j_i, d_id in enumerate(donors):
                    d_idx_map.setdefault(d_id, []).append(j_i)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    d_vals = np.vstack(
                        [np.nanmedian(vals[:, idxs], axis=1) for idxs in d_idx_map.values()]
                    ).T
                d_wts = np.array([p_wts[idxs].mean() for idxs in d_idx_map.values()])
            else:
                d_vals = vals
                d_wts = p_wts

            out_vec = np.empty(len(sdf), dtype=np.float64)
            for r_i in range(len(sdf)):
                row_v = d_vals[r_i]
                fin = np.isfinite(row_v)
                if not fin.any():
                    out_vec[r_i] = np.nan
                    continue
                rv = row_v[fin]
                rw = d_wts[fin]
                if stat_norm == "max":
                    out_vec[r_i] = float(rv.max())
                elif stat_norm == "host_exposure":
                    sv = np.sort(rv)
                    nd = sv.size
                    if nd >= 50:
                        p99 = float(np.percentile(sv, 99.0))
                        out_vec[r_i] = min(
                            float(sv[-1]), max(p99, float(sv[-2])) if nd < 100 else p99
                        )
                    elif nd >= 30:
                        out_vec[r_i] = float(sv[-2])
                    else:
                        out_vec[r_i] = float(sv[-1])
                elif spec.nnls_deconvolution:
                    out_vec[r_i] = bulk_qc.purity_weighted_harrell_davis_quantile(
                        rv, rw, quantile=q_map[stat_norm]
                    )
                else:
                    out_vec[r_i] = float(np.percentile(rv, 100.0 * q_map[stat_norm]))
            site_cols[site] = out_vec

        assert base_meta is not None
        return pd.concat([base_meta, pd.DataFrame(site_cols, index=base_meta.index)], axis=1)

    pins = reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"]
    up_sha = {f"{SHARD_SOURCE_ID}:v11-{s}": pins[f"v11-{s}"]["sha256"] for s in target_sites}
    key, manifest_payload = spec_cache_key(
        provider="direct",
        version=version,
        granularity="samples",
        stat=stat_norm,
        value_column="TPM",
        spec=spec,
        upstream_sha256=up_sha,
    )
    out = bulk_qc.read_or_compute_cached_parquet(
        _effective_derived_cache_dir(),
        key,
        manifest_payload,
        _compute_summary,
        use_cache=can_cache_full,
        spec_attr_name="gtex_spec",
    )
    out.attrs["provider"] = "direct"
    out.attrs["version"] = version
    out.attrs["statistic"] = stat_norm
    out.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
    out.attrs["spec"] = spec.to_dict()
    return out


__all__ = [
    "apply_summary_spec",
    "coverage_grid",
    "donor_sex",
    "sample_tpm",
    "tissue_rna",
    "tissue_summary",
    "tissues",
    "versions",
]
