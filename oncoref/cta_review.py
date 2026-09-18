# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

"""Comparable normal-tissue evidence for every CTA and watchlist candidate.

Atlas measurements, reviewed discovery exceptions, and unresolved evidence are
kept separate. Neither a zero RNA estimate nor negative IHC establishes absence
of peptide presentation. The fixed v23 atlas baseline is not silently mixed with
newer, individually reviewed measurements.
"""

from __future__ import annotations

import pandas as pd

from . import cta, hpa, reference_data
from .cta_tissues import (
    ALL_REPRODUCTIVE_TISSUES,
    NON_SOMATIC_TISSUES,
    PROTEIN_DETECTED_LEVELS,
    SAFETY_NTPM_THRESHOLD,
    SAFETY_TISSUE_GROUPS,
)
from .load_dataset import get_data

_VERSION = "v23"
# Reviewed modalities the summary surfaces. A row whose modality is absent here
# would load without error and never reach a caller, so the set is asserted
# against the reviewed-evidence table rather than left implicit.
REVIEWED_MODALITIES: tuple[str, ...] = (
    "donor_rna",
    "cell_resolved_rna",
    "mass_spectrometry",
    "isoform",
    "peptide_presentation",
    "clinical",
)
_SOURCES = {
    "bulk_rna": "hpa_rna_consensus",
    "ihc": "hpa_normal_tissue",
    "single_cell_rna": "hpa_single_cell",
}


def _strip_version(ids: pd.Series) -> pd.Series:
    """Ensembl IDs join without their version suffix, from every table.

    Applied to each curated table too, not only the universe: a reviewed row
    written as ``ENSG00000126890.14`` would otherwise merge onto nothing and
    read as *not reviewed* rather than as an error.
    """
    return ids.astype(str).str.split(".").str[0]


def _universe() -> pd.DataFrame:
    """Identity and provenance only, one row per candidate gene.

    Deliberately carries no measurement columns. Concatenating the curation
    tables whole would seat their own HPA numbers, of unstated version and
    different derivation, beside this module's version-pinned columns -- a
    ``rna_heart_max_ntpm`` next to our ``heart_rna_max_ntpm`` -- which is the
    confusion the summary exists to remove. Join :func:`cta.cta_evidence` on
    the gene ID for curation detail.
    """
    origins = {
        "cta_table": cta.cta_evidence(),
        "watchlist": cta.cta_candidate_references(),
        "clinical_references": cta.cta_clinical_target_references(),
    }
    frames = []
    for origin, frame in origins.items():
        rows = frame[["Symbol", "Ensembl_Gene_ID"]].copy()
        rows["candidate_origin"] = origin
        frames.append(rows)
    out = pd.concat(frames, ignore_index=True)
    out["Ensembl_Gene_ID"] = _strip_version(out["Ensembl_Gene_ID"])
    return out.drop_duplicates("Ensembl_Gene_ID").sort_values("Symbol").reset_index(drop=True)


def cta_reviewed_evidence() -> pd.DataFrame:
    """Source-anchored supplemental observations; absence means not reviewed.

    Unlike the systematic atlas summary, this table is not comprehensive for
    all genes. Each row states assay, scope, source version, and limitations.
    """
    return get_data("cta-reviewed-evidence").copy()


def _atlas_tables() -> dict[str, pd.DataFrame]:
    # Explicit version: future changes to the default must not silently relabel
    # old measurements or change the baseline of a reviewed exception.
    return {modality: hpa._read_hpa(source, _VERSION) for modality, source in _SOURCES.items()}


def _source_url(modality: str) -> str:
    return reference_data.REFERENCE_SOURCES[_SOURCES[modality]]["urls"][_VERSION]


def cta_normal_tissue_evidence() -> pd.DataFrame:
    """Long-form HPA v23 measurements for the complete CTA candidate universe.

    Includes all available tissues/cell types, not only the five safety groups.
    ``reported_zero`` is an RNA estimate; ``not_detected`` is an IHC annotation.
    Missing source rows are absent here and explicitly marked in the summary.
    Single-cell types are aggregated across organs: their tissue is left blank
    rather than assigning e.g. every fibroblast measurement to the heart.
    Downloads the three pinned atlas sources on first use if not cached.
    """
    ids = set(_universe()["Ensembl_Gene_ID"])
    frames = []
    for modality, raw in _atlas_tables().items():
        rows = raw.loc[raw["Gene"].isin(ids)].copy()
        out = rows.rename(
            columns={
                "Gene": "Ensembl_Gene_ID",
                "Gene name": "Symbol",
                "Tissue": "tissue",
                "Cell type": "cell_type",
                "Level": "ihc_level",
                "Reliability": "ihc_reliability",
                "nTPM": "rna_ntpm",
            }
        )
        out["modality"] = modality
        out["source_name"] = _SOURCES[modality]
        out["source_version"] = _VERSION
        out["source_url"] = _source_url(modality)
        out["unit"] = "IHC category" if modality == "ihc" else "nTPM"
        out["measurement_status"] = "unavailable"
        if modality == "ihc":
            out.loc[out["ihc_level"].eq("Not detected"), "measurement_status"] = "not_detected"
            out.loc[out["ihc_level"].isin(PROTEIN_DETECTED_LEVELS), "measurement_status"] = (
                "detected"
            )
        else:
            values = pd.to_numeric(out["rna_ntpm"], errors="coerce")
            out.loc[values.eq(0), "measurement_status"] = "reported_zero"
            out.loc[values.gt(0), "measurement_status"] = "positive_estimate"
        frames.append(out)
    return pd.concat(frames, ignore_index=True, sort=False).reindex(
        columns=[
            "Symbol",
            "Ensembl_Gene_ID",
            "modality",
            "tissue",
            "cell_type",
            "rna_ntpm",
            "ihc_level",
            "ihc_reliability",
            "measurement_status",
            "unit",
            "source_name",
            "source_version",
            "source_url",
        ]
    )


def _labels(values: pd.Series) -> str:
    return ";".join(sorted(set(values.dropna().astype(str)) - {""}))


def _rna_stats(rows: pd.DataFrame) -> dict:
    values = pd.to_numeric(rows["nTPM"], errors="coerce")
    measured = values[values.ge(0)]
    return {
        "status": "measured" if len(measured) else "unavailable",
        "max_ntpm": measured.max(),
        "measured_rows": len(measured),
    }


def _ihc_stats(rows: pd.DataFrame) -> dict:
    detected = rows["Level"].isin(PROTEIN_DETECTED_LEVELS)
    negative = rows["Level"].eq("Not detected")
    measured = detected | negative
    if detected.any():
        status = "detected"
    elif len(rows) and negative.all():
        status = "not_detected"
    elif measured.any():
        status = "incomplete"
    else:
        status = "unavailable"
    return {
        "status": status,
        "measured_rows": int(measured.sum()),
        "unavailable_rows": int((~measured).sum()),
        "detected_rows": int(detected.sum()),
        "not_detected_rows": int(negative.sum()),
        "detected_tissues": _labels(rows.loc[detected, "Tissue"]),
        "levels": _labels(rows["Level"]),
        "reliability": _labels(rows.loc[measured, "Reliability"]),
    }


def _add_stats(out: dict, prefix: str, stats: dict) -> None:
    out.update({f"{prefix}_{key}": value for key, value in stats.items()})


def _synthesize_atlas(universe: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Pure synthesis over supplied tables; no source fetches or tier promotion."""
    # iter(): a groupby's ``keys`` attribute is the grouping label, so a bare
    # dict() takes the mapping protocol path and tries to call that string.
    groups = {key: dict(iter(frame.groupby("Gene"))) for key, frame in tables.items()}
    resolutions = {
        group: hpa.resolve_safety_tissue_group(
            group,
            source_name="hpa_normal_tissue",
            source_version=_VERSION,
            require_complete=False,
        )
        for group in SAFETY_TISSUE_GROUPS
    }
    records = []
    for gene_id in universe["Ensembl_Gene_ID"]:
        gene = {key: groups[key].get(gene_id, tables[key].iloc[:0]) for key in tables}
        rna, ihc, single = gene["bulk_rna"], gene["ihc"], gene["single_cell_rna"]
        out = {"Ensembl_Gene_ID": gene_id, "atlas_source_version": _VERSION}
        for modality in tables:
            out[f"{modality}_source_url"] = _source_url(modality)
        somatic_rna = rna.loc[~rna["Tissue"].isin(NON_SOMATIC_TISSUES)]
        somatic_ihc = ihc.loc[~ihc["Tissue"].isin(ALL_REPRODUCTIVE_TISSUES)]
        _add_stats(out, "bulk_rna", _rna_stats(rna))
        _add_stats(out, "somatic_rna", _rna_stats(somatic_rna))
        _add_stats(out, "somatic_ihc", _ihc_stats(somatic_ihc))
        _add_stats(out, "single_cell_rna", _rna_stats(single))
        cardio = single.loc[single["Cell type"].str.casefold().eq("cardiomyocytes")]
        _add_stats(out, "cardiomyocyte_rna", _rna_stats(cardio))
        heart_ihc = ihc.loc[
            ihc["Tissue"].eq("heart muscle") & ihc["Cell type"].str.casefold().eq("myocytes")
        ]
        # HPA versions use either myocytes or cardiomyocytes for this cell type.
        heart_ihc = pd.concat(
            [
                heart_ihc,
                ihc.loc[
                    ihc["Tissue"].eq("heart muscle")
                    & ihc["Cell type"].str.casefold().eq("cardiomyocytes")
                ],
            ]
        )
        _add_stats(out, "cardiomyocyte_ihc", _ihc_stats(heart_ihc))
        warnings = []
        if out["somatic_ihc_detected_rows"]:
            warnings.append("somatic_protein_detected")
        if out["somatic_rna_max_ntpm"] > 0:
            warnings.append("somatic_rna_estimate")
            if out["somatic_ihc_status"] == "not_detected":
                warnings.append("rna_ihc_discordance")
        for group, tissues in SAFETY_TISSUE_GROUPS.items():
            subset = rna.loc[rna["Tissue"].isin(tissues)]
            _add_stats(out, f"{group}_rna", _rna_stats(subset))
            resolution = resolutions[group]
            _add_stats(
                out,
                f"{group}_ihc",
                _ihc_stats(ihc.loc[ihc["Tissue"].isin(resolution.source_tissues)]),
            )
            out[f"{group}_ihc_mapping_coverage"] = resolution.coverage_state
            out[f"{group}_ihc_mapped_tissues"] = ";".join(resolution.source_tissues)
            out[f"{group}_ihc_unmapped_tissues"] = ";".join(resolution.unavailable_tissues)
            if out[f"{group}_rna_max_ntpm"] >= SAFETY_NTPM_THRESHOLD:
                warnings.append(f"{group}_rna_ge_{SAFETY_NTPM_THRESHOLD:g}_ntpm")
        gaps = []
        for key in ("bulk_rna", "somatic_ihc", "single_cell_rna", "cardiomyocyte_ihc"):
            if out[f"{key}_status"] in {"unavailable", "incomplete"}:
                gaps.append(key)
        # A negative read on a partially mapped group is not group-wide absence:
        # HPA v23 stains 9 of the 14 requested brain regions, so "not_detected"
        # there excludes neither spinal cord nor thalamus. Positive detection
        # stands on its own, so only non-detections carry the coverage caveat.
        for group in SAFETY_TISSUE_GROUPS:
            if (
                out[f"{group}_ihc_mapping_coverage"] != "complete"
                and out[f"{group}_ihc_status"] != "detected"
            ):
                gaps.append(f"{group}_ihc_partial_mapping")
        out["atlas_warning_codes"] = ";".join(warnings)
        out["atlas_evidence_gaps"] = ";".join(gaps)
        records.append(out)
    return pd.DataFrame(records)


def cta_evidence_summary() -> pd.DataFrame:
    """One comparable evidence summary per CTA, watchlist, or clinical candidate.

    All candidates receive bulk RNA, all-tissue somatic IHC, five safety-group
    summaries, and cardiomyocyte RNA/IHC. Detailed donor, proteomics, isoform,
    peptide-presentation and clinical reviews are explicitly tracked separately.
    No negative assay result automatically promotes a gene into a broader set.
    """
    universe = _universe()
    out = universe.merge(_synthesize_atlas(universe, _atlas_tables()), on="Ensembl_Gene_ID")
    strict = cta.cta_gene_ids()
    filtered = cta.cta_filtered_gene_ids()
    refs = cta.cta_warning_references()
    refs["Ensembl_Gene_ID"] = _strip_version(refs["Ensembl_Gene_ID"])
    warnings = set(refs["Ensembl_Gene_ID"])
    out["discovery_tier"] = "excluded"
    out.loc[out["candidate_origin"].ne("cta_table"), "discovery_tier"] = "candidate"
    out.loc[out["Ensembl_Gene_ID"].isin(filtered - strict), "discovery_tier"] = "low_expression"
    out.loc[out["Ensembl_Gene_ID"].isin(warnings), "discovery_tier"] = "warning"
    out.loc[out["Ensembl_Gene_ID"].isin(strict), "discovery_tier"] = "strict"
    refs = refs.drop(columns="Symbol").rename(
        columns=lambda column: column if column == "Ensembl_Gene_ID" else f"warning_{column}"
    )
    out = out.merge(refs, on="Ensembl_Gene_ID", how="left", validate="one_to_one")
    reviewed = cta_reviewed_evidence()
    reviewed["Ensembl_Gene_ID"] = _strip_version(reviewed["Ensembl_Gene_ID"])
    for modality in REVIEWED_MODALITIES:
        subset = reviewed.loc[reviewed["modality"].eq(modality)]
        notes = subset.groupby("Ensembl_Gene_ID")["finding"].agg(" | ".join)
        out[f"{modality}_review_status"] = (
            out["Ensembl_Gene_ID"]
            .map(dict.fromkeys(notes.index, "reviewed"))
            .fillna("not_reviewed")
        )
        out[f"{modality}_review"] = out["Ensembl_Gene_ID"].map(notes)

    def summary(row):
        heart = row["heart_rna_max_ntpm"]
        value = f"{heart:g} nTPM" if pd.notna(heart) else "unavailable"
        return (
            f"{row['discovery_tier']}; HPA {_VERSION} heart RNA {value}; "
            f"cardiomyocyte IHC {row['cardiomyocyte_ihc_status']}; "
            f"somatic IHC {row['somatic_ihc_status']}. "
            "Atlas measurements do not establish peptide presentation or clinical safety."
        )

    out["evidence_summary"] = out.apply(summary, axis=1)
    return out.copy()
