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

from functools import cache, lru_cache

import pandas as pd

from . import cta, hpa, reference_data
from .cta_tissues import (
    ALL_REPRODUCTIVE_TISSUES,
    NON_SOMATIC_TISSUES,
    PROTEIN_DETECTED_LEVELS,
    SAFETY_NTPM_THRESHOLD,
    SAFETY_TISSUE_GROUPS,
)
from .load_dataset import _register_derived_cache, get_data

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
# HPA's IHC vocabulary is wider than the Low/Medium/High set the strict CTA
# filter scores. Ascending and Descending are gradient staining, so protein is
# present; omitting them would drop a somatic protein detection in the
# permissive direction. "Not representative" is an unusable annotation rather
# than a negative result, so it counts as neither detection nor non-detection.
_IHC_DETECTED_LEVELS: frozenset[str] = PROTEIN_DETECTED_LEVELS | {"Ascending", "Descending"}
_IHC_NEGATIVE_LEVEL = "Not detected"
# HPA writes either label for this cell type depending on release.
_CARDIOMYOCYTE_LABELS: frozenset[str] = frozenset({"myocytes", "cardiomyocytes"})


def _strip_version(ids: pd.Series) -> pd.Series:
    """Ensembl IDs join without their version suffix, from every table.

    Applied to each curated table too, not only the universe: a reviewed row
    written as ``ENSG00000126890.14`` would otherwise merge onto nothing and
    read as *not reviewed* rather than as an error.
    """
    return ids.astype(str).str.split(".").str[0]


@lru_cache(maxsize=1)
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


@lru_cache(maxsize=1)
def _atlas_tables() -> dict[str, pd.DataFrame]:
    # Explicit version: future changes to the default must not silently relabel
    # old measurements or change the baseline of a reviewed exception. The IHC
    # table goes through hpa's version-aware accessor so this module shares its
    # cache instead of holding a second copy; hpa_rna_consensus and
    # hpa_single_cell take no version, so pinning those needs the raw read.
    tables = {
        "ihc": hpa.hpa_normal_tissue(_VERSION),
        "bulk_rna": hpa._read_hpa(_SOURCES["bulk_rna"], _VERSION),
        "single_cell_rna": hpa._read_hpa(_SOURCES["single_cell_rna"], _VERSION),
    }
    # HPA marks "no antibody data" with a row whose tissue, cell type and level
    # are all null, carrying only a Reliability. Those are not measurements: a
    # null tissue passes every somatic filter, so keeping them would pad the
    # long-form output and let a formatting artifact turn a clean non-detection
    # into an incomplete one.
    ihc = tables["ihc"]
    tables["ihc"] = ihc.loc[ihc["Level"].notna()]
    return tables


def _source_url(modality: str) -> str:
    """The URL the cached artifact actually came from, not the configured one.

    Goes through ``provenance`` so a renamed source or version raises the
    module's own error, and so the URL reported here matches the one
    :func:`hpa.resolve_safety_tissue_group` records for the same artifact.
    """
    return reference_data.provenance(_SOURCES[modality], _VERSION)["url"]


@lru_cache(maxsize=1)
def _normal_tissue_frame() -> pd.DataFrame:
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
            out.loc[out["ihc_level"].eq(_IHC_NEGATIVE_LEVEL), "measurement_status"] = "not_detected"
            out.loc[out["ihc_level"].isin(_IHC_DETECTED_LEVELS), "measurement_status"] = "detected"
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


def _ntpm_phrase(value: float) -> str:
    """Render an nTPM figure without flattening the three cases into a number.

    A zero is an estimate of zero, not an observed absence, and an absent
    measurement is neither; prose that printed both as "0 nTPM" would discard
    the distinction the long-form table keeps.
    """
    if pd.isna(value):
        return "unavailable"
    if value == 0:
        return "estimated 0 nTPM"
    return f"{value:g} nTPM"


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
    detected = rows["Level"].isin(_IHC_DETECTED_LEVELS)
    negative = rows["Level"].eq(_IHC_NEGATIVE_LEVEL)
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


@cache
def _resolve(modality: str) -> dict[str, hpa.SafetyTissueResolution]:
    """Conceptual safety groups mapped to one source's own tissue labels.

    Resolved per modality rather than once for IHC: the RNA consensus and the
    IHC table cover different regions, so reusing one mapping for both would
    misreport whichever it was not built for.
    """
    return {
        group: hpa.resolve_safety_tissue_group(
            group,
            source_name=_SOURCES[modality],
            source_version=_VERSION,
            require_complete=False,
        )
        for group in SAFETY_TISSUE_GROUPS
    }


def _synthesize_atlas(universe: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Pure synthesis over supplied tables; no source fetches or tier promotion."""
    ids = set(universe["Ensembl_Gene_ID"])
    # Narrow before grouping: the atlas tables span ~20k genes each, and
    # grouping them whole materializes a frame per gene to then discard all but
    # these few hundred.
    groups = {
        key: dict(iter(frame.loc[frame["Gene"].isin(ids)].groupby("Gene")))
        # iter(): a groupby's ``keys`` attribute is the grouping label, so a bare
        # dict() takes the mapping protocol path and tries to call that string.
        for key, frame in tables.items()
    }
    resolutions = {"bulk_rna": _resolve("bulk_rna"), "ihc": _resolve("ihc")}
    source_urls = {f"{modality}_source_url": _source_url(modality) for modality in tables}
    records = []
    for gene_id in universe["Ensembl_Gene_ID"]:
        gene = {key: groups[key].get(gene_id, tables[key].iloc[:0]) for key in tables}
        rna, ihc, single = gene["bulk_rna"], gene["ihc"], gene["single_cell_rna"]
        out = {"Ensembl_Gene_ID": gene_id, "atlas_source_version": _VERSION}
        out.update(source_urls)
        somatic_rna = rna.loc[~rna["Tissue"].isin(NON_SOMATIC_TISSUES)]
        somatic_ihc = ihc.loc[~ihc["Tissue"].isin(ALL_REPRODUCTIVE_TISSUES)]
        _add_stats(out, "bulk_rna", _rna_stats(rna))
        _add_stats(out, "somatic_rna", _rna_stats(somatic_rna))
        _add_stats(out, "somatic_ihc", _ihc_stats(somatic_ihc))
        _add_stats(out, "single_cell_rna", _rna_stats(single))
        cardio = single.loc[single["Cell type"].str.casefold().isin(_CARDIOMYOCYTE_LABELS)]
        _add_stats(out, "cardiomyocyte_rna", _rna_stats(cardio))
        # Heart tissue comes from the same resolution as the safety groups, so a
        # release that renames the label raises there instead of quietly
        # reporting this gene's cardiomyocytes as unassayed.
        heart_ihc = ihc.loc[
            ihc["Tissue"].isin(resolutions["ihc"]["heart"].source_tissues)
            & ihc["Cell type"].str.casefold().isin(_CARDIOMYOCYTE_LABELS)
        ]
        _add_stats(out, "cardiomyocyte_ihc", _ihc_stats(heart_ihc))
        warnings = []
        if out["somatic_ihc_detected_rows"]:
            warnings.append("somatic_protein_detected")
        if out["somatic_rna_max_ntpm"] > 0:
            warnings.append("somatic_rna_estimate")
            if out["somatic_ihc_status"] == "not_detected":
                warnings.append("rna_ihc_discordance")
        for group in SAFETY_TISSUE_GROUPS:
            for modality, stats in (("rna", _rna_stats), ("ihc", _ihc_stats)):
                # Both modalities go through their own resolution. Matching the
                # conceptual names against the source directly would count only
                # the labels that happen to coincide and call that the group.
                resolution = resolutions["bulk_rna" if modality == "rna" else "ihc"][group]
                rows = rna if modality == "rna" else ihc
                prefix = f"{group}_{modality}"
                _add_stats(
                    out, prefix, stats(rows.loc[rows["Tissue"].isin(resolution.source_tissues)])
                )
                out[f"{prefix}_mapping_coverage"] = resolution.coverage_state
                out[f"{prefix}_mapped_tissues"] = ";".join(resolution.source_tissues)
                out[f"{prefix}_unmapped_tissues"] = ";".join(resolution.unavailable_tissues)
            if out[f"{group}_rna_max_ntpm"] >= SAFETY_NTPM_THRESHOLD:
                warnings.append(f"{group}_rna_ge_{SAFETY_NTPM_THRESHOLD:g}_ntpm")
        gaps = []
        for key in ("bulk_rna", "somatic_ihc", "single_cell_rna", "cardiomyocyte_ihc"):
            if out[f"{key}_status"] in {"unavailable", "incomplete"}:
                gaps.append(key)
        # Coverage is reported for what was assayed, not for what was concluded.
        # HPA v23 stains 9 of the 14 requested brain regions and measures RNA in
        # 10, so neither a detection nor a non-detection there speaks for spinal
        # cord or thalamus; gating the caveat on the verdict would let a positive
        # finding in one region imply the whole group had been looked at.
        for group in SAFETY_TISSUE_GROUPS:
            for modality in ("rna", "ihc"):
                prefix = f"{group}_{modality}"
                if out[f"{prefix}_mapping_coverage"] != "complete":
                    gaps.append(f"{prefix}_partial_mapping")
                if out[f"{prefix}_status"] in {"unavailable", "incomplete"}:
                    gaps.append(f"{prefix}_{out[f'{prefix}_status']}")
        out["atlas_warning_codes"] = ";".join(warnings)
        out["atlas_evidence_gaps"] = ";".join(gaps)
        records.append(out)
    return pd.DataFrame(records)


@lru_cache(maxsize=1)
def _evidence_summary_frame() -> pd.DataFrame:
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
        parts = [
            f"{row['discovery_tier']}; HPA {_VERSION} heart RNA {_ntpm_phrase(row['heart_rna_max_ntpm'])}; "
            f"cardiomyocyte IHC {row['cardiomyocyte_ihc_status']}; "
            f"somatic IHC {row['somatic_ihc_status']}."
        ]
        # The one-line field has to carry the caveats, not just the reassuring
        # measurements: a reader who sees only this must not conclude that the
        # atlas was complete or that nothing was flagged.
        if row["atlas_warning_codes"]:
            parts.append(f"Warnings: {row['atlas_warning_codes'].replace(';', ', ')}.")
        if row["atlas_evidence_gaps"]:
            parts.append(f"Coverage gaps: {row['atlas_evidence_gaps'].replace(';', ', ')}.")
        parts.append("Atlas measurements do not establish peptide presentation or clinical safety.")
        return " ".join(parts)

    out["evidence_summary"] = out.apply(summary, axis=1)
    return out.copy()


def cta_normal_tissue_evidence() -> pd.DataFrame:
    """Long-form HPA v23 measurements for the complete CTA candidate universe.

    See :func:`_normal_tissue_frame`. Returns a defensive copy, so a caller that
    mutates its frame in place cannot corrupt the shared result.
    """
    return _normal_tissue_frame().copy()


def cta_evidence_summary() -> pd.DataFrame:
    """One comparable evidence summary per CTA, watchlist, or clinical candidate.

    See :func:`_evidence_summary_frame`. Returns a defensive copy, so a caller
    that mutates its frame in place cannot corrupt the shared result.
    """
    return _evidence_summary_frame().copy()


# Synthesis is a pure function of the bundled CSVs and the pinned atlas files,
# so it is computed once per process. Registered for clearing so fixture swaps
# in tests still invalidate these alongside the dataset cache.
for _clear in (
    _universe.cache_clear,
    _atlas_tables.cache_clear,
    _resolve.cache_clear,
    _normal_tissue_frame.cache_clear,
    _evidence_summary_frame.cache_clear,
):
    _register_derived_cache(_clear)
