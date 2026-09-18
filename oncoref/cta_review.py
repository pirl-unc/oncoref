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

    Memoized: the returned frame is a shared read-only view, so callers must
    not mutate it in place. The two public frames copy before handing out.

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
    """The three pinned atlas tables, memoized as shared read-only views.

    Neither the dict nor its frames may be mutated in place: they back every
    later summary in the process, and a rebinding here would go unnoticed.
    """
    # Explicit version: future changes to the default must not silently relabel
    # old measurements or change the baseline of a reviewed exception. The IHC
    # table goes through hpa's version-aware public accessor rather than reaching
    # past it; hpa_rna_consensus and hpa_single_cell take no version, so pinning
    # those needs the raw read. The filter below copies the IHC frame, so this
    # does not share hpa's memory -- it shares its version resolution.
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


#: A label counts as routinely surveyed when at least this share of the
#: source's genes carry a measurement for it. The two populations are two
#: orders of magnitude apart in HPA v23 -- every routine label covers >=95% of
#: genes, every special-study label <1% -- so the cut is not delicate.
_ROUTINE_LABEL_SHARE = 0.5


def _routine_labels(table: pd.DataFrame, column: str) -> frozenset[str]:
    """Labels the source runs for most genes, not every label it ever used.

    HPA's IHC table mixes its standard tissue panel with special-study labels
    measured for a handful of genes (substantia nigra and sole of foot for one
    gene each, retina for 113 of 13,468). Treating those as part of a scope
    would mean no gene is ever fully surveyed in it.
    """
    genes = table["Gene"].nunique()
    per_label = table.groupby(column)["Gene"].nunique()
    return frozenset(per_label[per_label >= _ROUTINE_LABEL_SHARE * genes].index)


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


def _rna_stats(rows: pd.DataFrame, expected_rows: int | None) -> dict:
    """Summarize an RNA subset, counting what was measured against what exists.

    ``expected_rows`` is how many measurements the source holds for this scope,
    and is required: a gene observed in 3 of 10 brain regions must not report a
    brain maximum as though the region had been surveyed. Pass ``None`` only
    where no count defines the scope, so that choice is visible in review. The
    counts come from the source, so a release that drops a tissue lowers the
    denominator with it; there is no curated all-tissues list to check against.
    """
    values = pd.to_numeric(rows["nTPM"], errors="coerce")
    measured = values[values.ge(0)]
    if not len(measured):
        status = "unavailable"
    elif expected_rows is not None and len(measured) < expected_rows:
        status = "incomplete"
    else:
        status = "measured"
    return {
        "status": status,
        "max_ntpm": measured.max(),
        "measured_rows": len(measured),
        "expected_rows": expected_rows,
    }


def _ihc_stats(rows: pd.DataFrame, expected_tissues: int | None) -> dict:
    """Summarize an IHC subset, counting stained tissues against the scope.

    ``expected_tissues`` is how many tissues this scope is routinely surveyed
    in, so a gene stained in fewer of them reports ``incomplete`` rather than a
    clean negative. It must count what the source actually runs, not every
    label that appears somewhere in it: counting rarely-run special-study
    labels would put the threshold out of reach, making ``not_detected``
    arithmetically impossible and leaving a detection as the only way out of
    ``incomplete``. The shortfall between the tissues a group asks for and the
    ones the source maps at all is a separate axis, carried by
    ``mapping_coverage`` and :func:`cta_atlas_coverage`.
    """
    detected = rows["Level"].isin(_IHC_DETECTED_LEVELS)
    negative = rows["Level"].eq(_IHC_NEGATIVE_LEVEL)
    measured = detected | negative
    measured_tissues = rows.loc[measured, "Tissue"].nunique()
    surveyed = expected_tissues is None or measured_tissues >= expected_tissues
    if detected.any():
        status = "detected"
    elif len(rows) and negative.all() and surveyed:
        status = "not_detected"
    elif measured.any():
        status = "incomplete"
    else:
        status = "unavailable"
    return {
        "status": status,
        "measured_tissues": measured_tissues,
        "expected_tissues": expected_tissues,
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


def cta_atlas_coverage() -> pd.DataFrame:
    """How completely this release covers each safety group, tissue by tissue.

    A property of the release rather than of any gene: HPA v23 covers 8 of the
    14 requested brain regions by IHC -- 6 in full, basal ganglia by caudate
    alone and midbrain by two nuclei -- and 10 by RNA, so a brain result of any
    kind speaks for neither spinal cord nor thalamus. Stated once here rather
    than repeated into every candidate row, where it would crowd out per-gene
    facts.

    One row per requested tissue, not per group, so the counts reconcile and a
    region represented by a single substructure is not tallied as covered.
    ``modality`` matches the summary's column prefixes (``rna``, ``ihc``).
    """
    records = []
    for modality, prefix in (("bulk_rna", "rna"), ("ihc", "ihc")):
        for group, resolution in _resolve(modality).items():
            for mapping in resolution.mappings:
                records.append(
                    {
                        "safety_group": group,
                        "modality": prefix,
                        "requested_tissue": mapping.requested_tissue,
                        "coverage_level": mapping.coverage_level,
                        "mapping_kind": mapping.mapping_kind,
                        "source_tissues": ";".join(mapping.source_tissues),
                        "group_coverage_state": resolution.coverage_state,
                        "source_name": _SOURCES[modality],
                        "source_version": _VERSION,
                        "source_url": resolution.source_url,
                    }
                )
    return pd.DataFrame(records)


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
    # How many measurements the source itself holds for each RNA scope, so a
    # gene observed in only part of a scope is reported as incomplete rather
    # than summarized as if the whole scope had been surveyed.
    # Denominators come from what each source routinely surveys. Using every
    # label that appears anywhere would put them out of reach: no gene is
    # stained in HPA's rare special-study labels, so "not detected" would become
    # arithmetically impossible and -- because a detection short-circuits the
    # survey check -- detection would be the only escape from "incomplete",
    # inverting the gap field into a marker for the most concerning genes.
    all_tissues = _routine_labels(tables["bulk_rna"], "Tissue")
    all_cell_types = _routine_labels(tables["single_cell_rna"], "Cell type")
    ihc_tissues = _routine_labels(tables["ihc"], "Tissue")
    expected = {
        "bulk_rna": len(all_tissues),
        "somatic_rna": len(all_tissues - set(NON_SOMATIC_TISSUES)),
        "single_cell_rna": len(all_cell_types),
        "cardiomyocyte_rna": len(
            {c for c in all_cell_types if c.casefold() in _CARDIOMYOCYTE_LABELS}
        ),
        "somatic_ihc": len(ihc_tissues - set(ALL_REPRODUCTIVE_TISSUES)),
    }
    records = []
    for gene_id in universe["Ensembl_Gene_ID"]:
        gene = {key: groups[key].get(gene_id, tables[key].iloc[:0]) for key in tables}
        rna, ihc, single = gene["bulk_rna"], gene["ihc"], gene["single_cell_rna"]
        out = {"Ensembl_Gene_ID": gene_id, "atlas_source_version": _VERSION}
        out.update(source_urls)
        somatic_rna = rna.loc[~rna["Tissue"].isin(NON_SOMATIC_TISSUES)]
        somatic_ihc = ihc.loc[~ihc["Tissue"].isin(ALL_REPRODUCTIVE_TISSUES)]
        _add_stats(out, "bulk_rna", _rna_stats(rna, expected["bulk_rna"]))
        _add_stats(out, "somatic_rna", _rna_stats(somatic_rna, expected["somatic_rna"]))
        _add_stats(out, "somatic_ihc", _ihc_stats(somatic_ihc, expected["somatic_ihc"]))
        _add_stats(out, "single_cell_rna", _rna_stats(single, expected["single_cell_rna"]))
        cardio = single.loc[single["Cell type"].str.casefold().isin(_CARDIOMYOCYTE_LABELS)]
        _add_stats(out, "cardiomyocyte_rna", _rna_stats(cardio, expected["cardiomyocyte_rna"]))
        # Heart tissue comes from the same resolution as the safety groups, so a
        # release that renames the label raises there instead of quietly
        # reporting this gene's cardiomyocytes as unassayed.
        heart_ihc = ihc.loc[
            ihc["Tissue"].isin(resolutions["ihc"]["heart"].source_tissues)
            & ihc["Cell type"].str.casefold().isin(_CARDIOMYOCYTE_LABELS)
        ]
        _add_stats(
            out,
            "cardiomyocyte_ihc",
            _ihc_stats(heart_ihc, len(resolutions["ihc"]["heart"].source_tissues)),
        )
        warnings = []
        if out["somatic_ihc_detected_rows"]:
            warnings.append("somatic_protein_detected")
        if out["somatic_rna_max_ntpm"] > 0:
            warnings.append("somatic_rna_estimate")
            # Keyed on the absence of detection rather than on a clean negative:
            # somatic IHC is never surveyed across every tissue, so requiring
            # the literal "not_detected" would silence every discordance.
            if out["somatic_ihc_detected_rows"] == 0 and out["somatic_ihc_status"] != "unavailable":
                warnings.append("rna_ihc_discordance")
                # The token alone is not a completeness claim, so say when the
                # negative half of it rests on a partial survey.
                if out["somatic_ihc_status"] == "incomplete":
                    warnings.append("rna_ihc_discordance_partial_survey")
        for group in SAFETY_TISSUE_GROUPS:
            for modality, stats in (("rna", _rna_stats), ("ihc", _ihc_stats)):
                # Both modalities go through their own resolution. Matching the
                # conceptual names against the source directly would count only
                # the labels that happen to coincide and call that the group.
                resolution = resolutions["bulk_rna" if modality == "rna" else "ihc"][group]
                rows = rna if modality == "rna" else ihc
                prefix = f"{group}_{modality}"
                subset = rows.loc[rows["Tissue"].isin(resolution.source_tissues)]
                # Only the mapped labels the source routinely runs can be
                # expected of a gene; the rest are a mapping gap, not a
                # measurement this gene is missing.
                routine = all_tissues if modality == "rna" else ihc_tissues
                _add_stats(
                    out, prefix, stats(subset, len(set(resolution.source_tissues) & routine))
                )
                out[f"{prefix}_mapping_coverage"] = resolution.coverage_state
                out[f"{prefix}_mapped_tissues"] = ";".join(resolution.source_tissues)
                out[f"{prefix}_unmapped_tissues"] = ";".join(resolution.unavailable_tissues)
            if out[f"{group}_rna_max_ntpm"] >= SAFETY_NTPM_THRESHOLD:
                warnings.append(f"{group}_rna_ge_{SAFETY_NTPM_THRESHOLD:g}_ntpm")
        # Per-gene gaps only. The fixed mapping limitation of the release is the
        # same for every gene, so listing it here would leave the field never
        # empty and unable to distinguish a gene with missing data from one
        # measured everywhere; it is reported per release in
        # :func:`cta_atlas_coverage` and per row in ``atlas_coverage_limits``.
        # Derived from every status the record holds rather than a hand-kept
        # list, which has already fallen behind twice as scopes were added.
        # Named with the status, so "never assayed" stays distinguishable from
        # "assayed in part" without cross-reading a second column.
        gaps = [
            f"{key[: -len('_status')]}_{value}"
            for key, value in out.items()
            if key.endswith("_status") and value in {"unavailable", "incomplete"}
        ]
        limits = []
        for group in SAFETY_TISSUE_GROUPS:
            for modality in ("rna", "ihc"):
                prefix = f"{group}_{modality}"
                # Coverage describes what was assayed, not what was concluded,
                # so it is stated whether or not this gene was detected there.
                if out[f"{prefix}_mapping_coverage"] != "complete":
                    limits.append(f"{prefix}_{out[f'{prefix}_mapping_coverage']}_mapping")
        out["atlas_warning_codes"] = ";".join(warnings)
        out["atlas_evidence_gaps"] = ";".join(gaps)
        out["atlas_coverage_limits"] = ";".join(limits)
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
        columns=lambda column: (
            column
            if column == "Ensembl_Gene_ID" or column.startswith("warning_")
            else f"warning_{column}"
        )
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
            parts.append(f"Missing data: {row['atlas_evidence_gaps'].replace(';', ', ')}.")
        if row["atlas_coverage_limits"]:
            parts.append(
                f"Atlas covers only part of: {row['atlas_coverage_limits'].replace(';', ', ')}."
            )
        parts.append("Atlas measurements do not establish peptide presentation or clinical safety.")
        return " ".join(parts)

    out["evidence_summary"] = out.apply(summary, axis=1)
    return out.copy()


def cta_normal_tissue_evidence() -> pd.DataFrame:
    """Long-form HPA v23 measurements for the complete CTA candidate universe.

    One row per gene, modality and tissue or cell type, covering all available
    tissues rather than only the five safety groups. ``reported_zero`` is an RNA
    estimate and ``not_detected`` an IHC annotation; the two are never merged.
    Rows HPA does not provide are absent here and counted in the summary's
    ``atlas_evidence_gaps``. Single-cell types are aggregated across organs, so
    they carry no tissue rather than being attributed to one. Memoized, and
    returned as a copy so in-place mutation cannot corrupt the shared result.
    Downloads the three pinned atlas sources on first use if not cached.
    """
    return _normal_tissue_frame().copy()


def cta_evidence_summary() -> pd.DataFrame:
    """One comparable evidence summary per CTA, watchlist, or clinical candidate.

    Every candidate receives bulk RNA, all-tissue somatic IHC, the five
    safety-tissue groups and cardiomyocyte RNA/IHC on the pinned v23 baseline,
    plus ``discovery_tier``, ``atlas_warning_codes``, per-gene missing data in
    ``atlas_evidence_gaps``, the release's fixed mapping limits in
    ``atlas_coverage_limits`` (see :func:`cta_atlas_coverage`), a
    ``*_review_status`` per reviewed modality and a prose ``evidence_summary``.
    Detailed donor, proteomics, isoform, peptide-presentation and clinical
    reviews stay in :func:`cta_reviewed_evidence`, where absence of a row means
    not reviewed. No negative assay result promotes a gene into a broader set.
    Memoized, and returned as a copy so in-place mutation cannot corrupt the
    shared result.
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
