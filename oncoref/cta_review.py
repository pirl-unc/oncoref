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

import numpy as np
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


_LABEL_COLUMNS = {"bulk_rna": "Tissue", "ihc": "Tissue", "single_cell_rna": "Cell type"}


def _routine_labels_in(table: pd.DataFrame, column: str) -> frozenset[str]:
    """Labels the source runs for most genes, not every label it ever used.

    HPA's IHC table mixes its standard tissue panel with special-study labels
    measured for a handful of genes (substantia nigra and sole of foot for one
    gene each, retina for 113 of 13,468). Treating those as part of a scope
    would mean no gene is ever fully surveyed in it.

    Takes the table rather than a modality so :func:`_synthesize_atlas` stays a
    pure function of the tables it is handed; reading the cached atlas here
    would make synthesis ignore its own arguments.
    """
    genes = table["Gene"].nunique()
    per_label = table.groupby(column)["Gene"].nunique()
    return frozenset(per_label[per_label >= _ROUTINE_LABEL_SHARE * genes].index)


@cache
def _routine_labels(modality: str) -> frozenset[str]:
    """Memoized routine labels of one pinned source, for callers outside
    synthesis. Scanning the full tables costs ~0.36s, and
    :func:`cta_atlas_coverage` needs the same answer the denominators use."""
    return _routine_labels_in(_atlas_tables()[modality], _LABEL_COLUMNS[modality])


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


def _rna_stats(rows: pd.DataFrame, surveyed: frozenset[str], label_column: str) -> dict:
    """Summarize an RNA subset, counting what was measured against what exists.

    ``surveyed`` is the set of labels the source runs for most genes in this
    scope. Both sides of the comparison are counted over it: a gene observed in
    3 of 10 brain regions must not report a brain maximum as though the region
    had been surveyed, and measurements in labels outside the set must not make
    up the shortfall. The counts come from the source, so a release that drops
    a label lowers the denominator with it; there is no curated list to check
    against. ``max_ntpm`` deliberately spans every measurement, in scope or
    not, since the highest value anywhere is the safety-relevant one.
    """
    values = pd.to_numeric(rows["nTPM"], errors="coerce")
    measured = values[values.ge(0)]
    in_scope = measured[rows.loc[measured.index, label_column].isin(surveyed)]
    if not len(measured):
        status = "unavailable"
    elif len(in_scope) < len(surveyed):
        status = "incomplete"
    else:
        status = "measured"
    return {
        "status": status,
        "max_ntpm": measured.max(),
        "measured_rows": len(in_scope),
        "expected_rows": len(surveyed),
    }


def _ihc_stats(rows: pd.DataFrame, surveyed: frozenset[str]) -> dict:
    """Summarize an IHC subset, counting stained tissues against the scope.

    ``surveyed`` is the set of tissues this scope is routinely stained in, so a
    gene scored in fewer of them reports ``incomplete`` rather than a clean
    negative. Both sides are counted over that set. It must hold what the
    source actually runs, not every label appearing somewhere in it: a
    denominator inflated with rarely-run special-study labels puts the
    threshold out of reach, making ``not_detected`` arithmetically impossible
    and leaving a detection as the only way out of ``incomplete``; a numerator
    that counts them lets a gene scored only in special-study labels claim a
    clean negative having surveyed none of the routine ones. The shortfall
    between the tissues a group asks for and the ones the source maps at all is
    a separate axis, carried by ``mapping_coverage`` and
    :func:`cta_atlas_coverage`.
    """
    detected = rows["Level"].isin(_IHC_DETECTED_LEVELS)
    negative = rows["Level"].eq(_IHC_NEGATIVE_LEVEL)
    measured = detected | negative
    in_scope = measured & rows["Tissue"].isin(surveyed)
    measured_tissues = rows.loc[in_scope, "Tissue"].nunique()
    fully_surveyed = measured_tissues >= len(surveyed)
    if detected.any():
        status = "detected"
    elif len(rows) and negative.all() and fully_surveyed:
        status = "not_detected"
    elif measured.any():
        status = "incomplete"
    else:
        status = "unavailable"
    return {
        "status": status,
        "measured_tissues": measured_tissues,
        "expected_tissues": len(surveyed),
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


def _survey_state(mapping, surveyed: tuple[str, ...]) -> str:
    """Whether a requested tissue is unmapped, mapped-but-rarely-run, or run.

    A single boolean would put the first two together, which is the
    distinction this column exists to draw: thalamus has no HPA IHC label at
    all, while retina has an exact one that HPA stains for under 1% of genes.
    """
    if not mapping.source_tissues:
        return "not_mapped"
    if not surveyed:
        return "mapped_not_surveyed"
    return "surveyed"


def _surveyed_level(mapping, surveyed: tuple[str, ...]) -> str:
    """Restate a mapping's coverage counting only routinely surveyed labels.

    A region whose every label is a rare special-study one is not covered at
    all, however exactly it maps; one that keeps some of its labels is at best
    partial. Only a region whose labels are all routine keeps its mapped level.
    """
    if not surveyed:
        return "unavailable"
    if len(surveyed) < len(mapping.source_tissues):
        return "partial"
    return mapping.coverage_level


def _aggregate_level(levels) -> str:
    """Collapse per-tissue coverage into one group verdict, as hpa does.

    Kept beside the per-tissue levels so the group column cannot keep
    reporting the mapping figure after the rows below it gained a stricter
    one: for brain IHC that is the difference between 8 of 14 requested
    regions mapped and 4 of 14 actually surveyed.
    """
    distinct = set(levels)
    if distinct == {"complete"}:
        return "complete"
    if distinct == {"unavailable"}:
        return "unavailable"
    return "partial"


def cta_atlas_coverage() -> pd.DataFrame:
    """How completely this release covers each safety group, tissue by tissue.

    A property of the release rather than of any gene: HPA v23 maps 8 of the 14
    requested brain regions for IHC and routinely surveys only 4 of those, so a
    brain IHC result speaks for none of the other ten -- spinal cord and
    thalamus among them. RNA maps 10 and surveys all 10, missing thalamus,
    medulla oblongata, pons and white matter but covering spinal cord.
    ``group_mapped_regions`` and ``group_surveyed_regions`` carry both counts,
    since the aggregate verdict is ``partial`` either way and would hide the
    gap. Stated once here rather than repeated into every candidate row, where
    it would crowd out per-gene facts.

    One row per requested tissue, not per group, so the counts reconcile and a
    region represented by a single substructure is not tallied as covered.
    ``modality`` matches the summary's column prefixes (``rna``, ``ihc``).

    A label can be mapped and still not be part of the panel the source runs
    for most genes: HPA maps choroid plexus, dorsal raphe, hypothalamus, retina
    and substantia nigra for brain IHC but stains each for under 1% of genes.
    ``survey_state`` separates the three cases a single flag would blur --
    ``not_mapped``, ``mapped_not_surveyed``, ``surveyed`` -- so the population
    this column exists for is one filter rather than a conjunction, and
    ``source_tissues`` still names the labels behind a ``mapped_not_surveyed``
    region. ``surveyed_coverage_level`` restates the level counting only routine
    labels, which is the basis the summary's ``expected_tissues`` uses, so the
    two surfaces cannot disagree about whether a tissue was really covered.
    """
    records = []
    for modality, prefix in (("bulk_rna", "rna"), ("ihc", "ihc")):
        routine = _routine_labels(modality)
        for group, resolution in _resolve(modality).items():
            levels = {
                mapping.requested_tissue: _surveyed_level(
                    mapping, tuple(t for t in mapping.source_tissues if t in routine)
                )
                for mapping in resolution.mappings
            }
            group_surveyed = _aggregate_level(levels.values())
            # The aggregate verdict collapses to "partial" either way for
            # brain, so carry the counts too: 8 of 14 requested regions are
            # mapped but only 4 are routinely surveyed, and that gap is the
            # whole reason the denominators differ from the mapping.
            mapped_regions = sum(m.coverage_level != "unavailable" for m in resolution.mappings)
            surveyed_regions = sum(level != "unavailable" for level in levels.values())
            for mapping in resolution.mappings:
                surveyed = tuple(t for t in mapping.source_tissues if t in routine)
                surveyed_level = levels[mapping.requested_tissue]
                records.append(
                    {
                        "safety_group": group,
                        "modality": prefix,
                        "requested_tissue": mapping.requested_tissue,
                        "coverage_level": mapping.coverage_level,
                        "mapping_kind": mapping.mapping_kind,
                        "source_tissues": ";".join(mapping.source_tissues),
                        "survey_state": _survey_state(mapping, surveyed),
                        "surveyed_source_tissues": ";".join(surveyed),
                        "surveyed_coverage_level": surveyed_level,
                        "group_coverage_state": resolution.coverage_state,
                        "group_surveyed_coverage_state": group_surveyed,
                        "group_requested_regions": len(resolution.mappings),
                        "group_mapped_regions": mapped_regions,
                        "group_surveyed_regions": surveyed_regions,
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
    # Each scope's surveyed labels: what this source runs for most genes,
    # narrowed to the scope. Both sides of every completeness check are counted
    # over these. Using every label that appears anywhere would put the
    # denominators out of reach, making "not detected" arithmetically
    # impossible and leaving a detection as the only escape from "incomplete";
    # counting out-of-scope labels in the numerator would let a gene scored
    # only in special-study labels claim a clean negative instead.
    routine = {
        modality: _routine_labels_in(frame, _LABEL_COLUMNS[modality])
        for modality, frame in tables.items()
    }
    all_tissues = routine["bulk_rna"]
    all_cell_types = routine["single_cell_rna"]
    ihc_tissues = routine["ihc"]
    scope = {
        "bulk_rna": all_tissues,
        "somatic_rna": all_tissues - set(NON_SOMATIC_TISSUES),
        "single_cell_rna": all_cell_types,
        "cardiomyocyte_rna": frozenset(
            c for c in all_cell_types if c.casefold() in _CARDIOMYOCYTE_LABELS
        ),
        "somatic_ihc": ihc_tissues - set(ALL_REPRODUCTIVE_TISSUES),
        "cardiomyocyte_ihc": frozenset(resolutions["ihc"]["heart"].source_tissues) & ihc_tissues,
    }
    records = []
    for gene_id in universe["Ensembl_Gene_ID"]:
        gene = {key: groups[key].get(gene_id, tables[key].iloc[:0]) for key in tables}
        rna, ihc, single = gene["bulk_rna"], gene["ihc"], gene["single_cell_rna"]
        out = {"Ensembl_Gene_ID": gene_id, "atlas_source_version": _VERSION}
        out.update(source_urls)
        somatic_rna = rna.loc[~rna["Tissue"].isin(NON_SOMATIC_TISSUES)]
        somatic_ihc = ihc.loc[~ihc["Tissue"].isin(ALL_REPRODUCTIVE_TISSUES)]
        _add_stats(out, "bulk_rna", _rna_stats(rna, scope["bulk_rna"], "Tissue"))
        _add_stats(out, "somatic_rna", _rna_stats(somatic_rna, scope["somatic_rna"], "Tissue"))
        _add_stats(out, "somatic_ihc", _ihc_stats(somatic_ihc, scope["somatic_ihc"]))
        _add_stats(
            out, "single_cell_rna", _rna_stats(single, scope["single_cell_rna"], "Cell type")
        )
        cardio = single.loc[single["Cell type"].str.casefold().isin(_CARDIOMYOCYTE_LABELS)]
        _add_stats(
            out, "cardiomyocyte_rna", _rna_stats(cardio, scope["cardiomyocyte_rna"], "Cell type")
        )
        # Heart tissue comes from the same resolution as the safety groups, so a
        # release that renames the label raises there instead of quietly
        # reporting this gene's cardiomyocytes as unassayed.
        heart_ihc = ihc.loc[
            ihc["Tissue"].isin(resolutions["ihc"]["heart"].source_tissues)
            & ihc["Cell type"].str.casefold().isin(_CARDIOMYOCYTE_LABELS)
        ]
        # Same intersection as the group loop: a release that moved heart muscle
        # out of the routine panel must not leave this one scope unreachable.
        _add_stats(out, "cardiomyocyte_ihc", _ihc_stats(heart_ihc, scope["cardiomyocyte_ihc"]))
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
        unsurveyed: dict[str, set[str]] = {}
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
                surveyed = routine["bulk_rna" if modality == "rna" else "ihc"]
                mapped = set(resolution.source_tissues)
                in_scope = frozenset(mapped & surveyed)
                _add_stats(
                    out,
                    prefix,
                    stats(subset, in_scope, "Tissue")
                    if modality == "rna"
                    else stats(subset, in_scope),
                )
                out[f"{prefix}_mapping_coverage"] = resolution.coverage_state
                out[f"{prefix}_mapped_tissues"] = ";".join(resolution.source_tissues)
                out[f"{prefix}_unmapped_tissues"] = ";".join(resolution.unavailable_tissues)
                unsurveyed[prefix] = mapped - surveyed
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
                # Mapped but outside the panel the source runs for most genes,
                # so excluded from the denominator. A release-level fact, hence
                # a token here rather than ten near-empty per-gene columns;
                # cta_atlas_coverage names the labels per tissue.
                if unsurveyed[prefix]:
                    limits.append(f"{prefix}_unsurveyed_labels")
        out["atlas_warning_codes"] = ";".join(warnings)
        out["atlas_evidence_gaps"] = ";".join(gaps)
        out["atlas_coverage_limits"] = ";".join(limits)
        records.append(out)
    return pd.DataFrame(records)


def _review_notes(reviewed: pd.DataFrame, modality: str) -> pd.Series:
    """Per-gene reviewed findings for one modality, joined in curation order.

    A blank ``finding`` is rejected by name rather than joined: ``astype(str)``
    would turn it into the string "nan", and joining it raw raises deep inside
    the aggregation and takes down every gene rather than the one row at fault.
    """
    subset = reviewed.loc[reviewed["modality"].eq(modality)]
    if subset["finding"].isna().any():
        blank = sorted(subset.loc[subset["finding"].isna(), "Ensembl_Gene_ID"])
        raise ValueError(f"cta-reviewed-evidence rows with no finding: {blank}")
    return subset.groupby("Ensembl_Gene_ID")["finding"].agg(" | ".join)


@lru_cache(maxsize=1)
def _evidence_summary_frame() -> pd.DataFrame:
    """One comparable evidence summary per CTA, watchlist, or clinical candidate.

    All candidates receive bulk RNA, all-tissue somatic IHC, five safety-group
    summaries, and cardiomyocyte RNA/IHC. Detailed donor, proteomics, isoform,
    peptide-presentation and clinical reviews are explicitly tracked separately.
    No negative assay result automatically promotes a gene into a broader set.
    """
    universe = _universe()
    # Guarded like the warning merge below: a synthesis that dropped or
    # duplicated a candidate would otherwise shorten the frame in silence,
    # against a contract of one row per candidate.
    out = universe.merge(
        _synthesize_atlas(universe, _atlas_tables()),
        on="Ensembl_Gene_ID",
        how="left",
        validate="one_to_one",
    )
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
        notes = _review_notes(reviewed, modality)
        out[f"{modality}_review_status"] = np.where(
            out["Ensembl_Gene_ID"].isin(notes.index), "reviewed", "not_reviewed"
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
    return out


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
    _routine_labels.cache_clear,
    _resolve.cache_clear,
    _normal_tissue_frame.cache_clear,
    _evidence_summary_frame.cache_clear,
):
    _register_derived_cache(_clear)
del _clear
