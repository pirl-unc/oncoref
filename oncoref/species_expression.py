# Licensed under the Apache License, Version 2.0.
"""Native expression matrices with species contracts and auditable sample metadata."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from . import bulk_qc
from ._table_io import read_csv
from .expression_contract import (
    ExpressionContract,
    ExpressionMeasurement,
    MeasurementTransform,
    validate_expression_join,
)
from .expression_records import RECORD_KINDS, ExpressionRecords
from .expression_registry import (
    expression_source,
    expression_source_registry_entries,
    expression_source_registry_path,
)
from .reference_identity import _identity
from .species import resolve_species
from .species_assets import DEFAULT_MAX_BYTES, ensure_asset, species_asset

ALL_SPECIES_BULK_QC: tuple[str, ...] = (
    "pooled_libraries",
    "broad_background",
    "swaps",
    "splash",
)

SPECIES_BULK_QC_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_SPECIES_BULK_QC),
    "none": frozenset(),
    "raw": frozenset(),
}

ALL_SPECIES_BULK_CLEANUP: tuple[str, ...] = (
    "index_hopping",
    "deconvolution",
)

SPECIES_BULK_CLEANUP_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_SPECIES_BULK_CLEANUP),
    "none": frozenset(),
    "raw": frozenset(),
}

_SPECIES_OR_RE = re.compile(r"^(?:OR\d+[A-Z]+\d+|OLFR\d+)$", re.IGNORECASE)


def _tissue_to_species_panels(tissue_label: object) -> frozenset[str]:
    if tissue_label is None or pd.isna(tissue_label):
        return frozenset()
    low = str(tissue_label).strip().lower().replace("_", " ")
    matched: set[str] = set()
    if "liver" in low or "gallbladder" in low or "hepat" in low:
        matched.add("liver")
    if "pancrea" in low or "islet" in low:
        matched.add("pancreas")
    if "lung" in low or "bronch" in low or "trachea" in low:
        matched.add("lung")
    if "adrenal" in low or "suprarenal" in low:
        matched.add("adrenal_gland")
    elif "kidney" in low or re.search(r"\brenal\b", low):
        matched.add("kidney")
    if "thyroid" in low or "parathyroid" in low:
        matched.add("thyroid")
    if "pituitary" in low or "hypophys" in low:
        matched.add("pituitary")
    if (
        any(
            k in low
            for k in (
                "brain",
                "cerebell",
                "cortex",
                "hippocamp",
                "hypothalam",
                "amygdala",
                "spinal",
                "retina",
                "nerve",
            )
        )
        and "kidney" not in low
        and not re.search(r"\brenal\b", low)
        and "adrenal" not in low
        and "suprarenal" not in low
    ):
        matched.add("brain")
    if (
        any(k in low for k in ("muscle", "heart", "myocard", "ventricle", "atrium", "tongue"))
        and "smooth" not in low
    ):
        matched.add("striated_muscle")
    if "stomach" in low or "gastric" in low or "duodenum" in low:
        matched.add("stomach")
    if "salivary" in low or "parotid" in low or "submandibular" in low:
        matched.add("salivary_gland")
    if any(k in low for k in ("skin", "esophag", "mucosa", "tongue", "cervix", "vagina", "tonsil")):
        matched.add("squamous_epithelium")
    if "testis" in low or "testes" in low or "epididym" in low:
        matched.add("testis")
    if "ovary" in low or "ovaries" in low or "oviduct" in low or "fallopian" in low:
        matched.add("ovary")
    if "placenta" in low or "trophoblast" in low:
        matched.add("placenta")
    return frozenset(matched)


def _feature_symbols(features: pd.DataFrame, values_index: pd.Index) -> pd.Series:
    for col in ("GeneSymbol", "symbol", "gene_symbol", "Gene name", "gene_name"):
        if col in features.columns:
            s = features[col].reindex(values_index).fillna("").astype(str)
            return s.where(
                s.str.strip().ne(""), pd.Series(values_index.astype(str), index=values_index)
            )
    return pd.Series(values_index.astype(str), index=values_index)


def _species_native_organ_panels(
    values: pd.DataFrame,
    features: pd.DataFrame,
    samples: pd.DataFrame,
) -> tuple[dict[str, pd.Series], pd.DataFrame, dict[str, frozenset[str]]]:
    """Construct within-dataset species-native organ reference panels & profiles.

    Never borrows human expression numbers; derives reference medians exclusively
    from samples labeled with each organ inside ``(values, samples)``. When symbol
    columns are absent (e.g. Ensembl-ID-only Bgee tables), falls back to empirical
    within-dataset organ-restricted markers (median >= 50 and >= 20x off-organ P95).
    """
    syms = _feature_symbols(features, values.index).str.strip().str.upper()
    native_map = {
        str(sid): _tissue_to_species_panels(samples.loc[sid, "tissue"]) for sid in samples.index
    }
    panels: dict[str, pd.Series] = {}
    prof_cols: dict[str, pd.Series] = {}
    all_panel_names = sorted(bulk_qc.VERTEBRATE_ORGAN_ANCHOR_SYMBOLS)
    organ_medians: dict[str, pd.Series] = {}
    for pname in all_panel_names:
        native_samps = [
            sid for sid, pset in native_map.items() if pname in pset and sid in values.columns
        ]
        if native_samps:
            organ_medians[pname] = values.loc[:, native_samps].median(axis=1)

    for pname, organ_med in organ_medians.items():
        anchor_syms = bulk_qc.VERTEBRATE_ORGAN_ANCHOR_SYMBOLS[pname]
        anchor_mask = syms.isin({a.upper() for a in anchor_syms})
        cand = organ_med[anchor_mask & (organ_med >= 5.0)]
        if cand.empty:
            other_samps = [
                sid for sid in values.columns if pname not in native_map.get(str(sid), frozenset())
            ]
            if other_samps:
                off_p95 = values.loc[:, other_samps].quantile(0.95, axis=1).fillna(0.0)
                emp = organ_med[
                    (organ_med >= 50.0) & (organ_med >= 20.0 * np.maximum(off_p95, 1.0))
                ]
                cand = emp.sort_values(ascending=False).head(8)
        if cand.empty:
            continue
        panels[pname] = cand.astype(float)
        prof_cols[pname] = organ_med.fillna(0.0).clip(lower=0.0).astype(float)
    prof_df = (
        pd.DataFrame(prof_cols, index=values.index).fillna(0.0)
        if prof_cols
        else pd.DataFrame(index=values.index)
    )
    return panels, prof_df, native_map


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _numbers(frame):
    """Blank/NA measurements stay missing; malformed or negative values fail."""
    values = frame.replace({"": np.nan, "NA": np.nan, "NaN": np.nan}).apply(
        pd.to_numeric, errors="raise"
    )
    array = values.to_numpy(dtype=float)
    if np.isinf(array).any() or (array < 0).any():
        raise ValueError("Expression values must be finite nonnegative numbers or missing")
    return values.astype(float)


class ExpressionDataset:
    """One native feature-by-library matrix; accessors return defensive copies.

    Source byte pins preserve the raw records. Metadata retains source IDs, tissue
    labels, all mapping alternatives, pooling and unresolved donor/reference facts.
    """

    def __init__(self, values, features, samples, contract, *, records=None):
        if not isinstance(contract, ExpressionContract):
            raise ValueError("ExpressionContract required")
        self._contract = contract
        if records is not None and not isinstance(records, ExpressionRecords):
            raise TypeError("records must be ExpressionRecords or None")
        self._records = records
        self._values = _numbers(values.copy())
        self._features = pd.DataFrame(
            copy.deepcopy(features.to_dict("list")), index=features.index.copy()
        )
        self._samples = pd.DataFrame(
            copy.deepcopy(samples.to_dict("list")), index=samples.index.copy()
        )
        for index in (values.index, values.columns):
            if not index.is_unique or any(not isinstance(x, str) or not x.strip() for x in index):
                raise ValueError("Unique nonempty source feature/sample identifiers required")
        if not features.index.equals(values.index) or not samples.index.equals(values.columns):
            raise ValueError("Matrix feature/sample IDs must exactly match their metadata")
        required = {
            "taxon_id",
            "library_id",
            "run_ids",
            "donor_id",
            "donor_scope",
            "pooling",
            "tissue",
        }
        if not required <= set(samples):
            raise ValueError(f"Missing sample metadata: {sorted(required - set(samples))}")
        if not samples.taxon_id.map(lambda x: type(x) is int and x == contract.taxon_id).all():
            raise ValueError("Matrix sample taxon mismatch")
        if not samples.pooling.isin(["single", "pooled", "unknown"]).all():
            raise ValueError("Explicit single/pooled/unknown sample pooling required")
        seen_runs = set()
        for runs in samples.run_ids:
            if not isinstance(runs, list) or any(
                not isinstance(run, str) or not run for run in runs
            ):
                raise ValueError("run_ids must be explicit lists (empty means unresolved)")
            if len(set(runs)) != len(runs) or seen_runs.intersection(runs):
                raise ValueError(
                    "Run reused across source libraries; reconcile lineage before summarizing"
                )
            seen_runs.update(runs)

    @property
    def contract(self):
        return self._contract

    @property
    def values(self):
        return self._values.copy()

    @property
    def features(self):
        return pd.DataFrame(
            copy.deepcopy(self._features.to_dict("list")), index=self._features.index.copy()
        )

    @property
    def samples(self):
        return pd.DataFrame(
            copy.deepcopy(self._samples.to_dict("list")), index=self._samples.index.copy()
        )

    @property
    def provenance(self):
        """Native dataset identity and processing, or None for a matrix-only snapshot."""
        return self._records.provenance if self._records is not None else None

    def native_records(self, kind):
        """Stream source records with their original schema, IDs and uncertainty.

        Kinds: samples, assays, measurements, reference_mappings, mappings,
        relations, pairings. An empty iterator does not establish biological absence.
        Historical matrix-only snapshots have no native record archive.
        """
        if kind not in RECORD_KINDS:
            raise ValueError(f"Unknown expression record kind {kind!r}")
        if self._records is not None:
            yield from self._records.records(kind)

    def evidence_sources(self):
        """Inspect optional assignment sources without importing or downloading them."""
        from .expression_evidence import evidence_sources

        return evidence_sources(self)

    def with_evidence(
        self,
        source_id,
        *,
        feature_ids,
        offline=True,
        cache_root=None,
        local_overrides=None,
        max_bytes=0,
        max_records=50000,
    ):
        """Attach assignments for selected native features, retaining every candidate.

        Values, missingness and the measurement contract are unchanged. Evidence
        is same-study context, never another observation or isoform allocation.
        Acquisition requires explicit opt-in and a sufficient transfer budget.
        """
        from .expression_evidence import with_evidence

        return with_evidence(
            self,
            source_id,
            feature_ids=feature_ids,
            offline=offline,
            cache_root=cache_root,
            local_overrides=local_overrides,
            max_bytes=max_bytes,
            max_records=max_records,
        )

    def tissue_annotations(self):
        """Cited current tissue-scope annotations; retain the frozen sample metadata.

        No matrix download, feature remapping or cross-species expression inference.
        Unknown source versions return unresolved mappings.
        """
        from .normal_tissues import resolve_normal_tissue

        def optional(value):
            return None if value is None or pd.isna(value) else str(value)

        records = [
            resolve_normal_tissue(
                optional(row["tissue"]),
                species=self.contract.taxon_id,
                cell_type=optional(row.get("cell_type")),
                stage=optional(row.get("stage")),
                sex=optional(row.get("sex")),
                physiological_state=optional(row.get("physiological_state")),
                source_id=self.contract.source_id,
                source_version=self.contract.source_version,
                context=self.contract.context,
            ).as_dict()
            for row in self._samples.to_dict("records")
        ]
        return pd.DataFrame(records, index=self._samples.index.copy())

    def process(
        self,
        *,
        qc=None,
        cleanup=None,
        normalize: str | None = None,
    ) -> ExpressionDataset:
        """Apply unified ``(qc=, cleanup=, normalize=)`` bulk operations in native feature space (#780).

        - ``qc``: ``("pooled_libraries", "broad_background", "swaps", "splash")``, ``"all"``, or ``None``.
          Rejects pooled/corrupted/swapped libraries or ``NaN``-masks splashed ``(feature, sample)``
          cells without altering surviving observed values.
        - ``cleanup``: ``("index_hopping", "deconvolution")``, ``"all"``, or ``None``.
          Subtracts stoichiometric lane bleed and off-lineage necropsy shadow derived strictly from
          within-dataset species-native organ medians (protecting reproductive programs when
          ``contract.context == 'tumor'``).
        - ``normalize``: ``None`` (native unit) or ``"clean_tpm"`` (species-native 16/9/75 scale).
        """
        qc_set = bulk_qc.expand_ops(
            qc,
            allowed=ALL_SPECIES_BULK_QC,
            aliases=SPECIES_BULK_QC_ALIASES,
            param_name="qc",
            domain_label="species bulk",
        )
        cleanup_set = bulk_qc.expand_ops(
            cleanup,
            allowed=ALL_SPECIES_BULK_CLEANUP,
            aliases=SPECIES_BULK_CLEANUP_ALIASES,
            param_name="cleanup",
            domain_label="species bulk",
        )
        norm = bulk_qc.resolve_normalize_option(normalize, domain_label="species bulk")
        if not qc_set and not cleanup_set and norm is None:
            return ExpressionDataset(
                self.values, self.features, self.samples, self.contract, records=self._records
            )

        rna_linear_ops = bool(
            cleanup_set or norm == "clean_tpm" or (qc_set & {"broad_background", "swaps", "splash"})
        )
        if rna_linear_ops and (
            self._contract.measurement.feature_level not in {"gene", "transcript"}
            or self._contract.measurement.scale != "linear"
            or self._contract.measurement.unit
            not in {
                "TPM",
                "TMM-CPM",
                "CPM",
                "FPKM",
            }
        ):
            raise ValueError(
                f"Bulk RNA processing requires raw linear gene/transcript RNA measurement semantics; "
                f"got feature_level={self._contract.measurement.feature_level!r}, "
                f"unit={self._contract.measurement.unit!r}, scale={self._contract.measurement.scale!r}"
            )
        if norm == "clean_tpm" and self._contract.measurement.unit not in {
            "TPM",
            "TMM-CPM",
            "CPM",
        }:
            raise ValueError(
                f"normalize='clean_tpm' requires native linear TPM/CPM measurement semantics; "
                f"got unit={self._contract.measurement.unit!r}"
            )

        vals = self.values
        feats = self.features
        samps = self.samples
        syms = _feature_symbols(feats, vals.index)
        prot_panels = (
            bulk_qc.REPRODUCTIVE_ORGAN_PANELS if self._contract.context == "tumor" else frozenset()
        )

        keep_samps = list(samps.index)
        if "pooled_libraries" in qc_set:
            keep_samps = [sid for sid in keep_samps if samps.loc[sid, "pooling"] != "pooled"]
        if "broad_background" in qc_set and keep_samps:
            or_mask = syms.map(lambda s: bool(_SPECIES_OR_RE.match(s.strip()))).to_numpy()
            if or_mask.any():
                or_sub = vals.loc[or_mask, keep_samps]
                meas = or_sub.notna().sum(axis=0)
                share = or_sub.ge(5.0).sum(axis=0) / meas.replace(0, np.nan)
                keep_samps = [
                    sid
                    for sid in keep_samps
                    if not (int(meas.loc[sid]) >= 20 and float(share.loc[sid]) > 0.20)
                ]

        panels, prof_df, native_map = _species_native_organ_panels(
            vals.loc[:, keep_samps], feats, samps.loc[keep_samps]
        )

        if "swaps" in qc_set and keep_samps and panels:
            basis_genes = sorted({g for p in panels for g in panels[p].index if g in vals.index})
            if basis_genes:
                basis = prof_df.reindex(index=basis_genes, columns=sorted(panels)).fillna(0.0)
                w_mat = bulk_qc.solve_nonnegative_archetype_mixture(
                    basis, vals.loc[basis_genes, keep_samps]
                )
                surviving: list[str] = []
                for sid in keep_samps:
                    nat_p = native_map.get(sid, frozenset())
                    foreign_p = [p for p in w_mat.index if p not in nat_p and p not in prot_panels]
                    f_max = float(w_mat.loc[foreign_p, sid].max()) if foreign_p else 0.0
                    n_max = (
                        float(w_mat.loc[list(nat_p & set(w_mat.index)), sid].max())
                        if (nat_p & set(w_mat.index))
                        else 0.0
                    )
                    if f_max >= 0.35 and f_max >= 2.0 * max(n_max, 1e-6):
                        continue
                    surviving.append(sid)
                keep_samps = surviving
                panels, prof_df, native_map = _species_native_organ_panels(
                    vals.loc[:, keep_samps], feats, samps.loc[keep_samps]
                )

        vals = vals.loc[:, keep_samps]
        samps = samps.loc[keep_samps]
        sub_native = {sid: native_map.get(sid, frozenset()) for sid in keep_samps}

        splash_modified = False
        if "splash" in qc_set and not vals.empty and panels:
            pre_na = vals.isna().to_numpy()
            vals = bulk_qc.mask_panel_splash_cells(
                vals,
                panels,
                sample_native_panels=sub_native,
                protected_panels=prot_panels,
                full_profiles=prof_df,
            )
            splash_modified = not np.array_equal(pre_na, vals.isna().to_numpy())

        transforms = list(self._contract.measurement.transformations)
        curr_unit = str(self._contract.measurement.unit)
        n_feat = len(vals.index)

        if cleanup_set and not vals.empty and panels:
            if "index_hopping" in cleanup_set:
                vals = bulk_qc.subtract_panel_index_hopping(
                    vals,
                    panels,
                    sample_native_panels=sub_native,
                    protected_panels=prot_panels,
                    full_profiles=prof_df,
                )
            if "deconvolution" in cleanup_set:
                vals = bulk_qc.deconvolve_panel_admixture(
                    vals,
                    panels,
                    sample_native_panels=sub_native,
                    protected_panels=prot_panels,
                    full_profiles=prof_df,
                )
            transforms.append(
                MeasurementTransform(
                    operation="bulk_cleanup:" + "+".join(sorted(cleanup_set)),
                    input_unit=curr_unit,
                    output_unit=curr_unit,
                    assumptions=(
                        "within-dataset species-native organ anchor panels",
                        "non-negative residual clamp preserving missingness",
                    ),
                    input_features=n_feat,
                    output_features=n_feat,
                )
            )

        is_proxy = self._contract.measurement.is_proxy
        unfilled_attrs = None
        if norm == "clean_tpm":
            if not vals.empty:
                vals = bulk_qc.species_compartment_clean_scale(vals, symbols=syms)
                unfilled_attrs = vals.attrs.get("unfilled_compartments")
            else:
                # Validate identifiable compartments even when 0 columns survive QC
                bulk_qc.species_compartment_clean_scale(
                    pd.DataFrame({"_check": np.ones(len(vals.index))}, index=vals.index),
                    symbols=syms,
                )
            out_unit = "clean-TPM" if "TPM" in curr_unit else "clean-CPM"
            transforms.append(
                MeasurementTransform(
                    operation="species_compartment_clean_scale",
                    input_unit=curr_unit,
                    output_unit=out_unit,
                    assumptions=(
                        "16% ribosomal / 9% technical / 75% biological compartment scale",
                        "conserved symbol-classified compartments without human ID mapping",
                    ),
                    input_features=n_feat,
                    output_features=n_feat,
                )
            )
            curr_unit = out_unit

        new_contract = self._contract
        if transforms != list(self._contract.measurement.transformations):
            new_meas = ExpressionMeasurement(
                feature_level=self._contract.measurement.feature_level,
                unit=curr_unit,
                scale=self._contract.measurement.scale,
                is_proxy=is_proxy,
                transformations=tuple(transforms),
            )
            new_contract = replace(self._contract, measurement=new_meas)

        untouched_records = (
            self._records
            if (
                keep_samps == list(self._samples.index)
                and not splash_modified
                and not cleanup_set
                and norm is None
            )
            else None
        )
        ds = ExpressionDataset(vals, feats, samps, new_contract, records=untouched_records)
        if unfilled_attrs:
            ds._values.attrs["unfilled_compartments"] = dict(unfilled_attrs)
        return ds

    def tissue_summary(self, *, qc=None, cleanup=None, normalize: str | None = None):
        """Within-source native values; measured libraries and donor coverage remain separate.

        Optional ``(qc=, cleanup=, normalize=)`` arguments run :meth:`process` prior to
        summarizing (#780). When omitted, summarizes the untouched native dataset.
        """
        if qc is not None or cleanup is not None or normalize is not None:
            return self.process(qc=qc, cleanup=cleanup, normalize=normalize).tissue_summary()
        rows = []
        group_columns = ["tissue", "stage"] if "stage" in self._samples else ["tissue"]
        for keys, metadata in self._samples.groupby(group_columns, dropna=False, sort=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            context = {
                col: None if pd.isna(value) else value for col, value in zip(group_columns, keys)
            }
            values = self._values.loc[:, metadata.index]
            measured = values.notna()
            individual = metadata[metadata.pooling.eq("single") & metadata.donor_id.notna()]
            # Count each known donor once per feature, regardless of technical replicates.
            donor_counts = pd.Series(0, index=values.index, dtype=int)
            for _, donor_samples in individual.groupby(["donor_scope", "donor_id"], dropna=False):
                donor_counts += measured.loc[:, donor_samples.index].any(axis=1).astype(int)
            summary = pd.DataFrame(
                {
                    "feature_id": values.index,
                    **context,
                    "native_mean": values.mean(axis=1).to_numpy(),
                    "native_max": values.max(axis=1).to_numpy(),
                    "measured_libraries": measured.sum(axis=1).to_numpy(),
                    "total_libraries": len(metadata),
                    "known_individual_donors_measured": donor_counts.to_numpy(),
                    "libraries_with_unresolved_or_pooled_donors": len(metadata) - len(individual),
                    "unit": self.contract.measurement.unit,
                }
            )
            rows.append(summary)
        empty_cols = [
            "feature_id",
            *group_columns,
            "native_mean",
            "native_max",
            "measured_libraries",
            "total_libraries",
            "known_individual_donors_measured",
            "libraries_with_unresolved_or_pooled_donors",
            "unit",
        ]
        result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=empty_cols)
        result.attrs["expression_contract"] = self.contract.as_dict()
        return result

    def save(self, destination):
        """Write a new portable artifact atomically; refuse to replace existing results."""
        destination = Path(destination)
        if destination.exists():
            raise FileExistsError(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".oncoref-expression-", dir=destination.parent
        ) as temporary:
            staged = Path(temporary) / "artifact"
            staged.mkdir()
            self._values.to_parquet(staged / "values.parquet")
            version = 2 if self._records is not None else 1
            if version == 1:
                self._features.to_parquet(staged / "features.parquet")
                feature_file = "features.parquet"
            else:
                feature_file = "features.json"
                # Keep the index separate from arbitrary native feature columns.
                features = self._features.to_dict("tight")
                (staged / feature_file).write_text(json.dumps(features, allow_nan=False) + "\n")
            # JSON retains nested original sample facts and provenance without stringifying them.
            records = copy.deepcopy(
                self._samples.astype(object).where(self._samples.notna(), None).to_dict("records")
            )
            for sample_id, row in zip(self._samples.index, records):
                row["sample_id"] = sample_id
            (staged / "samples.json").write_text(json.dumps(records, allow_nan=False) + "\n")
            names = ["values.parquet", feature_file, "samples.json"]
            if self._records is not None:
                names.extend(self._records.write(staged))
            payload = {
                "schema": "oncoref.expression-dataset",
                "schema_version": version,
                "contract": self.contract.as_dict(),
                "files": {
                    name: {
                        "sha256": _sha256(staged / name),
                        "size_bytes": (staged / name).stat().st_size,
                    }
                    for name in names
                },
            }
            if version == 2:
                payload["provenance"] = self.provenance
            record = {**payload, "artifact_id": _identity("expression-dataset", version, payload)}
            (staged / "metadata.json").write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n"
            )
            os.rename(staged, destination)
        return destination

    @classmethod
    def load(cls, directory):
        """Offline replay validates the metadata identity and every artifact file."""
        directory = Path(directory)
        record = json.loads((directory / "metadata.json").read_text())
        version = record.get("schema_version")
        fields = {"schema", "schema_version", "contract", "files", "artifact_id"}
        if version == 2:
            fields.add("provenance")
        if (
            set(record) != fields
            or record["schema"] != "oncoref.expression-dataset"
            or type(record["schema_version"]) is not int
            or version not in (1, 2)
        ):
            raise ValueError("Unsupported expression artifact schema")
        payload = {key: value for key, value in record.items() if key != "artifact_id"}
        if record["artifact_id"] != _identity("expression-dataset", version, payload):
            raise ValueError("Expression artifact metadata digest mismatch")
        expected = {"values.parquet", "features.parquet", "samples.json"}
        if version == 2:
            expected = {"values.parquet", "features.json", "samples.json"}
            expected.update(f"{kind}.jsonl.gz" for kind in RECORD_KINDS)
        if set(record["files"]) != expected:
            raise ValueError("Unexpected expression artifact file inventory")
        for name, info in record["files"].items():
            path = directory / name
            if (
                set(info) != {"sha256", "size_bytes"}
                or path.stat().st_size != info["size_bytes"]
                or _sha256(path) != info["sha256"]
            ):
                raise ValueError(f"Expression artifact checksum mismatch: {name}")
        raw_samples = json.loads((directory / "samples.json").read_text())
        if raw_samples:
            samples = pd.DataFrame(raw_samples).set_index("sample_id")
        else:
            samples = pd.DataFrame(
                columns=[
                    "sample_id",
                    "taxon_id",
                    "library_id",
                    "run_ids",
                    "donor_id",
                    "donor_scope",
                    "pooling",
                    "tissue",
                ]
            ).set_index("sample_id")
        return cls(
            pd.read_parquet(directory / "values.parquet"),
            pd.read_parquet(directory / "features.parquet")
            if version == 1
            else pd.DataFrame.from_dict(
                json.loads((directory / "features.json").read_text()), orient="tight"
            ),
            samples,
            ExpressionContract.from_dict(record["contract"]),
            records=ExpressionRecords.from_directory(directory, record["provenance"])
            if version == 2
            else None,
        )


PROTEIN_SAMPLE_COLUMNS = (
    "Run",
    "Channel",
    "BioReplicate",
    "Condition",
    "TechRepMixture",
    "Mixture",
    "sample.type",
)


def protein_sample_id(record):
    return _identity("protein-sample", 1, {key: record[key] for key in PROTEIN_SAMPLE_COLUMNS})


def read_species_expression(
    source_id,
    *,
    species,
    offline=False,
    cache_root=None,
    local_overrides=None,
    max_bytes=DEFAULT_MAX_BYTES,
    qc=None,
    cleanup=None,
    normalize: str | None = None,
):
    """Load a registered source in its native space; optional ``(qc=, cleanup=, normalize=)`` (#780)."""
    qc_set = bulk_qc.expand_ops(
        qc,
        allowed=ALL_SPECIES_BULK_QC,
        aliases=SPECIES_BULK_QC_ALIASES,
        param_name="qc",
        domain_label="species bulk",
    )
    cleanup_set = bulk_qc.expand_ops(
        cleanup,
        allowed=ALL_SPECIES_BULK_CLEANUP,
        aliases=SPECIES_BULK_CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="species bulk",
    )
    norm = bulk_qc.resolve_normalize_option(normalize, domain_label="species bulk")

    selected = resolve_species(species)
    source = expression_source(source_id, species=selected)
    if source is None:
        raise KeyError(f"Unknown {selected.name} expression source {source_id!r}")
    entry = next(
        row
        for row in expression_source_registry_entries(species=selected)
        if row["id"] == source_id
    )
    adapter = entry.get("adapter")
    if adapter not in {"bgee", "barkbase", "orfdata"}:
        raise ValueError(
            f"{selected.name}/{source_id} is unavailable ({source.acquisition_status or 'reader_not_registered'}): {source.notes or 'no native reader registered'}"
        )
    contract = source.expression_contract
    expected = {
        "bgee": ("gene", "TPM"),
        "barkbase": ("gene", "TMM-CPM"),
    }.get(adapter, (contract.measurement.feature_level, contract.measurement.unit))
    if (
        (contract.measurement.feature_level, contract.measurement.unit) != expected
        or contract.measurement.scale != "linear"
        or contract.measurement.is_proxy is not False
    ):
        raise ValueError("Native reader and declared measurement semantics disagree")
    filename = entry["sample_metadata"]
    if Path(filename).name != filename:
        raise ValueError("Sample metadata must be a package-local filename")
    metadata_path = expression_source_registry_path(species=selected).parent / filename
    if _sha256(metadata_path) != contract.content_provenance["samples"]:
        raise ValueError("Sample metadata checksum mismatch")
    metadata = json.loads(metadata_path.read_text())
    if (
        set(metadata) != {"schema", "schema_version", "taxon_id", "source_id", "samples"}
        or metadata["schema"] != "oncoref.expression-samples"
        or type(metadata["schema_version"]) is not int
        or metadata["schema_version"] != 1
        or metadata["taxon_id"] != selected.taxon_id
        or metadata["source_id"] != source_id
    ):
        raise ValueError("Expression sample metadata source/schema mismatch")
    samples = pd.DataFrame(metadata["samples"]).set_index("sample_id")
    if adapter == "orfdata":
        from .managed_expression import read_managed_expression

        ds = read_managed_expression(
            entry,
            samples,
            contract,
            offline=offline,
            cache_root=cache_root,
            local_overrides=local_overrides,
            max_bytes=max_bytes,
        )
        if qc_set or cleanup_set or norm is not None:
            return ds.process(qc=qc_set, cleanup=cleanup_set, normalize=norm)
        return ds
    asset = species_asset(entry["matrix_asset"], species=selected)
    if asset.sha256 != contract.content_provenance["matrix"]:
        raise ValueError("Matrix asset and expression contract disagree")
    path = ensure_asset(
        asset,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )
    if adapter == "bgee":
        chunks = read_csv(
            path,
            sep="\t",
            usecols=["Experiment ID", "Library ID", "Gene ID", "TPM"],
            dtype=str,
            keep_default_na=False,
            chunksize=100000,
        )
        frame = pd.concat(chunks, ignore_index=True)
        if set(frame["Experiment ID"]) != {source.accession}:
            raise ValueError("Bgee matrix contains an unexpected study")
        if frame.duplicated(["Gene ID", "Library ID"]).any():
            raise ValueError("Duplicate Bgee gene/library measurements")
        frame["TPM"] = _numbers(frame[["TPM"]])["TPM"]
        values = frame.pivot(index="Gene ID", columns="Library ID", values="TPM")
        features = pd.DataFrame({"source_feature_id": values.index}, index=values.index)
    else:
        frame = read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        required = {"StringTieID", "EnsemblGeneID", "GeneSymbol"}
        if not required <= set(frame):
            raise ValueError("BarkBase native feature columns are missing")
        frame = frame.set_index("StringTieID")
        # StringTie identity is primary: retain unmapped and multi-mapped feature rows.
        features = frame[["EnsemblGeneID", "GeneSymbol"]].copy()
        features.insert(0, "source_feature_id", features.index)
        values = _numbers(frame.drop(columns=["EnsemblGeneID", "GeneSymbol"]))
    if set(samples.index) != set(values.columns):
        raise ValueError("Quantitative libraries and sample metadata disagree")
    samples = samples.loc[values.columns]
    values.index.name = features.index.name = "feature_id"
    values.columns.name = samples.index.name = "sample_id"
    ds = ExpressionDataset(values, features, samples, contract)
    if qc_set or cleanup_set or norm is not None:
        return ds.process(qc=qc_set, cleanup=cleanup_set, normalize=norm)
    return ds


def combine_expression(datasets):
    """Concatenate compatible native libraries without zero filling or double counting.

    Unknown run lineage cannot establish independent libraries. Source contracts
    must already carry verified reference identities; this is not a crosswalk.
    """
    datasets = tuple(datasets)
    if not datasets:
        raise ValueError("At least one expression dataset required")
    first = datasets[0]
    seen = set()
    for dataset in datasets:
        validate_expression_join(first.contract, dataset.contract, combine_values=True)
        for runs in dataset._samples.run_ids:
            if not runs:
                raise ValueError("Cannot combine libraries with unresolved run lineage")
            if seen.intersection(runs):
                raise ValueError("Shared runs would count reused libraries more than once")
            seen.update(runs)
    # Returning component contracts retains source provenance for all contributors.
    values = pd.concat([d._values for d in datasets], axis=1, verify_integrity=True)
    return values, tuple(d.contract.as_dict() for d in datasets)


def expression_lineage_overlaps(species):
    """Source inventory overlap evidence; independent-donor counts remain unresolved."""
    path = expression_source_registry_path(species=species).parent / "lineage_overlaps.json"
    return json.loads(path.read_text()) if path.exists() else None
