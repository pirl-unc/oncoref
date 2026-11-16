#!/usr/bin/env python3
"""Review the tenfold cancer rule: fold thresholds, near misses, sample background, CTpedia.

Reads the released catalog, the rebuild's per-sample tumor matrix and normal panel, and
the pinned tumor matrices. Writes ``docs/audits/cta-threshold-review-20261116/`` and
``docs/ctpedia-filtering.md``. It changes no catalog decision: the background check is
a proposal, and the threshold sweep reclassifies a copy of the catalog.
"""

from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import rebuild_target_catalog as rb  # noqa: E402
from build_target_filter_audit import FAMILIAR  # noqa: E402
from cta_report_common import sha256, write_csv  # noqa: E402

from oncoref import gene_knowledge, target_filters, target_selection  # noqa: E402

OUT = ROOT / "docs/audits/cta-threshold-review-20261116"
CTPEDIA_PAGE = ROOT / "docs/ctpedia-filtering.md"
REBUILD = ROOT / "docs/audits/target-selection-20261110/rebuild"
TUMOR_IDENTITY = ROOT / "docs/audits/cta-selection-20261028/tumor-identity-repair"
THRESHOLDS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 100]
#: Background check (proposed, not applied): share of olfactory receptor genes at
#: >= 5 TPM in one tumor library, evaluated only when >= 100 receptors are measured.
OR_PATTERN = re.compile(r"^OR\d+[A-Z]+\d+$")
BACKGROUND_TPM = 5.0
BACKGROUND_SHARE = 0.20
BACKGROUND_MIN_GENES = 100
BACKGROUND_SWEEP = [0.20, 0.10, 0.05, 0.02]
NEAR_MISS = (8.0, 10.0)
CHECK_LABELS = {
    "normal_screen_pass": "normal-screen hold",
    "required_coverage_complete": "required normal coverage incomplete",
    "identity_translation_pass": "identity or translation not qualified",
    "evidence_resolved": "evidence unresolved",
    "expression_bounds_pass": "expression bounds fail",
    "sequence_pass": "sequence",
}
LAYER_LABELS = {
    "off_target": "off-target TPM",
    "curation": "curated tier",
    "sets": "set",
    "papers": "papers",
    "exclusions": "name exclusion",
}
SWEEP_HEADINGS = {
    "fold": "Fold",
    "cta": "CTA",
    "cancer_antigens": "Cancer Antigens",
    "familiar_cta": "Familiar CTA",
    "ctpedia_cta": "CTpedia CTA",
    "median_best_tumor_tpm": "Median best tumor TPM",
    "pct_one_sample": "% one sample",
    "pct_flagged_sample_only": "% flagged library only",
    "pct_two_cancer_types": "% ≥2 cancer types",
    "pct_published_nomination": "% nominated",
    "pct_mapping_risk": "% mapping risk",
}
NEAR_HEADINGS = {
    "protein": "Protein",
    "group": "Group",
    "tumor_fold": "Fold",
    "best_tumor_tpm": "Best tumor TPM",
    "best_tumor_cancer": "Cancer",
    "samples_ge_5x": "Samples ≥5x",
    "cancer_types_ge_5x": "Cancer types ≥5x",
    "normal_ceiling_tpm": "Normal ceiling TPM",
    "normal_ceiling_context": "Ceiling context",
    "datasets_at_half_ceiling": "Datasets ≥ half ceiling",
    "fold_without_that_dataset": "Fold without that dataset",
    "mapping_risk": "Mapping risk",
}
IMPACT_HEADINGS = {
    "olfactory_share_above": "Drop libraries above",
    "samples_dropped": "Libraries dropped",
    "cohorts": "Cohorts",
    "cancer_antigens_lost": "Cancer Antigens lost",
    "of_which_cta": "of which CTA",
    "proteins": "Proteins",
}
COLORS = {
    "surface": "#fcfcfb",
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "grid": "#e4e3de",
    "ref": "#8a8984",
    "s1": "#2a78d6",
    "s2": "#eb6834",
}


def dataset_family(source_id: str) -> str:
    source = source_id.split(":")[0]
    if source.startswith("dice_"):
        return "DICE immune"
    return {
        "hpa_immune_monaco": "Monaco immune",
        "rna_tissue_hpa": "HPA tissue",
        "rna_tissue_gtex": "GTEx",
        "gtex_tissue_rna": "GTEx",
        "rna_brain_gtex": "GTEx",
        "rna_brain_hpa": "HPA brain",
    }[source]


def fmt(value, digits=1) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return ""
    return f"{value:,.{digits}f}"


def flag(value) -> bool:
    """A boolean catalog column read from CSV: True, False or missing."""
    return bool(value) if pd.notna(value) else False


def markdown_table(frame: pd.DataFrame) -> str:
    cells = frame.astype(object).where(frame.notna(), "").map(str)
    cells = cells.map(lambda text: text.replace("|", "/"))
    lines = ["| " + " | ".join(cells.columns) + " |"]
    lines.append("|" + "|".join("---:" if is_num else "---" for is_num in numeric(frame)) + "|")
    lines += ["| " + " | ".join(r) + " |" for r in cells.to_numpy()]
    return "\n".join(lines)


def numeric(frame):
    return [pd.api.types.is_numeric_dtype(frame[c]) for c in frame.columns]


def load():
    catalog = target_selection.results()
    recorded = json.loads((REBUILD / "comparison.json").read_text())["tumor_matrix_cache_sha256"]
    cache = rb.CACHE / "tumor-comparison.parquet"
    if sha256(cache) != recorded:
        raise ValueError("The cached tumor matrix is not the one the released catalog used")
    matrix = pd.read_parquet(cache)
    samples = pd.read_csv(REBUILD / "tumor-samples.csv.gz").drop_duplicates("sample_key")
    cancer_code = samples.set_index("sample_key").cancer_code.reindex(matrix.index)
    if cancer_code.isna().any():
        raise ValueError("A tumor sample has no cancer code")
    normal = pd.read_csv(REBUILD / "normal-comparison.csv.gz")
    normal = normal[normal.qc_pass & normal.identity_resolved & normal.value.notna()].copy()
    normal["family"] = normal.source_id.map(dataset_family)
    return catalog, matrix, cancer_code, normal


def sample_background() -> pd.DataFrame:
    """Olfactory receptor background and detected genes for every pinned tumor library."""
    manifest = json.loads((TUMOR_IDENTITY / "manifest.json").read_text())
    paths = {Path(p).stem: (Path(p), digest) for p, digest in manifest["matrices"].items()}
    available = pd.read_csv(TUMOR_IDENTITY / "cohort-availability.csv")
    identity = pd.read_csv(TUMOR_IDENTITY / "sample-identity.csv.gz")
    rows = []
    for job in available.itertuples():
        path, digest = paths[job.cancer_code]
        if sha256(path) != digest:
            raise ValueError(f"Changed tumor matrix: {path}")
        ids = identity[identity.cancer_code.eq(job.cancer_code)]
        columns = ids.sample_id.astype(str).tolist()
        m = pq.read_table(path).to_pandas()
        values = m[columns].apply(pd.to_numeric, errors="coerce")
        receptors = values[m.Symbol.astype(str).str.match(OR_PATTERN)]
        high = receptors.ge(BACKGROUND_TPM)
        alone = high.sum(axis=1).eq(1)
        for sample in columns:
            n_high = int(high[sample].sum())
            rows.append(
                {
                    "sample_key": f"{ids.namespace.iloc[0]}:{sample}",
                    "cancer_code": job.cancer_code,
                    "olfactory_receptors_measured": int(receptors[sample].notna().sum()),
                    "olfactory_share_ge5": float((receptors[sample] >= BACKGROUND_TPM).mean()),
                    "olfactory_share_ge1": float((receptors[sample] >= 1).mean()),
                    # Receptors at 5 TPM here and in no other library of the cohort.
                    "olfactory_private_share": (
                        float((high[sample] & alone).sum() / n_high) if n_high else np.nan
                    ),
                    "genes_detected": int((values[sample] > 0).sum()),
                }
            )
    out = pd.DataFrame(rows).drop_duplicates("sample_key").reset_index(drop=True)
    out["evaluable"] = out.olfactory_receptors_measured.ge(BACKGROUND_MIN_GENES)
    out["flagged"] = out.evaluable & out.olfactory_share_ge5.gt(BACKGROUND_SHARE)
    return out


def reclassify(catalog: pd.DataFrame, fold: float) -> pd.DataFrame:
    """The catalog's own classify(), with cancer evidence re-read at another fold."""
    f = catalog.drop(
        columns=[
            "category",
            "subtype",
            "classification",
            "is_cancer_antigen",
            "is_reproductive_antigen",
            "is_cta",
            "is_target_candidate",
        ]
    )
    f["unsafe"] = catalog.category.eq("Unsafe")
    for check in target_selection.SHARED_CHECKS:
        f[check] = f[check].astype(bool)
    f["home_tissue_unmeasured"] = f.home_tissue.isin(target_selection.UNMEASURED_HOME_TISSUES)
    assessed = f.cancer_evidence.isin(["qualifies", "does_not_qualify"])
    f.loc[assessed, "cancer_evidence"] = np.where(
        f.loc[assessed, "tumor_fold"].ge(fold), "qualifies", "does_not_qualify"
    )
    return target_selection.classify(f)


def label_lookup(catalog: pd.DataFrame, name: str):
    hit = catalog[catalog.protein_group_label.map(lambda x: name in (x, x.split("/")[0]))]
    return None if hit.empty else hit.entity_id.iloc[0]


def ctpedia_entities(catalog: pd.DataFrame) -> pd.DataFrame:
    table = pd.read_csv(ROOT / "oncoref/data/cancer-testis-antigens.csv")
    table = table[table.source_databases.astype(str).str.contains("CTpedia")]
    gene_to_entity = {
        g: r.entity_id for r in catalog.itertuples() for g in json.loads(r.member_gene_ids)
    }
    table = table.assign(entity_id=table.Ensembl_Gene_ID.map(gene_to_entity))
    if table.entity_id.isna().any():
        raise ValueError("A CTpedia gene is missing from the catalog")
    return (
        table.groupby("entity_id")
        .Symbol.agg(lambda s: ", ".join(sorted(s)))
        .rename("ctpedia_genes")
    )


def sweep(catalog, matrix, cancer_code, flagged, ctpedia):
    entity = catalog.set_index("entity_id")
    base = entity.loc[matrix.columns, "single_tumor_threshold_tpm"] / target_selection.TUMOR_FOLD
    familiar = {
        name: label_lookup(catalog, name) for names in FAMILIAR.values() for name in names.split()
    }
    familiar_ids = {e for e in familiar.values() if e is not None}
    mapping_risk = entity.paralog_ambiguous.eq(True) | entity.host_overlapped.eq(True)
    rows = []
    for fold in THRESHOLDS:
        g = reclassify(catalog, fold).set_index("entity_id")
        cta = g.index[g.is_cta]
        hits = matrix[cta].ge(base[cta] * fold)
        samples = hits.sum()
        cancers = hits.apply(lambda col: cancer_code[col].nunique())
        clean = hits.loc[~matrix.index.isin(flagged)].any()
        rows.append(
            {
                "fold": fold,
                "cta": len(cta),
                "cancer_antigens": int(g.is_cancer_antigen.sum()),
                "familiar_cta": len(familiar_ids & set(cta)),
                "familiar_named": len(familiar),
                "ctpedia_cta": int(g.loc[ctpedia.index, "is_cta"].sum()),
                "ctpedia_proteins": len(ctpedia),
                "pct_one_sample": 100 * samples.eq(1).mean(),
                "pct_two_cancer_types": 100 * cancers.ge(2).mean(),
                "pct_flagged_sample_only": 100 * (~clean).mean(),
                "median_best_tumor_tpm": float(matrix[cta].max().median()),
                "pct_published_nomination": 100 * entity.loc[cta, "published_nomination"].mean(),
                "pct_mapping_risk": 100 * mapping_risk[cta].mean(),
            }
        )
    return pd.DataFrame(rows)


def sweep_figure(table: pd.DataFrame, path: Path):
    panels = [
        ("Proteins in the set", [("cta", "CTA"), ("cancer_antigens", "Cancer Antigens")], ""),
        (
            f"Familiar CTAs that are CTA (of {table.familiar_named.iloc[0]} named)",
            [("familiar_cta", None)],
            "",
        ),
        (
            f"CTpedia proteins that are CTA (of {table.ctpedia_proteins.iloc[0]})",
            [("ctpedia_cta", None)],
            "",
        ),
        ("Median best tumor TPM", [("median_best_tumor_tpm", None)], ""),
        ("% resting on one tumor sample", [("pct_one_sample", None)], "%"),
        ("% with qualifying tumors in ≥2 cancer types", [("pct_two_cancer_types", None)], "%"),
        ("% with a published nomination", [("pct_published_nomination", None)], "%"),
        ("% with paralog or host-gene mapping risk", [("pct_mapping_risk", None)], "%"),
    ]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.edgecolor": COLORS["grid"],
            "axes.labelcolor": COLORS["ink2"],
            "xtick.color": COLORS["ink2"],
            "ytick.color": COLORS["ink2"],
        }
    )
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.4), facecolor=COLORS["surface"])
    x = table.fold.to_numpy()
    at10 = table.fold.eq(target_selection.TUMOR_FOLD)
    for ax, (title, series, unit) in zip(axes.flat, panels):
        ax.set_facecolor(COLORS["surface"])
        ax.axvline(10, color=COLORS["ref"], lw=1, ls=(0, (3, 3)), zorder=1)
        for i, (column, name) in enumerate(series):
            color = COLORS["s1"] if i == 0 else COLORS["s2"]
            y = table[column].to_numpy(dtype=float)
            ax.plot(x, y, color=color, lw=2, marker="o", ms=4.5, zorder=3, label=name)
            value = float(table.loc[at10, column].iloc[0])
            ax.plot([10], [value], "o", ms=7, color=color, mec=COLORS["surface"], mew=1.5, zorder=4)
            text = f"{value:.0f}{unit}"
            if name:
                text = f"{name} {text}"
            # Two series: the lower one is labelled below-left of its point, clear of
            # both lines; the upper one above-right.
            below = len(series) > 1 and i == 0
            ax.annotate(
                text,
                (10, value),
                xytext=(-6, -16) if below else (6, 6),
                textcoords="offset points",
                ha="right" if below else "left",
                fontsize=8.5,
                color=COLORS["ink"],
            )
        ax.set_xscale("log")
        ax.set_xticks([2, 3, 5, 10, 20, 50, 100])
        ax.set_xticklabels(["2x", "3x", "5x", "10x", "20x", "50x", "100x"])
        ax.minorticks_off()
        top = max(float(table[c].max()) for c, _ in series)
        ax.set_ylim(0, top * 1.18)
        ax.grid(axis="y", color=COLORS["grid"], lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_title(title, loc="left", fontsize=9.5, color=COLORS["ink"])
        if len(series) > 1:
            ax.legend(frameon=False, fontsize=8, loc="lower left", labelcolor=COLORS["ink2"])
    fig.suptitle(
        "The CTA set as the tumor-fold threshold moves (dashed line: the catalog's 10x)",
        x=0.01,
        ha="left",
        fontsize=11,
        color=COLORS["ink"],
    )
    fig.supxlabel(
        "Best tumor sample / max(non-reproductive normal, Disposable, 1 TPM)",
        fontsize=9,
        color=COLORS["ink2"],
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=COLORS["surface"])
    plt.close(fig)


def near_misses(catalog, matrix, cancer_code, normal, background, families):
    """Reproductive Antigens whose best tumor sample is 8-10x: CTAs at a lower fold."""
    entity = catalog.set_index("entity_id")
    low, high = NEAR_MISS
    chosen = entity[
        entity.is_reproductive_antigen.astype(bool)
        & entity.tumor_fold.ge(low)
        & entity.tumor_fold.lt(high)
    ]
    share = background.set_index("sample_key").olfactory_share_ge5
    flagged = set(background.loc[background.flagged, "sample_key"])
    heme = families[families.str.startswith("heme")].index
    rows = []
    for e, r in chosen.iterrows():
        base = r.single_tumor_threshold_tpm / target_selection.TUMOR_FOLD
        folds = (matrix[e].dropna() / base).sort_values(ascending=False)
        top = folds.index[0]
        ref = normal[normal.entity_id.eq(e) & normal.scope.isin(["normal", "disposable"])]
        ref = ref.sort_values("value", ascending=False)
        ceiling = max(r.normal_reference_tpm, r.disposable_max_tpm)
        if ceiling < target_selection.DETECTION_FLOOR:
            group = "below 10 TPM"
            replicated = ""
            without = np.nan
        else:
            first = ref.iloc[0]
            near = ref[ref.value.ge(0.5 * ceiling)]
            replicated = " + ".join(sorted(near.family.unique()))
            others = ref[ref.family.ne(first.family)].value.max()
            others = 0.0 if pd.isna(others) else others
            without = folds.iloc[0] * base / max(others, target_selection.DETECTION_FLOOR)
            if near.family.nunique() >= 2:
                group = "normal in two datasets"
            elif first.family in ("DICE immune", "Monaco immune"):
                group = "one immune dataset"
            elif first.family == "HPA brain":
                group = "one brain dataset"
            else:
                group = "one other dataset"
        top_ref = ref.iloc[0] if len(ref) else None
        rows.append(
            {
                "protein": r.protein_group_label,
                "classification": r.classification,
                "group": group,
                "tumor_fold": round(float(folds.iloc[0]), 2),
                "best_tumor_tpm": round(float(folds.iloc[0] * base), 1),
                "best_tumor_cancer": cancer_code[top],
                "best_sample_flagged": top in flagged,
                "best_sample_olfactory_share": round(float(share.get(top, np.nan)), 3),
                "second_sample_fold": round(float(folds.iloc[1]), 1) if len(folds) > 1 else np.nan,
                "samples_ge_5x": int(folds.ge(5).sum()),
                "cancer_types_ge_5x": int(cancer_code[folds.index[folds.ge(5)]].nunique()),
                "blood_cancer_vs_immune_cells": bool(
                    cancer_code[top] in heme
                    and top_ref is not None
                    and top_ref.family.endswith("immune")
                ),
                "normal_ceiling_tpm": round(float(ceiling), 2),
                "normal_ceiling_context": (
                    f"{top_ref.family}: {top_ref.source_context}" if top_ref is not None else ""
                ),
                "datasets_at_half_ceiling": replicated,
                "fold_without_that_dataset": round(float(without), 1),
                "patients_above_normal": r.independent_patients_above_normal,
                "mapping_risk": flag(r.paralog_ambiguous) or flag(r.host_overlapped),
                "published_nomination": bool(r.published_nomination),
                "reproductive_routes": r.reproductive_routes,
            }
        )
    return pd.DataFrame(rows).sort_values("tumor_fold", ascending=False).reset_index(drop=True)


def ctpedia_comparison(entities) -> pd.DataFrame:
    """The rebuild's own tumor comparison for every CTpedia protein, germ cell tumors kept.

    Uses rebuild_target_catalog's normal panel and tumor reader unchanged; its outputs
    go to a temporary folder. Germ cell tumors are reported separately, never as evidence.
    """
    facts = rb.shared_facts()
    subset = facts[facts.entity_id.isin(entities)].copy()
    normal, required = rb.normal_panel(subset, write=False)
    real = (rb.OUT, rb.CACHE, target_selection.GERM_CELL_TUMORS)
    with tempfile.TemporaryDirectory() as tmp:
        rb.OUT = rb.CACHE = Path(tmp)
        target_selection.GERM_CELL_TUMORS = ()
        try:
            matrix, meta = rb.tumors(subset)
        finally:
            rb.OUT, rb.CACHE, target_selection.GERM_CELL_TUMORS = real
    meta = meta.set_index("sample_key").loc[matrix.index]
    germ = meta.cancer_code.isin(target_selection.GERM_CELL_TUMORS)
    grouped = dict(iter(normal.groupby("entity_id")))
    rows = []
    for row in subset.itertuples():
        values = matrix.get(row.entity_id, pd.Series(np.nan, index=matrix.index))
        somatic = values[~germ]
        tumor = pd.DataFrame(
            {
                "entity_id": row.entity_id,
                "source_id": meta.namespace[~germ],
                "sample_id": meta.sample_id[~germ],
                "participant_id": meta.participant_id[~germ],
                "value": somatic,
                "comparison_id": rb.COMPARISON_ID,
                "unit": "TPM",
                "qc_pass": somatic.notna(),
                "identity_resolved": somatic.notna(),
            }
        ).reset_index(drop=True)
        reference = grouped.get(row.entity_id, normal.iloc[:0])
        receipt = target_selection.assess_cancer_rna(
            tumor,
            reference,
            required_normal_contexts=required,
            comparison_id=rb.COMPARISON_ID,
            comparability_reference="comparison.json",
            unit="TPM",
        )
        ok = reference[reference.qc_pass & reference.identity_resolved & reference.value.notna()]
        top_normal = ok[ok.scope.eq("normal")].sort_values("value", ascending=False).head(1)
        base = (receipt["single_sample_threshold"] or np.nan) / target_selection.TUMOR_FOLD
        measured = somatic.notna().any()
        rows.append(
            {
                "entity_id": row.entity_id,
                "rule_status": receipt["status"],
                "rule_tumor_max_tpm": somatic.max() if measured else np.nan,
                "rule_tumor_cancer": meta.cancer_code[somatic.idxmax()] if measured else "",
                "rule_fold": somatic.max() / base if measured else np.nan,
                "rule_normal_tpm": receipt["normal_maximum"],
                "rule_normal_context": (
                    " / ".join(
                        dict.fromkeys(
                            [top_normal.tissue.iloc[0], top_normal.source_context.iloc[0]]
                        )
                    )
                    if len(top_normal)
                    else ""
                ),
                "germ_cell_tumor_fold": values[germ].max() / base
                if values[germ].notna().any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows).set_index("entity_id")


def unsafe_evidence(catalog, entities) -> pd.Series:
    """Highest Unsafe-qualifying observation per evidence layer, one line per protein.

    Bulk and cluster layers come from the replayed 2026-10-28 normal-peak table (per
    locus); rare normal-cell RNA and published expression come from the catalog row.
    """
    entity = catalog.set_index("entity_id")
    gene_to_entity = {g: e for e in entities for g in json.loads(entity.at[e, "member_gene_ids"])}
    peaks = pd.read_csv(
        rb.BASE / "normal-peak-evidence.csv.gz",
        usecols=["gene_id", "tissue", "cell_type", "value", "unit", "assay", "unsafe_evidence"],
        low_memory=False,
    )
    peaks = peaks[peaks.unsafe_evidence.astype(str).eq("True")].copy()
    peaks["entity_id"] = peaks.gene_id.map(gene_to_entity)
    peaks = peaks.dropna(subset=["entity_id"]).sort_values("value", ascending=False)
    top = peaks.drop_duplicates(["entity_id", "assay"])
    out = {}
    for e in entities:
        r = entity.loc[e]
        if r.category != "Unsafe":
            continue
        layers = r.Unsafe_reasons.split(" | ")
        parts = []
        for row in top[top.entity_id.eq(e) & top.assay.isin(layers)].itertuples():
            place = row.tissue
            if isinstance(row.cell_type, str) and row.cell_type not in ("", "bulk"):
                place += f" / {row.cell_type}"
            parts.append(f"{row.assay}: {place} {row.value:.1f} {row.unit}")
        if rb.RARE_UNSAFE in layers and isinstance(r.finding_tissue, str):
            parts.append(
                f"{rb.RARE_UNSAFE}: {r.finding_tissue} / {r.finding_top_cell_type}, "
                f"{r.finding_top_multi_molecule_pct:.1f}% of cells with >=2 molecules, "
                f"{r.finding_donors_multi_molecule:.0f} donors"
            )
        if rb.LITERATURE_UNSAFE in layers:
            parts.append(f"{rb.LITERATURE_UNSAFE}: {r.literature_decision_tissue}")
        out[e] = "; ".join(parts)
    return pd.Series(out, dtype=object)


def ctpedia_table(catalog, ctpedia, comparison) -> pd.DataFrame:
    entity = catalog.set_index("entity_id").loc[ctpedia.index]
    evidence = unsafe_evidence(catalog, list(ctpedia.index))
    unsafe_symbols = {
        symbol
        for label in catalog.loc[catalog.category.eq("Unsafe"), "protein_group_label"]
        for symbol in label.split("/")
    }
    knowledge = gene_knowledge.findings().set_index("id")
    modes = {
        name: target_filters.apply(f).set_index("entity_id")
        for name, f in target_filters.FILTER_MODES.items()
        if name != "all"
    }
    rows = []
    for e, r in entity.iterrows():
        c = comparison.loc[e]
        assessed = r.cancer_evidence != "unassessed"
        # The rerun must reproduce the released catalog exactly.
        disagrees = c.rule_status != r.cancer_evidence or not np.isclose(
            c.rule_fold, r.tumor_fold, rtol=1e-6
        )
        if assessed and disagrees:
            raise ValueError(f"Rerun disagrees with the catalog for {r.protein_group_label}")
        ways = []
        unsafe = r.category == "Unsafe"
        if unsafe:
            ways.append(f"Unsafe: {r.Unsafe_reasons.replace(' | ', ', ')}")
        for check, text in CHECK_LABELS.items():
            if check == "expression_bounds_pass" and r.expression_bounds_status == "unassessed":
                continue
            # For Unsafe proteins the normal-screen holds restate normal expression, and
            # the sequence check measures overlap with Unsafe proteins, which they are.
            if unsafe and check in ("normal_screen_pass", "sequence_pass"):
                continue
            if not flag(r[check]):
                if check == "sequence_pass":
                    parts = []
                    if flag(r.readthrough_locus):
                        parts.append(f"read-through of {r.readthrough_parents}")
                    if r.unsafe_sequence_overlap_pct >= 95:
                        partners = re.findall(r"([^;()]+?) \(", str(r.sequence_partners))
                        named = list(dict.fromkeys(p.strip() for p in partners))
                        named = [p for p in named if p in unsafe_symbols]
                        parts.append(
                            f"{r.unsafe_sequence_overlap_pct:.0f}% of 9-mers in Unsafe"
                            + (f" ({', '.join(named)})" if named else "")
                        )
                    text = "sequence: " + (", ".join(parts) or "nested microprotein")
                if check == "normal_screen_pass" and isinstance(r.legacy_exclusion_reasons, str):
                    held = [
                        knowledge.loc[i]
                        for i in json.loads(r.knowledge_finding_ids or "[]")
                        if i in knowledge.index and knowledge.loc[i, "effect"] == "review"
                    ]
                    details = [
                        f"{k.tissue}, {k.assay}, interpretation {k.interpretation_status}"
                        for k in held
                    ]
                    held_by = r.legacy_exclusion_reasons.replace(" | ", ", ")
                    text += f" ({held_by}: {' / '.join(details)})" if details else f" ({held_by})"
                ways.append(text)
        fold, germ = c.rule_fold, c.germ_cell_tumor_fold
        prefix = "" if assessed else "if assessed, "
        if np.isfinite(fold) and fold < target_selection.TUMOR_FOLD:
            if np.isfinite(germ) and germ >= target_selection.TUMOR_FOLD:
                ways.append(
                    f"{prefix}only germ cell tumors reach 10x (best other {fold:.1f}x, "
                    f"germ cell {germ:.0f}x)"
                )
            else:
                ways.append(f"{prefix}no tumor at 10x (best {fold:.1f}x)")
        elif not np.isfinite(fold):
            ways.append(f"{prefix}no tumor measurement")
        if r.subtype == "Tissue":
            ways.append(f"home tissue not measured ({r.home_tissue}) and no reproductive evidence")
        if r.reproductive_evidence == "does_not_qualify":
            ways.append("no reproductive evidence")
        elif r.reproductive_evidence == "unassessed":
            ways.append("reproductive evidence not measured")
        removed = {}
        for name, d in modes.items():
            if not r.is_target_candidate:
                removed[name] = "not a candidate"
            elif bool(d.at[e, "keep"]):
                removed[name] = "kept"
            else:
                failed = [
                    LAYER_LABELS[k] for k in target_filters.LAYERS if not d.at[e, f"pass_{k}"]
                ]
                removed[name] = ", ".join(failed)
        rows.append(
            {
                "protein": r.protein_group_label,
                "ctpedia_genes": ctpedia[e],
                "classification": r.classification,
                "ways_out": len(ways),
                "reasons": "; ".join(ways),
                "unsafe_evidence": evidence.get(e, ""),
                "tumor_fold": round(float(fold), 1) if np.isfinite(fold) else np.nan,
                "fold_assessed_by_catalog": assessed,
                "best_tumor_tpm": round(float(c.rule_tumor_max_tpm), 1),
                "best_tumor_cancer": c.rule_tumor_cancer,
                "germ_cell_tumor_fold": round(float(germ), 1) if np.isfinite(germ) else np.nan,
                "normal_ceiling_tpm": round(float(c.rule_normal_tpm), 2),
                "normal_ceiling_context": c.rule_normal_context,
                "reproductive_routes": r.reproductive_routes
                if isinstance(r.reproductive_routes, str)
                else "",
                "papers_counted": int(modes["loose"].at[e, "papers_counted"]),
                **{f"{name}_filter": value for name, value in removed.items()},
            }
        )
    return pd.DataFrame(rows).sort_values(["classification", "protein"]).reset_index(drop=True)


def background_impact(catalog, matrix, background):
    entity = catalog.set_index("entity_id")
    base = entity.loc[matrix.columns, "single_tumor_threshold_tpm"]
    hits = matrix.ge(base, axis=1)
    cancer = entity.loc[matrix.columns, "is_cancer_antigen"].astype(bool)
    cancer = cancer[cancer].index
    keyed = background.set_index("sample_key").reindex(matrix.index)
    rows = []
    for cut in BACKGROUND_SWEEP:
        drop = keyed.evaluable & keyed.olfactory_share_ge5.gt(cut)
        kept = hits.loc[~drop.to_numpy(), cancer].any()
        lost = kept.index[~kept]
        rows.append(
            {
                "olfactory_share_above": f"{cut:.0%}",
                "samples_dropped": int(drop.sum()),
                "cohorts": ", ".join(
                    f"{k} {v}" for k, v in keyed[drop].cancer_code.value_counts().items()
                ),
                "cancer_antigens_lost": len(lost),
                "of_which_cta": int(entity.loc[lost, "is_cta"].astype(bool).sum()),
                "proteins": ", ".join(sorted(entity.loc[lost, "protein_group_label"])),
            }
        )
    return pd.DataFrame(rows)


def threshold_page(sweep_table, near, background, impact):
    s = sweep_table.set_index("fold")
    flagged = background[background.flagged]
    b = background[background.evaluable].olfactory_share_ge5
    hl = flagged.iloc[0] if len(flagged) else None
    groups = near.group.value_counts()
    immune = near[near.group.eq("one immune dataset")]
    names = lambda d: ", ".join(d.protein)  # noqa: E731
    broad = immune[immune.samples_ge_5x.ge(10)]
    risky = broad[broad.mapping_risk]
    broad_text = (
        f"Those with 10 or more samples at 5x ({names(broad)}) "
        + (
            "all carry paralog or host-gene mapping risk, so their tumor values are as "
            "uncertain as the normal ones."
            if len(risky) == len(broad)
            else f"include {len(broad) - len(risky)} without mapping risk."
        )
        if len(broad)
        else "None has 10 or more samples at 5x."
    )
    cols = [
        "protein",
        "group",
        "tumor_fold",
        "best_tumor_tpm",
        "best_tumor_cancer",
        "samples_ge_5x",
        "cancer_types_ge_5x",
        "normal_ceiling_tpm",
        "normal_ceiling_context",
        "datasets_at_half_ceiling",
        "fold_without_that_dataset",
        "mapping_risk",
    ]
    lines = [
        "# Tenfold rule review: thresholds, near misses and tumor-sample background",
        "",
        "Generated by `scripts/build_cta_threshold_review.py` from the 1.38.0 catalog "
        f"(policy `{target_selection.POLICY_ID}`). Nothing here changes a decision: the "
        "threshold sweep reclassifies a copy of the catalog with its own `classify()`, and "
        "the background check is a proposal.",
        "",
        "## How the CTA set moves with the fold threshold",
        "",
        "![CTA set metrics across fold thresholds](threshold-sweep.png)",
        "",
        "Only the fold changes; every other gate, the normal reference and the 10 TPM minimum "
        "(`max(N, Disposable, 1)`) stay fixed. Values: [threshold-sweep.csv](threshold-sweep.csv).",
        "",
        markdown_table(sweep_table[list(SWEEP_HEADINGS)].round(1).rename(columns=SWEEP_HEADINGS)),
        "",
        f"- From 10x down to 5x the set grows from {s.at[10, 'cta']} to {s.at[5, 'cta']} CTAs "
        f"and CTpedia proteins from {s.at[10, 'ctpedia_cta']} to {s.at[5, 'ctpedia_cta']}; "
        f"familiar CTAs: {s.at[10, 'familiar_cta']} at 10x, {s.at[5, 'familiar_cta']} at 5x, "
        f"{s.at[2, 'familiar_cta']} at 2x. The additions are lower-expressed (median best "
        f"tumor {s.at[10, 'median_best_tumor_tpm']:.0f} TPM at 10x, "
        f"{s.at[5, 'median_best_tumor_tpm']:.0f} at 5x) and more of them rest on one sample.",
        f"- From 10x up, familiar CTAs hold at {s.at[15, 'familiar_cta']} through 15x and "
        f"fall to {s.at[20, 'familiar_cta']} at 20x and {s.at[50, 'familiar_cta']} at 50x, "
        f"while the set shrinks to {s.at[20, 'cta']} at 20x.",
        f"- At 10x, {s.at[10, 'pct_one_sample']:.0f}% of CTAs rest on one tumor sample and "
        f"{s.at[10, 'pct_flagged_sample_only']:.0f}% rest only on the flagged high-background "
        "library below.",
        "- Lowering the fold adds proteins with less tumor support; raising it past 15x "
        "removes familiar CTAs (the lowest, AKAP4, is at 16x). 10x sits 1.6-fold below the "
        "first familiar loss.",
        "",
        "## Near misses: Reproductive Antigens at 8-10x",
        "",
        f"{len(near)} Reproductive Antigens have a best tumor sample at 8-10x; at 8x they "
        "would be CTAs. Each is placed in one group by what sets its normal ceiling "
        "(the highest non-reproductive normal value, or 1 TPM):",
        "",
    ]
    explain = {
        "below 10 TPM": "normal is under 1 TPM, so the bar is 10 TPM and the best tumor "
        "sample is under it. These fail on tumor level, not specificity.",
        "normal in two datasets": "at least two independent datasets put a normal value "
        "at half the ceiling or more. The ceiling is real.",
        "one brain dataset": "the ceiling comes from HPA brain regions only. Low brain "
        "expression is where a conservative ceiling matters most.",
        "one immune dataset": "the ceiling comes from one sorted-immune-cell dataset "
        "(DICE or Monaco) and the other does not show it. This is the group whose placement "
        "rests most on one source.",
        "one other dataset": "the ceiling comes from one tissue dataset.",
    }
    for group, count in groups.items():
        members = near[near.group.eq(group)]
        lines.append(f"- **{group}** ({count}): {explain[group]} {names(members)}.")
    lines += [
        "",
        f"In the one-immune-dataset group, the best tumor is a blood cancer for "
        f"{int(immune.blood_cancer_vs_immune_cells.sum())} of {len(immune)}: a leukemia or "
        "lymphoma compared with the normal immune cells it comes from. Without that dataset, "
        f"{int(immune.fold_without_that_dataset.ge(10).sum())} of {len(immune)} would reach "
        "10x. Their tumor support is narrow: a median of "
        f"{immune.samples_ge_5x.median():.0f} samples at 5x or more, in a median of "
        f"{immune.cancer_types_ge_5x.median():.0f} cancer type(s). " + broad_text,
        "",
        f"{int(near.best_sample_flagged.sum())} near misses have their best sample in the "
        "flagged high-background library below. Full table: [near-miss.csv](near-miss.csv).",
        "",
        markdown_table(near[cols].rename(columns=NEAR_HEADINGS)),
        "",
        "## Tumor-sample background (proposed check, not applied)",
        "",
        "A library contaminated with genomic DNA, or with broad unspliced background, gives "
        "low TPM to genes no tumor expresses. Olfactory receptors are a large panel of such "
        "genes. The proposed check: drop a tumor library from the cancer comparison when "
        f"more than {BACKGROUND_SHARE:.0%} of its olfactory receptor genes are at "
        f"{BACKGROUND_TPM:g} TPM or more, evaluated only when at least "
        f"{BACKGROUND_MIN_GENES} receptors are measured. Dropping a sample can only remove "
        "evidence.",
        "",
        f"- Across {len(b):,} evaluable tumor libraries the share is 0 at the median, "
        f"{b.quantile(0.99):.1%} at the 99th percentile and {b.quantile(0.999):.1%} at the "
        f"99.9th; the next-highest library after the flagged one is at "
        f"{b.sort_values().iloc[-2]:.0%}.",
    ]
    for hl in flagged.itertuples():
        peers = background[
            background.cancer_code.eq(hl.cancer_code) & background.sample_key.ne(hl.sample_key)
        ].genes_detected
        lines.append(
            f"- Flagged: `{hl.sample_key}` ({hl.cancer_code}), at "
            f"{hl.olfactory_share_ge5:.0%}, with {hl.genes_detected:,} genes detected; the "
            f"other {len(peers)} libraries in its cohort detect {peers.min():,}-"
            f"{peers.max():,}. {hl.olfactory_private_share:.0%} of its receptors at "
            f"{BACKGROUND_TPM:g} TPM are at that level in no other library of its cohort."
        )
    merkel = background[
        background.cancer_code.eq("NEC_MERKEL") & background.olfactory_share_ge5.gt(0.05)
    ]
    lines += [
        f"- {int((~background.evaluable).sum())} libraries cannot be evaluated (fewer than "
        f"{BACKGROUND_MIN_GENES} receptors in the matrix).",
        f"- Lower cutoffs reach into the Merkel cell carcinoma cohort. Its {len(merkel)} "
        "libraries above 5% share their receptors: a median of "
        f"{merkel.olfactory_private_share.median():.0%} of the receptors at "
        f"{BACKGROUND_TPM:g} TPM in one of them are absent from every other Merkel library. "
        "That reads as expression common to the tumor type, not library background, so "
        "the share alone should not drop them:",
        "",
        markdown_table(impact.rename(columns=IMPACT_HEADINGS)),
        "",
        "Per-library values: [sample-background.csv.gz](sample-background.csv.gz).",
        "",
    ]
    (OUT / "index.md").write_text("\n".join(lines))


def ctpedia_page(table: pd.DataFrame):
    counts = table.classification.value_counts()
    multi = table[table.ways_out.ge(2)]
    out = table[~table.classification.eq("CTA")]
    reason_kinds = out.reasons.str.split("; ").explode().str.replace(r" \(.*", "", regex=True)
    for prefix in ("Unsafe", "sequence", "home tissue not measured"):
        reason_kinds = reason_kinds.where(~reason_kinds.str.startswith(prefix), prefix)
    unsafe_kinds = (
        table.loc[table.classification.eq("Unsafe"), "reasons"]
        .str.extract(r"Unsafe: ([^;]*)")[0]
        .str.split(", ")
        .explode()
    )
    lines = [
        "# How CTpedia genes are classified and filtered",
        "",
        "Generated by `scripts/build_cta_threshold_review.py` from the 1.38.0 catalog "
        f"(policy `{target_selection.POLICY_ID}`, filters `{target_filters.FILTER_POLICY_ID}`).",
        "",
        f"oncoref's CTA reference tags {table.ctpedia_genes.str.split(', ').str.len().sum()} "
        f"genes as CTpedia. They are {len(table)} catalog proteins, because loci with "
        "identical selected sequences count as one protein. CTpedia tags have no "
        "primary-paper status ([sources](cta-landscape-sources.md)); they are a reference "
        "list, not evidence the catalog uses.",
        "",
        "Every protein below lists **every** check it fails, not only the one that decides "
        "its category. A protein that is out for several independent reasons says so. "
        "For Unsafe and otherwise ineligible proteins the catalog never runs the tumor "
        "comparison; this page runs the catalog's own comparison for them and marks the "
        "result *if assessed*.",
        "",
        "## Outcome",
        "",
        markdown_table(counts.rename_axis("Classification").reset_index(name="CTpedia proteins")),
        "",
        f"{len(multi)} of the {len(out)} proteins that are not CTA fail two or more checks.",
        "",
        "How often each check removes a non-CTA CTpedia protein (a protein can count under "
        "several):",
        "",
        markdown_table(
            reason_kinds.value_counts().rename_axis("Check").reset_index(name="Proteins")
        ),
        "",
        "The normal evidence behind the Unsafe calls (a protein can have several):",
        "",
        markdown_table(
            unsafe_kinds.value_counts().rename_axis("Unsafe evidence").reset_index(name="Proteins")
        ),
        "",
        "For Unsafe proteins the normal-screen holds, which restate normal expression, and "
        "the sequence check, which measures overlap with Unsafe proteins, are not listed "
        "again.",
        "",
        "## The checks",
        "",
        "- **Unsafe**: replicated normal-tissue RNA, rare normal-cell RNA, published normal "
        "expression or protein in a non-reproductive tissue "
        "([target selection](target-selection.md)).",
        "- **No tumor at 10x**: no non-germ-cell tumor sample reaches 10 x max(highest "
        "non-reproductive normal, Disposable, 1 TPM). Germ cell tumors are reported, never "
        "counted.",
        "- **No reproductive evidence**: no adult reproductive, reproductive-cell, IHC, "
        "reviewed-protein or developmental RNA route.",
        "- **Sequence**, **identity**, **coverage**, **evidence unresolved**: the shared "
        "eligibility checks.",
        "- **Filters**: `loose`, `standard` and `strict` remove Target candidates by layer; "
        "every failing layer is listed ([filters](target-filters.md)).",
        "",
        "## Filter modes",
        "",
    ]
    for mode in ("loose", "standard", "strict"):
        col = table[f"{mode}_filter"]
        kept = int(col.eq("kept").sum())
        layers = col[~col.isin(["kept", "not a candidate"])].str.split(", ").explode()
        detail = ", ".join(f"{k} {v}" for k, v in layers.value_counts().items())
        lines.append(
            f"- `{mode}` keeps {kept} CTpedia proteins; layers that remove candidates: "
            f"{detail or 'none'}."
        )
    lines += [
        "",
        "Full table: [ctpedia-genes.csv](audits/cta-threshold-review-20261116/ctpedia-genes.csv).",
        "",
    ]
    order = [
        "CTA",
        "Cancer Antigen only",
        "Reproductive Antigen only",
        "Other:Expression",
        "Other:Tissue",
        "Other:Review",
        "Other:Sequence",
        "Unsafe",
    ]
    headings = {
        "protein": "Protein",
        "ctpedia_genes": "CTpedia genes",
        "ways_out": "Ways out",
        "reasons": "Reasons (all that apply)",
        "unsafe_evidence": "Highest Unsafe evidence per layer",
        "tumor_fold": "Best tumor fold",
        "best_tumor_tpm": "Best tumor TPM",
        "best_tumor_cancer": "Cancer",
        "normal_ceiling_tpm": "Normal ceiling TPM",
        "normal_ceiling_context": "Normal ceiling context",
        "germ_cell_tumor_fold": "Germ cell tumor fold",
        "loose_filter": "loose",
        "standard_filter": "standard",
        "strict_filter": "strict",
    }
    tumor = ["tumor_fold", "best_tumor_tpm", "best_tumor_cancer", "normal_ceiling_tpm"]
    filters = ["loose_filter", "standard_filter", "strict_filter"]
    shown = {
        "CTA": [
            "protein",
            "ctpedia_genes",
            *tumor,
            "normal_ceiling_context",
            "germ_cell_tumor_fold",
            *filters,
        ],
        "Unsafe": ["protein", "ctpedia_genes", "ways_out", "reasons", "unsafe_evidence", *tumor],
    }
    other = [
        "protein",
        "ctpedia_genes",
        "ways_out",
        "reasons",
        *tumor,
        "normal_ceiling_context",
        "germ_cell_tumor_fold",
        *filters,
    ]
    notes = {
        "CTA": "Kept by the catalog. The filter columns list every layer that removes the "
        "protein in that mode.",
        "Unsafe": "Excluded by normal expression. Tumor columns are *if assessed*: the "
        "catalog's own comparison, which it does not run for Unsafe proteins.",
    }
    for name in order + sorted(set(table.classification) - set(order)):
        part = table[table.classification.eq(name)]
        if part.empty:
            continue
        lines += [f"## {name} ({len(part)})", ""]
        if name in notes:
            lines += [notes[name], ""]
        columns = shown.get(name, other)
        lines += [markdown_table(part[columns].rename(columns=headings)), ""]
    CTPEDIA_PAGE.write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    catalog, matrix, cancer_code, normal = load()
    families = pd.read_csv(ROOT / "oncoref/data/cancer-type-registry.csv").set_index("code").family
    background = sample_background()
    write_csv(background, OUT / "sample-background.csv.gz")
    flagged = set(background.loc[background.flagged, "sample_key"])
    ctpedia = ctpedia_entities(catalog)
    table = sweep(catalog, matrix, cancer_code, flagged, ctpedia)
    if not reclassify(catalog, 10.0).classification.equals(catalog.classification):
        raise ValueError("Reclassifying at 10x does not reproduce the catalog")
    write_csv(table.round(2), OUT / "threshold-sweep.csv")
    sweep_figure(table, OUT / "threshold-sweep.png")
    near = near_misses(catalog, matrix, cancer_code, normal, background, families)
    write_csv(near, OUT / "near-miss.csv")
    impact = background_impact(catalog, matrix, background)
    write_csv(impact, OUT / "background-cutoffs.csv")
    comparison = ctpedia_comparison(set(ctpedia.index))
    genes = ctpedia_table(catalog, ctpedia, comparison)
    write_csv(genes, OUT / "ctpedia-genes.csv")
    threshold_page(table, near, background, impact)
    ctpedia_page(genes)
    print(table.to_string(), flush=True)
    print(genes.classification.value_counts().to_string(), flush=True)


if __name__ == "__main__":
    main()
