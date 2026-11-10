#!/usr/bin/env python3
"""Counterfactual target maxima and Disposable caps; never rewrite the catalog."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import target_selection, tissue_hazards  # noqa: E402
from oncoref.gene_matching import (  # noqa: E402
    gene_identity_catalog,
    resolve_long_measurement_ids,
)

REBUILD = ROOT / "docs/audits/target-selection-20261109/rebuild"
OUT = REBUILD.parent / "expression-maxima"
CAPS = (1, 4, 8, 16, 32, 64, 128)


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def evaluate(target, normal, disposable, *, cap=None):
    """Strict maximum, inclusive cap, missing observations remain unassessed.

    Inputs must already be on one declared measurement scale. This function
    cannot establish source comparability or completeness of tissue coverage.
    """
    for value in (target, normal, disposable):
        if pd.notna(value) and (not math.isfinite(value) or value < 0):
            raise ValueError("Expression must be nonnegative and finite")
    if cap is not None and (not math.isfinite(cap) or cap < 0):
        raise ValueError("Cap must be nonnegative and finite")
    measured = all(pd.notna(v) for v in (target, normal, disposable))
    maximum = (
        ("passes" if target >= 1 and target > max(normal, disposable) else "fails")
        if measured
        else "unassessed"
    )
    capped = (
        "not_applied"
        if cap is None
        else ("unassessed" if pd.isna(disposable) else "passes" if disposable <= cap else "fails")
    )
    return {"maximum_status": maximum, "cap_status": capped}


def maxima(frame, scope):
    d = frame[frame.scope.eq(scope)]
    values = d.groupby("entity_id").value.max()
    witnesses = d[d.value.eq(d.entity_id.map(values))].copy()
    witnesses["witness"] = witnesses.apply(
        lambda r: f"{r.source_id}: {r.source_context} ({r.statistic})", axis=1
    )
    return pd.DataFrame(
        {
            f"{scope}_max_tpm": values,
            f"{scope}_witnesses": witnesses.groupby("entity_id").witness.agg(
                lambda v: json.dumps(sorted(set(v)))
            ),
        }
    )


def audit_cells(path, candidates):
    """Retain source observations, identities and read support on native nTPM."""
    contexts, selected = {}, []
    for chunk in pd.read_csv(path, sep="\t", chunksize=300_000):
        keys = list(zip(chunk.Tissue, chunk["Cell type"]))
        for tissue, cell in set(keys) - contexts.keys():
            contexts[(tissue, cell)] = tissue_hazards.assess(tissue, cell, resolution="cell")
        keep = [contexts[key]["tier"] == "disposable" for key in keys]
        selected.append(chunk.loc[keep])
    inventory = pd.DataFrame(
        {"tissue": t, "cell_type": c, "tier": result["tier"], "reason": result["reason"]}
        for (t, c), result in sorted(contexts.items())
    )
    inventory.to_csv(OUT / "hpa-cell-contexts.csv", index=False)
    raw = pd.concat(selected, ignore_index=True).fillna({"Gene name": ""})
    crosswalk, decisions = resolve_long_measurement_ids(
        raw,
        id_column="Gene",
        symbol_column="Gene name",
        context_columns=["Tissue", "Cluster", "Cell type"],
        value_columns=["Read count", "nTPM"],
        source="hpa_single_cell_tissue",
        source_version="v23",
    )
    decisions.to_csv(
        OUT / "hpa-disposable-identity-decisions.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    raw = raw.merge(crosswalk, left_on="Gene", right_on="source_gene_id", validate="many_to_one")
    members = candidates[["entity_id", "protein", "category", "member_gene_ids"]].copy()
    members["gene_id"] = members.member_gene_ids.map(json.loads)
    members = members.explode("gene_id").drop(columns="member_gene_ids")
    observed = raw.merge(
        members, left_on="canonical_gene_id", right_on="gene_id", validate="many_to_one"
    )
    observed["unit"] = "nTPM"
    observed["read_support_pass"] = observed["Read count"].ge(20)
    observed.to_csv(OUT / "hpa-disposable-cell-observations.csv", index=False)
    missing = members[~members.gene_id.isin(observed.gene_id)]
    missing.to_csv(OUT / "hpa-disposable-cell-missing-members.csv", index=False)
    supported = observed[observed.read_support_pass]
    return {
        "source": "hpa_single_cell_tissue:v23",
        "source_sha256": digest(path),
        "identity_catalog_sha256": gene_identity_catalog().identity_catalog_sha256.iloc[0],
        "unit": "nTPM",
        "disposable_contexts": inventory[inventory.tier.eq("disposable")].to_dict("records"),
        "proteins_with_any_mapped_observation": observed.entity_id.nunique(),
        "missing_member_loci": len(missing),
        "proteins_with_supported_positive_observation": supported[
            supported.nTPM.gt(0)
        ].entity_id.nunique(),
        "unsupported_positive_member_observations": len(
            observed[observed.nTPM.gt(0) & ~observed.read_support_pass]
        ),
        "strongest_observed_signal": observed.sort_values("nTPM", ascending=False)
        .head(1)[["protein", "Tissue", "Cell type", "Read count", "nTPM", "read_support_pass"]]
        .to_dict("records"),
        "supported_cap_exceedances": {
            str(cap): supported[supported.nTPM.gt(cap)].entity_id.nunique() for cap in CAPS
        },
    }


def build(cell_path):
    OUT.mkdir(parents=True, exist_ok=True)
    catalog = target_selection.results()
    x = (
        catalog[catalog.is_target_candidate][
            [
                "entity_id",
                "protein_group_label",
                "member_gene_ids",
                "category",
                "tumor_max_tpm",
                "tumor_top_sample",
                "tumor_top_cancer_code",
                "developmental_mean_tpm",
            ]
        ]
        .rename(columns={"protein_group_label": "protein"})
        .set_index("entity_id")
    )
    raw = pd.read_csv(REBUILD / "normal-comparison.csv.gz")
    if set(raw.unit) != {"TPM"}:
        raise ValueError("Do not pool expression scales")
    panel = raw[raw.qc_pass & raw.identity_resolved & raw.value.notna()]
    for scope in ("normal", "reproductive", "disposable"):
        x = x.join(maxima(panel, scope))
    positive = ["tumor_max_tpm", "reproductive_max_tpm", "developmental_mean_tpm"]
    x["target_max_tpm"] = x[positive].max(axis=1)
    x["target_max_source"] = x[positive].idxmax(axis=1)
    x["other_max_tpm"] = x[["normal_max_tpm", "disposable_max_tpm"]].max(axis=1)
    status = pd.DataFrame(
        [
            evaluate(r.target_max_tpm, r.normal_max_tpm, r.disposable_max_tpm)
            for r in x.itertuples()
        ],
        index=x.index,
    )
    x["maximum_status"] = status.maximum_status
    x = x.reset_index()
    x.to_csv(OUT / "protein-maxima.csv", index=False)
    scenarios = []
    losses = []
    for cap in (None, *CAPS):
        statuses = pd.DataFrame(
            [
                evaluate(r.target_max_tpm, r.normal_max_tpm, r.disposable_max_tpm, cap=cap)
                for r in x.itertuples()
            ]
        )
        lost = statuses.eq("fails").any(axis=1)
        unassessed = ~lost & statuses.eq("unassessed").any(axis=1)
        passed = ~lost & ~unassessed
        counts = x.loc[passed, "category"].value_counts()
        scenarios.append(
            {
                "disposable_cap_tpm": cap,
                "moved_to_other": int(lost.sum()),
                "unassessed_in_measured_panel": int(unassessed.sum()),
                **{
                    c: int(counts.get(c, 0))
                    for c in ["CTA", "Cancer Antigen only", "Reproductive Antigen only"]
                },
            }
        )
        losses.append(x[lost].assign(disposable_cap_tpm=cap))
    scenarios = pd.DataFrame(scenarios)
    scenarios.to_csv(OUT / "cap-scenarios.csv", index=False)
    pd.concat(losses, ignore_index=True).to_csv(OUT / "proteins-moved-to-other.csv", index=False)
    inventory = tissue_hazards.anatomy_inventory()
    inventory.to_csv(OUT / "anatomy-inventory.csv", index=False)
    cells = audit_cells(cell_path, x)
    summary = {
        "policy_id": target_selection.POLICY_ID,
        "catalog_sha256": digest(ROOT / "oncoref/data/target-selection.csv.gz"),
        "normal_panel_sha256": digest(REBUILD / "normal-comparison.csv.gz"),
        "candidate_proteins": len(x),
        "maximum_status": x.maximum_status.value_counts().to_dict(),
        "disposable_tpm_contexts": sorted(
            panel.loc[panel.scope.eq("disposable"), "tissue"].unique()
        ),
        "anatomy_selection_categories": inventory.selection_category.fillna("Needs localization")
        .value_counts()
        .to_dict(),
        "specific_registered_contexts_needing_hazard_review": inventory[
            inventory.requires_hazard_review
        ].concept_id.tolist(),
        "broad_location_labels": inventory[inventory.requires_localization].concept_id.tolist(),
        "hpa_cell_audit": cells,
        "cap_scenarios": scenarios.astype(object).where(scenarios.notna(), None).to_dict("records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    cap_table = "\n".join(
        "| "
        + " | ".join(
            str(v)
            for v in [
                "None" if pd.isna(row.disposable_cap_tpm) else int(row.disposable_cap_tpm),
                int(row.moved_to_other),
                int(row.CTA),
                int(row["Cancer Antigen only"]),
                int(row["Reproductive Antigen only"]),
            ]
        )
        + " |"
        for _, row in scenarios.iterrows()
    )
    (OUT / "index.md").write_text(f"""# Expression maxima and Disposable caps

**Requiring cancer or reproductive/developmental expression to exceed every measured nonreproductive value removes zero of the current {len(x)} Target candidates.** This is an additional filter on the current catalog. The original cancer gate, reproductive association and shared eligibility rules still apply.

For each protein, compare the largest tumor-sample TPM, adult reproductive-tissue TPM or supported developmental-lineage mean TPM against the largest nonreproductive normal or Disposable value. The target value must be at least 1 TPM and strictly greater; ties fail. One tumor sample can establish the maximum. Protein equivalence records retain the rebuild's within-context member aggregation; maxima from different samples are never added.

## Absolute Disposable cap

The cap is inclusive: a value equal to it passes. Exceeding it moves a candidate to Other, without creating new Unsafe evidence. These are sensitivity analyses; no cap has been adopted.

| Disposable cap (native TPM) | Move to Other | CTA remaining | Cancer Antigen only remaining | Reproductive Antigen only remaining |
|---:|---:|---:|---:|---:|
{cap_table}

LCE3C has the highest Disposable value, **21.4 TPM in tonsil**. NLRP4 is next at **4.6 TPM in tonsil**. The [complete protein table](protein-maxima.csv) retains each maximum, source, tissue and tumor sample. [Every changed protein at every cap](proteins-moved-to-other.csv) is listed separately.

## What this comparison measures

The TPM panel's Disposable tissues are **appendix, gallbladder and tonsil**. It contains no isolated hair-follicle or cutaneous-melanocyte TPM measurement. Thus zero losses means zero in this measured panel, not a complete all-tissue result. Normal values include tissue means and donor medians; this is not a maximum over every healthy individual. The underlying [tumor-source comparability concerns](../scrutiny/index.md) remain.

[HPA v23](https://v23.proteinatlas.org/about/download) has a **skin-melanocyte cluster** on a separate nTPM scale. {cells["proteins_with_any_mapped_observation"]} candidates have a mapped member observation there; {cells["missing_member_loci"]} member loci are missing. No candidate has a positive signal meeting the existing 20-read requirement. The largest observed signal is **TSGA10IP: 24 nTPM, 14 reads**; it is retained as insufficiently supported, not converted to zero. HPA's source labels contain no hair-follicle cluster. [Cell observations and read support](hpa-disposable-cell-observations.csv), [missing members](hpa-disposable-cell-missing-members.csv) and [all source cell labels](hpa-cell-contexts.csv) are retained. TPM and nTPM are not pooled or treated as interchangeable. The cell result adds no supported cap exceedance on its own scale and does not fill the missing hair-follicle panel.

The proposed rule has two independent parts: **target expression must be highest**, and **Disposable expression must remain below a fixed ceiling**. A large tumor value cannot excuse unlimited Disposable expression. Required but missing normal measurements belong in Other:Review; an observed cap failure belongs in Other:Expression. The current required-coverage panel has not been expanded by this audit.

## Reproductive contexts and location labels

**Reproductive includes adult reproductive organs, placenta/trophoblast, extraembryonic membranes and structures, embryo and fetus.** The taxonomy already contains this union. The inventory now records selection scope separately from injury tier and uses the correct embryonic/fetal stage. Fourteen of the previously unresolved hazard entries belong to this Reproductive scope. Reproductive membership does not exempt the same immune, vascular or other shared cells at nonreproductive sites.

**Soft tissue, head and foot are location labels requiring localization.** They receive no new hazard tier and no Disposable allowance. Retain their observations and resolve the actual tissue/cell identity before using them for clearance. Among the currently registered specific human concepts, none remains without either a Reproductive selection scope or an adult hazard tier. This is a registry inventory, not a complete list of human tissues or source coverage.

Three concrete additions need separate curation:

| Compartment | Why the existing label is insufficient | Source to support mapping and expression review |
|---|---|---|
| Olfactory sensory neurons, sustentacular cells and basal progenitors | Nasopharynx and generic respiratory epithelium do not identify the olfactory neuroepithelium. | [Durante et al., 2020](https://www.nature.com/articles/s41593-020-0587-9), [GSE139522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522) |
| Mineralized bone; osteoblasts, osteocytes and osteoclasts | Bone marrow and cartilage do not cover these cells. | [Human osteoblast PANX3 study](https://pubmed.ncbi.nlm.nih.gov/32982413/) provides a candidate-specific lead; its cell-line evidence does not establish an adult tissue hazard tier. |
| Dental pulp, odontoblasts and periodontal cells | Oral mucosa does not identify the tooth's living compartments. Adult cells belong in normal-tissue assessment; fetal specimens also retain developmental scope. | [Pagella et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33997688/) |

These studies support distinct normal compartments; they do not establish tolerance of immune-mediated depletion. New hazard calls require the registry's functional-injury rationale and evidence.

[Current anatomy inventory](anatomy-inventory.csv) · [Scenario counts](cap-scenarios.csv) · [Summary](summary.json) · [Provenance manifest](manifest.json).

Reproduce: `python scripts/audit_target_expression_maxima.py --hpa-cell-file PATH_TO_V23_RNA_SINGLE_CELL_TYPE_TISSUE_TSV`.
""")
    manifest = {
        "implementation_sha256": {
            p: digest(ROOT / p)
            for p in [
                "scripts/audit_target_expression_maxima.py",
                "oncoref/tissue_hazards.py",
                "oncoref/gene_matching.py",
            ]
        },
        "inputs_sha256": {
            str(p.relative_to(ROOT)): digest(p)
            for p in [
                ROOT / "oncoref/data/target-selection.csv.gz",
                REBUILD / "normal-comparison.csv.gz",
                ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
                ROOT / "oncoref/data/normal-tissue-hazards.yaml",
            ]
        },
        "external_input": {
            "source": "hpa_single_cell_tissue:v23",
            "path": str(cell_path),
            "sha256": digest(cell_path),
        },
        "outputs_sha256": {
            p.name: digest(p)
            for p in sorted(OUT.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hpa-cell-file", type=Path, required=True)
    build(parser.parse_args().hpa_cell_file)
