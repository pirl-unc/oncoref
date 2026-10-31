#!/usr/bin/env python3
"""Audit released Neither boundaries without changing classification policy.

Reads committed 1.24.0 evidence tables. Writes review queues and a numerical
report, not replacement classifications. TPM, nTPM and CP10K remain separate;
the current policy's CP10K floor is used only for a labelled diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from oncoref.cta_selection_policy import CANDIDATE_ONLY_LAYERS

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/audits/cta-selection-20261028"
PROFILES = ROOT / "docs/audits/cta-unassigned-20261029"
DEFAULT_OUTPUT = ROOT / "docs/audits/cta-neither-20261030"


def markdown(frame: pd.DataFrame) -> str:
    """Small dependency-free Markdown table with explicit missing values."""

    def cell(value: object) -> str:
        if pd.isna(value):
            return "not measured / unavailable"
        if isinstance(value, float):
            return f"{value:.4g}"
        return str(value).replace("|", ";").replace("\n", " ")

    rows = [list(frame.columns), ["---"] * len(frame.columns)]
    rows.extend(frame.itertuples(index=False, name=None))
    return "\n".join("| " + " | ".join(map(cell, row)) + " |" for row in rows)


def build(output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    inputs = [
        DEFAULT_OUTPUT / "baseline-1.24.0.csv.gz",
        SOURCE / "normal-peak-evidence.csv.gz",
        SOURCE / "hpa-protein-evidence.csv.gz",
        PROFILES / "hpa-v23-protein-profiles.csv.gz",
    ]
    f, peaks, ihc, hpa = [pd.read_csv(path, low_memory=False) for path in inputs]
    if not f.protein_equivalence_group_id.is_unique:
        raise ValueError("Protein identities must be unique")
    policy = f.policy.unique().tolist()
    if policy != ["oncoref.cta-selection.2026-10-29.4"]:
        raise ValueError(f"This audit is for the released 1.24.0 evidence: {policy}")
    safety_columns = [c for c in f if c.startswith("hold:") and c[5:] not in CANDIDATE_ONLY_LAYERS]
    f["normal_screen_hold"] = f[safety_columns].any(axis=1)
    f["normal_panel_pass"] = f.required_coverage & ~f.unresolved_evidence & ~f.normal_screen_hold
    f["recorded_category"] = f.category.fillna("Unassigned")
    unassigned = f.category.isna()
    clear = unassigned & f.normal_panel_pass
    neither = f.category.eq("Neither")
    independent_vetoes = [
        "unsafe:HPA normal protein",
        "unsafe:Accepted normal-protein evidence",
        "unsafe:Clinical safety exclusion",
        "unsafe:External validated normal protein",
    ]
    f["independent_protein_or_clinical_veto"] = f[independent_vetoes].any(axis=1)
    worst = peaks.loc[peaks.unsafe_evidence].groupby("gene_id").ratio.max()
    f["largest_unsafe_rna_ceiling_ratio"] = f.representative_gene_id.map(worst)
    boundary = (
        f.category.eq("Unsafe")
        & ~f.independent_protein_or_clinical_veto
        & f.largest_unsafe_rna_ceiling_ratio.between(1, 1.25, inclusive="right")
    )
    f["nonpermitted_normal_above_30_TPM_or_nTPM"] = f.healthy_maximum_ntpm.gt(
        30
    ) | f.healthy_maximum_tpm.gt(30)
    # Positive reproductive IHC is evidence for review, not automatic admission.
    reproductive = ihc.category.eq("reproductive")
    positive = ihc.Level.isin(["Low", "Medium", "High"])
    reliable = ihc.Reliability.isin(["Approved", "Supported", "Enhanced"])
    strong = ihc.Level.isin(["Medium", "High"]) & ihc.Reliability.isin(["Supported", "Enhanced"])
    measured = set(ihc.loc[reproductive & positive & reliable, "gene_id"])
    strongly_measured = set(ihc.loc[reproductive & strong, "gene_id"])
    members = f.member_gene_ids.map(json.loads)
    f["positive_reproductive_ihc_all_members"] = members.map(lambda ids: set(ids) <= measured)
    f["strong_reproductive_ihc_all_members"] = members.map(
        lambda ids: set(ids) <= strongly_measured
    )

    selected = [
        "protein_equivalence_group_id",
        "protein_group_label",
        "representative_gene_id",
        "member_gene_ids",
        "recorded_category",
        "assessment_status",
        "neither_profile",
        "required_coverage",
        "unresolved_evidence",
        "normal_screen_hold",
        "reproductive_ntpm",
        "developmental_mean_tpm",
        "developmental_all_lineage_max_tpm",
        "healthy_maximum_ntpm",
        "healthy_maximum_tpm",
        "healthy_maximum_cp10k",
        "oncoref_maximum_tpm",
        "oncoref_maximum_sample",
        "oncoref_cancer_code",
        "hpa_tumor_group_lower_bound_ptpm",
        "tumor_union_ge2",
        "tumor_profile_complete",
        "optional_DICE_complete",
        "optional_cornea_complete",
        "any_detected_protein",
        "manageable_tissue_count",
        "hpa_diffuseness",
        "positive_reproductive_ihc_all_members",
        "strong_reproductive_ihc_all_members",
        "independent_protein_or_clinical_veto",
        "largest_unsafe_rna_ceiling_ratio",
        "nonpermitted_normal_above_30_TPM_or_nTPM",
        "exclusion_reasons",
        "Unsafe_reasons",
    ]

    def export(name: str, frame: pd.DataFrame) -> None:
        frame.to_csv(output / name, index=False, float_format="%.9g")

    export("neither-individual-data.csv", f.loc[neither, selected])
    export("unassigned-normal-pass.csv", f.loc[clear, selected])
    export("unassigned-all.csv", f.loc[unassigned, selected])
    export(
        "reproductive-ihc-review.csv",
        f.loc[clear & f.positive_reproductive_ihc_all_members, selected],
    )
    export("unsafe-rna-boundaries.csv", f.loc[boundary, selected])
    family = f.protein_group_label.str.match(r"^(USP17L|FAM90A|TAF11L|KRT|LCE2|CRYG)")
    export("family-neighbors.csv", f.loc[family, selected])
    peak_ids = set(f.loc[neither | boundary | family, "representative_gene_id"])
    export("review-normal-peaks.csv.gz", peaks.loc[peaks.gene_id.isin(peak_ids)])
    bulk = f.loc[clear, selected].merge(
        hpa[["entity_id", "tonsil", "appendix", "gallbladder"]],
        left_on="protein_equivalence_group_id",
        right_on="entity_id",
        validate="one_to_one",
    )
    bulk = bulk[bulk[["tonsil", "appendix", "gallbladder"]].max(axis=1).ge(2)]
    export("disposable-bulk-examples.csv", bulk)
    examples = f[f.protein_group_label.isin(["USP17L3", "FAM90A9"])]
    export("tumor-positive-low-normal-examples.csv", examples[selected])

    counts = f.recorded_category.value_counts().to_dict()
    clear_count = int(clear.sum())
    receipt = {
        "policy": policy[0],
        "scope": "diagnostic only; no released category or threshold changed",
        "inputs": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "current_counts": counts,
        "unassigned_normal_panel_pass": clear_count,
        "unassigned_normal_pass_tumor_positive": int((clear & f.tumor_union_ge2).sum()),
        "unassigned_normal_pass_reproductive_ihc": int(
            (clear & f.positive_reproductive_ihc_all_members).sum()
        ),
        "unassigned_normal_pass_strong_reproductive_ihc": int(
            (clear & f.strong_reproductive_ihc_all_members).sum()
        ),
        "unassigned_normal_pass_strong_reproductive_ihc_tumor_positive": int(
            (clear & f.strong_reproductive_ihc_all_members & f.tumor_union_ge2).sum()
        ),
        "neither_with_normal_hold": int((neither & f.normal_screen_hold).sum()),
        "neither_nonpermitted_normal_above_30_TPM_or_nTPM": int(
            (neither & f.nonpermitted_normal_above_30_TPM_or_nTPM).sum()
        ),
        "unsafe_rna_only_within_25_percent_of_ceiling": int(boundary.sum()),
        "fallback_counterfactual": {
            **counts,
            "Neither": counts["Neither"] + clear_count,
            "Unassigned": counts["Unassigned"] - clear_count,
        },
    }
    if sum(counts.values()) != len(f) or sum(receipt["fallback_counterfactual"].values()) != len(f):
        raise ValueError("Counterfactual must preserve the protein partition")
    if (clear & f[[c for c in f if c.startswith("unsafe:")]].any(axis=1)).any():
        raise ValueError("The fallback queue must contain no recorded Unsafe evidence")
    (output / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    count_table = pd.DataFrame(
        {"Current": counts, "Fallback only": receipt["fallback_counterfactual"]}
    )
    count_table.index.name = "Category"
    report = [
        "# Numerical audit of Neither and Unassigned\n",
        f"Evidence: OncoRef 1.24.0, `{policy[0]}`. Reproduce with "
        "`.venv/bin/python scripts/audit_cta_neither_boundaries.py`. "
        "This report does not change the released assignments.\n",
        "## Structural fallback counterfactual\n",
        markdown(count_table.reset_index()),
        f"\n{clear_count} Unassigned proteins pass required normal coverage, have no unresolved "
        "evidence and have no normal-screen hold. 757 fail only reproductive association; "
        "six additionally lack a qualifying cancer sample. This is a mechanical counterfactual, "
        "not a claim that these proteins are clinically safe or all belong in Neither.\n",
        "## Important audit counts\n",
        markdown(
            pd.DataFrame(
                [
                    (
                        "Current Neither with a normal-screen hold",
                        receipt["neither_with_normal_hold"],
                    ),
                    (
                        "Current Neither with non-permitted normal >30 TPM or nTPM",
                        receipt["neither_nonpermitted_normal_above_30_TPM_or_nTPM"],
                    ),
                    (
                        "Normal-pass Unassigned with positive reproductive IHC in every member",
                        receipt["unassigned_normal_pass_reproductive_ihc"],
                    ),
                    (
                        "Subset with Medium/High, Supported/Enhanced reproductive IHC",
                        receipt["unassigned_normal_pass_strong_reproductive_ihc"],
                    ),
                    (
                        "Unsafe with RNA-only exclusions <=25% over every Unsafe ceiling",
                        receipt["unsafe_rna_only_within_25_percent_of_ceiling"],
                    ),
                ],
                columns=["Finding", "Proteins"],
            )
        ),
        "\nThe IHC rows are review evidence, not a new admission rule. Boundary proteins may "
        "still have Manageable excess or other admission failures; removing one RNA exclusion "
        "would not necessarily make them candidates.\n",
        "## All current Neither proteins\n",
        markdown(
            f.loc[
                neither,
                [
                    "protein_group_label",
                    "neither_profile",
                    "reproductive_ntpm",
                    "healthy_maximum_ntpm",
                    "healthy_maximum_tpm",
                    "oncoref_maximum_tpm",
                    "normal_screen_hold",
                ],
            ]
        ),
        "\nNormal maxima are supported non-permitted normal compartments, not maxima over "
        "every donor or every tissue in the body. Source units remain separate. The tumor "
        "column is the OncoRef tumor maximum; HPA pTPM is retained separately in the CSV.\n",
        "## Tumor-positive, low-normal examples\n",
        markdown(
            examples[
                [
                    "protein_group_label",
                    "healthy_maximum_ntpm",
                    "healthy_maximum_tpm",
                    "reproductive_ntpm",
                    "developmental_all_lineage_max_tpm",
                    "oncoref_maximum_tpm",
                    "oncoref_cancer_code",
                    "hpa_tumor_group_lower_bound_ptpm",
                    "optional_DICE_complete",
                    "optional_cornea_complete",
                ]
            ]
        ),
        "\nNeither example has a measured HPA tissue-IHC row in the pinned group-IHC table. "
        "No positive IHC is therefore not equivalent to negative IHC. Tumor maxima do not "
        "establish prevalence, malignant-cell origin, or protein production.\n",
        "## Disposable bulk examples\n",
        markdown(
            bulk[
                [
                    "protein_group_label",
                    "tonsil",
                    "appendix",
                    "gallbladder",
                    "healthy_maximum_ntpm",
                    "healthy_maximum_tpm",
                ]
            ]
        ),
        "\nThese are examples with permitted bulk RNA >=2 nTPM, not an exhaustive "
        "Disposable-only set: hair-follicle and melanocyte compartments require their "
        "specific cell data. Skin as a whole is Manageable.\n",
        "## Interpretation and limitations\n",
        "- `low_expression` unnecessarily requires low tumor expression and measured developmental "
        "profiles; those requirements do not establish adult-normal safety.\n"
        "- `diffuse_manageable` counts >=2 Manageable tissue labels with RNA >=2. It can allow "
        "above-ceiling expression and counts thymus, unlike the separate organic diffuseness score.\n"
        "- Current reproductive admission is RNA-only (>=5); positive reproductive IHC is not "
        "an independent admission route.\n"
        "- Required-panel completeness is not exhaustive biological coverage. Rare cells, activation "
        "states and paralog-specific read assignment remain material limitations.\n"
        "- Low observed expression is not proof of absence. Optional-source missingness remains explicit.\n"
        "- Review queues are overlapping diagnostics, not disjoint replacement categories.\n",
    ]
    (output / "numbers.md").write_text("\n".join(report))
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build(args.output)
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in {"inputs", "script_sha256"}}, indent=2
        )
    )
