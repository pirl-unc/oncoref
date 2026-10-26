#!/usr/bin/env python3
"""Replay pinned measurement gates, apply reviewed findings, rebuild current sets.

No new download or literature search occurs during a build. Use a fresh output
directory. --publish-dataset updates the packaged API only after verification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from target_sets_report import proteins  # noqa: E402

from oncoref import cta_policy as policy  # noqa: E402

BASE = ROOT / "docs/audits/cta-neither-review-20261024"
EXPRESSION = ROOT / "docs/audits/normal-expression-categories-20261021-r2"
LITERATURE = ROOT / "docs/audits/cta-literature-20261025"
OUT = ROOT / "docs/audits/cta-policy-20261025"
ADJUDICATION = ROOT / "docs/audits/cta-adjudication-20261025"
RELEASE = "https://github.com/pirl-unc/oncoref/releases/tag/v1.20.0"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return pd.read_csv(path, keep_default_na=False, low_memory=False)


def write(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


def dump(value, path):
    path.write_text(json.dumps(value, indent=2) + "\n")


def verify_report(folder):
    m = json.loads((folder / "run-manifest.json").read_text())
    for name, expected in m["outputs"].items():
        if digest(folder / name) != expected:
            raise ValueError(f"Pinned report changed: {folder / name}")
    return m


def decision_facts(base, categories, reference):
    """Adapt frozen measurement outputs; never use their final Category as a rule."""
    c = categories.set_index("Ensembl_Gene_ID").reindex(base["Gene ID"]).reset_index(drop=True)
    if c.category.isna().any():
        raise ValueError("Missing reference-bound expression row")
    f = base[["Gene", "Gene ID", "Reproductive RNA (nTPM)"]].copy()
    f["Reproductive RNA (nTPM)"] = pd.to_numeric(f["Reproductive RNA (nTPM)"], errors="coerce")
    f["Taxonomy ID"], f["Reference"] = 9606, reference
    f["Nominated"] = base["Previous CTA screen"].ne("Not nominated")
    f["Restriction passed"] = base["Previous CTA screen"].eq("Retained")
    f["Coverage adequate"] = (
        c.rna_routine_missing.eq(0)
        & c.quality_flags.eq("")
        & f["Reproductive RNA (nTPM)"].notna()
        & pd.to_numeric(base["Maximum RNA (nTPM)"], errors="coerce").notna()
        & base["Normal RNA assessed"]
    )
    safety = read(ROOT / "oncoref/data/cta-external-safety-reviews.csv")
    f["Safety veto"] = f["Gene ID"].isin(safety.loc[safety.action.eq("hold"), "Ensembl_Gene_ID"])
    f["Normal IHC positive"] = base["Normal protein positive"]
    f["Normal RNA failure"] = (
        base["Normal RNA above ceiling"] | base["Normal donor RNA above ceiling"]
    )
    f["High normal RNA"] = base["High normal RNA"]
    f["Variable"] = base["Previous CTA screen"].eq("Variable")
    f["Antibody concern"] = base["Antibody review"]
    f["Cell RNA concern"] = base[
        ["Cell RNA review", "Immune RNA review", "Corneal RNA review"]
    ].any(axis=1)
    f["Reproductive protein"] = (
        base["Reproductive positive control"] | base["Additional reproductive protein"]
    )
    return f


def markdown(frame, columns):
    def esc(x):
        return str(x).replace("|", "/").replace("\n", " ")

    return "\n".join(
        [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join("---" for _ in columns) + " |",
            *[
                "| " + " | ".join(esc(x) for x in row) + " |"
                for row in frame[columns].itertuples(index=False, name=None)
            ],
        ]
    )


def build(out, assets, previous_assets, publish=False):
    for folder in (out, assets):
        if folder.exists() and any(folder.iterdir()):
            raise ValueError(f"Use an empty output directory: {folder}")
        folder.mkdir(parents=True, exist_ok=True)
    verify_report(BASE)
    reference = json.loads(
        (ROOT / "docs/audits/normal-sequence-sets-20261021/run-manifest.json").read_text()
    )["reference"]["reference_id"]
    base = read(BASE / "gene-results.csv.gz")
    categories = read(EXPRESSION / "gene-categories.csv.gz")
    em = json.loads((EXPRESSION / "run-manifest.json").read_text())
    if digest(EXPRESSION / "gene-categories.csv.gz") != em["outputs"]["gene-categories.csv.gz"]:
        raise ValueError("Expression-category checksum mismatch")
    facts = decision_facts(base, categories, reference)
    reviews = read(ROOT / "oncoref/data/cta-policy-findings.csv")
    replay = policy.classify(facts, reviews.loc[reviews["Review ID"].str.startswith("legacy-")])
    if not replay.Category.equals(base.Category):
        raise AssertionError("The explicit measurement facts do not reproduce the frozen baseline")
    decisions = policy.classify(facts, reviews)
    adjudications = read(ROOT / "oncoref/data/cta-antibody-adjudications.csv")
    draft = read(ADJUDICATION / "baseline-gene-results.csv.gz")
    raw_ids = set(read(LITERATURE / "retained-primary-normal-signals.csv").gene_id)
    required = raw_ids & set(draft.loc[draft.Category.eq("CTA"), "Gene ID"])
    if set(adjudications["Gene ID"]) != required or adjudications["Gene ID"].duplicated().any():
        raise AssertionError(
            "Every retained primary-antibody discrepancy needs exactly one adjudication"
        )
    if not adjudications["Review complete"].all():
        raise AssertionError("Unfinished antibody review")
    resolved = adjudications.merge(decisions, on="Gene ID", validate="one_to_one")
    if not resolved.Decision.eq(resolved.Category).all():
        raise AssertionError("Evidence-derived outcome differs from the documented adjudication")
    if not set(adjudications["Review ID"]) <= set(reviews["Review ID"]):
        raise AssertionError("Adjudication is missing its source-specific finding")
    if raw_ids & set(decisions.loc[decisions.Category.eq("CTA"), "Gene ID"]):
        raise AssertionError("Unresolved raw normal staining remains in the CTA set")
    genes = base.copy()
    for col in ("Category", "Detail", "Policy", "Reason"):
        genes["Previous " + col.lower()] = base[col]
    for col in decisions:
        genes[col] = decisions[col]
    # Compatibility columns must describe the current findings, not stale history.
    genes["Published normal protein"] = genes["Accepted normal protein"]
    genes["Publication support"] = genes["Reviewed affirmative support"]
    genes["Publication sources"] = genes["Evidence sources"]
    genes["Change"] = "Unchanged"
    changed = genes.Category.ne(base.Category)
    genes.loc[changed, "Change"] = (
        base.loc[changed, "Category"] + " → " + genes.loc[changed, "Category"]
    )
    candidates = genes.loc[genes.Nominated]
    extra = read(BASE / "sequence-gene-reviews.csv.gz")
    extra["Previous policy"], extra["Policy"] = extra.Policy, policy.POLICY_ID
    maxima = read(BASE / "rna-context-maxima.csv.gz")
    high_ids = set(maxima.loc[~maxima.permitted & maxima.maximum.gt(30), "gene_id"])
    all_genes = pd.concat([genes, extra], ignore_index=True)
    if high_ids - set(all_genes.loc[all_genes.Category.eq("Unsafe"), "Gene ID"]):
        raise AssertionError("Measured high normal RNA escaped Unsafe")
    protein_table, sequence_manifest, sequence_names, incomplete, missing = proteins(
        assets, previous_assets, all_genes
    )
    recipe = json.loads((assets / "exclusion-basis.json").read_text())
    recipe["policy"] = policy.POLICY_ID
    dump(recipe, assets / "exclusion-basis.json")
    notes = read(BASE / "follow-up.csv")[
        ["Gene ID", "Review question", "Reconsideration evidence", "Limitations"]
    ]
    new_notes = []
    for gid, rs in reviews.loc[~reviews["Review ID"].str.startswith("legacy-")].groupby("Gene ID"):
        new_notes.append(
            {
                "Gene ID": gid,
                "Review question": "Resolve the accepted normal-protein finding"
                if rs.Status.eq("Accepted").any()
                else "Resolve the material normal-expression/identity conflict",
                "Reconsideration evidence": "Review full original methods, endogenous protein specificity and exact isoforms in the named normal cells. A rescue requires an explicit evidence revision or independently assessed protein/peptide restriction; tumor expression alone is insufficient.",
                "Limitations": " | ".join(sorted(set(rs.Limitations))),
            }
        )
    notes = pd.concat([notes, pd.DataFrame(new_notes)], ignore_index=True).drop_duplicates(
        "Gene ID", keep="last"
    )
    queue = notes.merge(genes, on="Gene ID", validate="one_to_one")
    queue["Review date"] = "2026-10-25"
    queue["Evidence report"] = "docs/audits/cta-policy-20261025/index.md"
    benchmark_ids = set(read(BASE / "important-ctas.csv")["Gene ID"])
    benchmark = genes.loc[genes["Gene ID"].isin(benchmark_ids)]
    frames = {
        "decision-facts.csv.gz": facts,
        "reviewed-findings.csv": reviews,
        "antibody-adjudications.csv": adjudications,
        "gene-results.csv.gz": genes,
        "candidate-results.csv.gz": candidates,
        "changes.csv": genes.loc[changed],
        "follow-up.csv": queue,
        "sequence-gene-reviews.csv.gz": extra,
        "important-ctas.csv": benchmark,
    }
    for group in ("CTA", "Neither", "Unsafe"):
        frames[f"{group.lower()}-genes.csv.gz"] = genes.loc[genes.Category.eq(group)]
    for name, frame in frames.items():
        write(frame, out / name)
    counts = genes.Category.value_counts().to_dict()
    summary = {
        "policy": policy.POLICY_ID,
        "genome": counts,
        "nominated": candidates.Category.value_counts().to_dict(),
        "proteins": protein_table.Group.value_counts().to_dict(),
        "changes": genes.loc[changed]
        .groupby(["Previous category", "Category"])
        .size()
        .rename_axis(["from", "to"])
        .reset_index(name="genes")
        .to_dict("records"),
        "baseline_replayed_exactly": True,
        "antibody_reviews_completed": len(adjudications),
        "antibody_review_outcomes": adjudications.Decision.value_counts().to_dict(),
        "changes_from_local_draft": int(decisions.Category.ne(draft.Category).sum()),
        "followup_genes": len(queue),
        "neither_or_cta_above_30": 0,
        "extra_sequence_loci": len(extra),
        "incomplete_translations": incomplete,
        "genes_without_complete_proteins": missing,
        "literature_review_limit": "32 genes newly read in the initial primary-paper review; 46 have specific notes including earlier work. Search hits are not adjudications.",
        "remaining_primary_ihc_queue": len(
            set(read(LITERATURE / "retained-primary-normal-signals.csv").gene_id)
            & set(genes.loc[genes.Category.eq("CTA"), "Gene ID"])
        ),
    }
    dump(summary, out / "summary.json")
    asset_hashes = {p.name: digest(p) for p in sorted(assets.iterdir()) if p.is_file()}
    dump(
        {
            "release": RELEASE,
            "download_base": RELEASE.replace("/tag/", "/download/") + "/",
            "sha256": asset_hashes,
        },
        out / "protein-assets.json",
    )
    for name in ("exclusion-basis.json", "reference.json"):
        shutil.copyfile(assets / name, out / name)
    m = {
        "schema": policy.SCHEMA,
        "policy": policy.POLICY_ID,
        "baseline_policy": policy.BASE_POLICY_ID,
        "reference": reference,
        "taxonomy_id": 9606,
        "scope": {
            "canonical_genes": 20097,
            "source_nominations": 2540,
            "additional_sequence_loci": len(extra),
            "legacy_cta_api": "Unchanged historical 624-gene selection; not this current policy",
            "other_species": "Anatomy mappings available; numerical policy is human-only",
            "clinical_products": "Product-specific evidence does not confer gene-wide clearance",
            "protein_precedence": ["Unsafe", "Neither", "CTA"],
            "missing_evidence": "Unknown; never zero or proof of absence",
        },
        "measurement_policy": "oncoref.cta-screen.2026-10-23.1",
        "evidence_order": [
            "accepted normal exclusion",
            "unresolved or insufficient",
            "qualified CTA",
        ],
        "rna_admission_ntpm": 5,
        "hard_normal_rna_ceiling": {
            "value": 30,
            "comparison": ">",
            "units_separately": ["TPM", "nTPM"],
        },
        "ihc_columns": [
            "missing_or_uncertain",
            "Approved_negative",
            "Supported_negative",
            "Enhanced_negative",
        ],
        "rna_ceiling_ntpm": {
            "vital": [2, 4, 6, 12],
            "protected": [6, 8, 12, 24],
            "manageable": [10, 20, 30, 60],
        },
        "selectivity_setting": {
            "vital": [0.80, 0.75, 0.70, 0.65],
            "protected": [0.75, 0.70, 0.65, 0.60],
            "manageable": [0.70, 0.65, 0.60, 0.55],
        },
        "bulk_rule": "N <= C AND N <= R*(1-s)/s + 0.1 AND N < R; independent >30 rule also applies",
        "ihc_relaxation": "Complete local routine negative IHC at the weakest confidence, plus reliable reproductive-positive curated profile; not paired specimens or sensitivity calibration",
        "normal_cell_rules": "Pinned native-unit triggers in normal_cell_screen.py; no unit pooling or antibody relaxation of unpaired cell RNA",
        "normal_cell_triggers": {
            "hpa_single_cell": {
                "minimum_ntpm": 5,
                "minimum_reads": 20,
                "comparison": ">=",
                "relative_rule": "At least the reproductive-cell maximum, or missing reproductive coverage/unresolved anatomy",
            },
            "dice_distributions": {
                "minimum_tpm": 5,
                "minimum_reported_values": 3,
                "comparison": ">=",
                "identity_limit": "Distribution values are not stable cross-gene donor identities",
            },
            "hpa_monaco_immune": {
                "minimum": 5,
                "units_separately": ["TPM", "nTPM"],
                "comparison": ">=",
            },
            "cornea": {
                "minimum_mean_cp10k": 1,
                "minimum_positive_cells": 5,
                "minimum_counts": 20,
                "minimum_sample_groups": 2,
                "comparison": ">=",
            },
            "outcome": "Neither unless a stronger exclusion applies; no unpaired IHC relaxation",
        },
        "counts": summary,
        "forbidden_categories": ["Unsafe"],
        "report": "docs/audits/cta-policy-20261025/index.md",
        "report_url": "https://pirl-unc.github.io/oncoref/audits/cta-policy-20261025/",
        "reports_at_release": RELEASE.replace("/releases/tag/", "/tree/")
        + "/docs/audits/cta-policy-20261025",
        "decision_artifact_sha256": {
            name: digest(out / name)
            for name in [
                "decision-facts.csv.gz",
                "reviewed-findings.csv",
                "candidate-results.csv.gz",
                "changes.csv",
            ]
        },
        "sequence_artifacts": {
            "release": RELEASE,
            "download_base": RELEASE.replace("/tag/", "/download/") + "/",
            "sha256": asset_hashes,
        },
        "data_sha256": {
            "target-gene-sets.csv.gz": digest(out / "gene-results.csv.gz"),
            "cta-policy-findings.csv": digest(ROOT / "oncoref/data/cta-policy-findings.csv"),
            "cta-antibody-adjudications.csv": digest(
                ROOT / "oncoref/data/cta-antibody-adjudications.csv"
            ),
            "sequence-gene-rna-reviews.csv.gz": digest(out / "sequence-gene-reviews.csv.gz"),
            "cta-followup-reviews.csv": digest(out / "follow-up.csv"),
        },
    }
    dump(m, out / "policy.json")
    changes_table = markdown(
        benchmark[["Gene", "Category", "Detail"]], ["Gene", "Category", "Detail"]
    )
    (out / "index.md").write_text(f"""# CTA definitions — 25 October 2026

**{counts["CTA"]} CTA genes**, {counts["Neither"]} Neither and {counts["Unsafe"]} Unsafe across 20,097 canonical human genes.
Among 2,540 mapped source nominations: **{summary["nominated"]["CTA"]} CTA / {summary["nominated"]["Neither"]} Neither / {summary["nominated"]["Unsafe"]} Unsafe**.

[Plain-language rules and complete threshold table](../../target-sets.md) · [All candidates](candidate-results.csv.gz) · [All genes](gene-results.csv.gz) · [CTA](cta-genes.csv.gz) · [Neither](neither-genes.csv.gz) · [Unsafe](unsafe-genes.csv.gz)

The updated policy applies reviewed literature as biological findings, not gene-name overrides. The unchanged measurement facts plus the previously accepted findings reproduce all 20,097 old categories exactly. New reviewed findings change {int(changed.sum())} assignments from the released 399-gene CTA set (6 to Unsafe, 180 to Neither). RNA thresholds and tissue mappings are unchanged.

## Named targets — retained and excluded

{changes_table}

Every changed gene has a [review record](reviewed-findings.csv) with source, species, reference, assay, identity/specificity judgement, review depth, limitations and linked issue. Abstract-level decisions are marked. Unsafe means exclusion from unrestricted gene-level CTA targeting, not demonstrated clinical injury. Neither conflicts have a completed curation decision; biological uncertainty remains. They need evidence that resolves the conflict before reconsideration.

## Evidence and follow-up

[Decision facts](decision-facts.csv.gz) · [Reviewed literature](reviewed-findings.csv) · [Changed rows](changes.csv) · [Follow-up questions](follow-up.csv) · [Named CTA checks](important-ctas.csv) · [Literature search and review coverage](../cta-literature-20261025/index.md)

CT55 remains CTA through affirmative reproductive and tumor-recognition evidence below the RNA floor; normal gates still apply. CTAG1A/B and XAGE1A/B remain CTA. CTAG2 remains Neither / Variable and FOXR2 Neither / Antibody. MAGEA4 is now Neither / Conflict; PRAME is Unsafe because an independent human adrenal-protein observation is accepted. This does not relabel their selected therapeutic products as clinically unsafe.

[All 178 completed antibody adjudications, evidence and reconsideration requirements](../cta-adjudication-20261025/index.md). PLCZ1's pancreatic paper measures its RNA, while its immunostaining detects CFTR; the erroneous protein interpretation is rejected and the qualitative RNA finding is held separately.

The literature audit is incomplete: 32 genes received new primary-paper assessment, with 46 specific notes including prior reviews. All **178** raw-antibody discrepancies remaining in the local 391-gene CTA draft received documented adjudications: **177 Neither and 1 Unsafe**. The two other genes in the original 180-gene queue were already addressed by the eight-gene review. **{summary["remaining_primary_ihc_queue"]} unresolved raw-antibody cases remain retained.** Review completion does not mean the molecular identity of every stain was resolved. The retained CTA set remains a research selection, not validated clinical safety. No new continuous breadth gate or calibrated RNA uncertainty model is claimed.

No remaining CTA or Neither has an audited non-permitted/unsited TPM or nTPM observation >30. Missing observations remain gaps; this is not proof of absence in unmeasured cells. Original source versions, native units and donor limitations are preserved in the [RNA audit](../cta-neither-review-20261024/index.md).

## Exact proteins

Complete exact-sequence groups: **{summary["proteins"]["CTA"]:,} CTA / {summary["proteins"]["Neither"]:,} Neither / {summary["proteins"]["Unsafe"]:,} Unsafe**. Priority for a protein encoded by multiple genes is Unsafe → Neither → CTA. The additional 1,721 reference loci retain their measured decisions, reissued under this policy. Neither is excluded from the default forbidden union.

[Sequence recipe](exclusion-basis.json) · [Reference](reference.json) · [Release asset inventory and SHA-256](protein-assets.json) · [Additional loci](sequence-gene-reviews.csv.gz)

Use the matching **[v1.20.0 release assets]({RELEASE})** and their recorded checksums. Earlier FASTAs retain their dated policies and must not be combined with the new gene list. The historical 624-gene CTA API is unchanged; downstream applications must explicitly choose the current `target_sets` policy.

## Reproduction

```sh
.venv/bin/python scripts/build_cta_policy.py --output /tmp/cta-policy-report --assets /tmp/cta-policy-proteins --sequence-input /path/to/pinned-original-sequence-assets
```

Use empty output directories. Add `--publish-dataset` only to update the repository's packaged current data after all invariants pass. The builder verifies pinned input and original sequence hashes, replays the baseline, evaluates every reason, and checks the >30 invariant before publishing.

[Policy](policy.json) · [Counts](summary.json) · [Input/code/output hashes](run-manifest.json). Historical reports remain frozen.
""")
    inputs = [
        BASE / "run-manifest.json",
        BASE / "gene-results.csv.gz",
        BASE / "rna-context-maxima.csv.gz",
        BASE / "sequence-gene-reviews.csv.gz",
        BASE / "follow-up.csv",
        BASE / "important-ctas.csv",
        EXPRESSION / "run-manifest.json",
        EXPRESSION / "gene-categories.csv.gz",
        LITERATURE / "retained-primary-normal-signals.csv",
        ROOT / "oncoref/data/cta-policy-findings.csv",
        ROOT / "oncoref/data/cta-antibody-adjudications.csv",
        ADJUDICATION / "baseline-gene-results.csv.gz",
        ADJUDICATION / "primary-observations.csv.gz",
        ROOT / "oncoref/data/cta-external-safety-reviews.csv",
        ROOT / "oncoref/cta_policy.py",
        Path(__file__),
        ROOT / "scripts/target_sets_report.py",
        ROOT / "oncoref/target_sets.py",
        ROOT / "docs/audits/normal-sequence-sets-20261021/run-manifest.json",
    ]
    dump(
        {
            "policy": policy.POLICY_ID,
            "implementation_and_inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
            "sequence_inputs": {n: sequence_manifest["outputs"][n] for n in sequence_names},
            "outputs": {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()},
            "forbidden_categories": ["Unsafe"],
        },
        out / "run-manifest.json",
    )
    if publish:
        for source, destination in [
            ("gene-results.csv.gz", "target-gene-sets.csv.gz"),
            ("follow-up.csv", "cta-followup-reviews.csv"),
            ("sequence-gene-reviews.csv.gz", "sequence-gene-rna-reviews.csv.gz"),
            ("policy.json", "cta-policy.json"),
        ]:
            shutil.copyfile(out / source, ROOT / "oncoref/data" / destination)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--sequence-input", type=Path, required=True)
    parser.add_argument("--publish-dataset", action="store_true")
    args = parser.parse_args()
    build(
        args.output.resolve(),
        args.assets.resolve(),
        args.sequence_input.resolve(),
        args.publish_dataset,
    )
