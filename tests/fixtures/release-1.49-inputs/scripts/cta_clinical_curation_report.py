#!/usr/bin/env python3
"""Join publication intake, clinical target evidence and genome-wide expression groups."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_recalibration_report import markdown_table  # noqa: E402
from cta_scope_report import sha256, write_csv  # noqa: E402

from oncoref import figure_style  # noqa: E402


def read_verified(directory, name):
    manifest = json.loads((directory / "run-manifest.json").read_text())
    assert sha256(directory / name) == manifest["outputs"][name]
    return pd.read_csv(directory / name, keep_default_na=False)


def clinical_disposition(row):
    if row.external_safety_hold:
        return "external_safety_hold"
    if row.retained:
        return "retained"
    reasons = set(row.failure_reasons.split(";"))
    if "positive_ihc_outside_permitted" in reasons:
        return "normal_protein_counterevidence"
    if "rna_above_tissue_ceiling" in reasons:
        return "rna_threshold_review"
    if "reproductive_selectivity_below_threshold" in reasons:
        return "selectivity_reference_review"
    return "coverage_or_curation_review"


def build(out, categories, clinical_inputs):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    baseline = ROOT / "docs/audits/cta-recalibration-20261021"
    original = read_verified(baseline, "candidate-comparison.csv")
    intake = read_verified(baseline, "intake-decisions.csv")
    genes = read_verified(categories, "gene-categories.csv.gz").set_index("Ensembl_Gene_ID")
    clinical_path = ROOT / "oncoref/data/cta-trial-evidence.csv"
    safety_path = ROOT / "oncoref/data/cta-external-safety-reviews.csv"
    clinical = pd.read_csv(clinical_path, keep_default_na=False)
    safety = pd.read_csv(safety_path, keep_default_na=False)
    assert clinical.Ensembl_Gene_ID.is_unique and clinical.Ensembl_Gene_ID.ne("").all()
    assert safety.Ensembl_Gene_ID.is_unique
    current = original.copy()
    added = []
    for row in clinical.itertuples():
        if row.Ensembl_Gene_ID in set(current.Ensembl_Gene_ID):
            continue
        gene = genes.loc[row.Ensembl_Gene_ID]
        reason = ";".join(filter(None, [gene.failure_reasons, gene.quality_flags]))
        if not gene.reproductive_positive:
            reason = ";".join(filter(None, [reason, "insufficient_positive_reproductive_evidence"]))
        added.append(
            {
                **gene.to_dict(),
                "Ensembl_Gene_ID": row.Ensembl_Gene_ID,
                "evidence_pass": gene.category == "Reproductive",
                "family_eligible": True,
                "nomination_review_eligible": True,
                "nomination_review_reason": "Primary clinical intervention nomination; general expression gates apply",
                "retained": gene.category == "Reproductive",
                "failure_reasons": reason,
                "source_dois": "",
                "clinical_source_url": row.source_url,
            }
        )
    current = pd.concat([current, pd.DataFrame(added)], ignore_index=True).fillna("")
    current["pre_safety_retained"] = current.retained.astype(bool)
    current["external_safety_hold"] = current.Ensembl_Gene_ID.isin(
        safety.loc[safety.action.eq("hold"), "Ensembl_Gene_ID"]
    )
    current["retained"] = current.pre_safety_retained & ~current.external_safety_hold
    for index in current.index[current.external_safety_hold]:
        current.at[index, "failure_reasons"] = ";".join(
            filter(None, [current.at[index, "failure_reasons"], "external_safety_hold"])
        )
    current["expression_category"] = current.Ensembl_Gene_ID.map(genes.category)
    current["clinical_evidence_available"] = current.Ensembl_Gene_ID.isin(clinical.Ensembl_Gene_ID)
    clinical = clinical.merge(
        current, on=["Symbol", "Ensembl_Gene_ID"], how="left", validate="one_to_one"
    )
    assert clinical.retained.notna().all(), "Every clinical benchmark locus must have a decision"
    clinical["clinical_disposition"] = [clinical_disposition(r) for r in clinical.itertuples()]
    clinical = clinical.merge(
        safety[["Ensembl_Gene_ID", "action", "rationale", "source_url"]].rename(
            columns={"source_url": "safety_source_url", "rationale": "safety_rationale"}
        ),
        on="Ensembl_Gene_ID",
        how="left",
        validate="one_to_one",
    ).fillna("")
    new_ids = set(current.Ensembl_Gene_ID) - set(intake.Ensembl_Gene_ID)
    new_intake = [
        {
            "intake_identity": gid,
            "Symbol": genes.loc[gid, "Symbol"],
            "Ensembl_Gene_ID": gid,
            "identity_status": "canonical_mapped",
            "screened": True,
            "protein_candidate_eligible": True,
        }
        for gid in sorted(new_ids)
    ]
    intake = pd.concat([intake, pd.DataFrame(new_intake)], ignore_index=True).fillna("")
    decisions = current.set_index("Ensembl_Gene_ID")
    intake["retained"] = intake.Ensembl_Gene_ID.map(decisions.retained).fillna(False).astype(bool)
    intake["failure_reasons"] = [
        decisions.loc[gid, "failure_reasons"] if gid in decisions.index else old
        for gid, old in zip(intake.Ensembl_Gene_ID, intake.failure_reasons)
    ]
    intake["screen_status"] = [
        "retained" if r.retained else "excluded" if r.screened else r.screen_status
        for r in intake.itertuples()
    ]
    assert intake.intake_identity.is_unique
    failures = current.failure_reasons.str.split(";").map(set)
    mask = pd.Series(True, index=current.index)
    funnel = [
        {"stage": "Publication and clinical intake identities", "remaining": len(intake)},
        {"stage": "Mapped coding candidates", "remaining": len(current)},
    ]
    for label, rejected in [
        ("Positive reproductive evidence", {"insufficient_positive_reproductive_evidence"}),
        (
            "Resolved evidence and routine RNA coverage",
            {
                "incomplete_routine_rna",
                "unusable_measurement",
                "unresolved_positive_anatomy",
                "reviewed_unreliable_ihc",
                "unresolved_reproductive_selectivity",
            },
        ),
        ("No positive outside IHC", {"positive_ihc_outside_permitted"}),
        ("Absolute RNA ceilings", {"rna_above_tissue_ceiling"}),
        ("Pairwise selectivity", {"reproductive_selectivity_below_threshold"}),
        (
            "Family and nomination review",
            {"family_exclusion", "existing_nomination_or_manual_review"},
        ),
        ("External normal-tissue safety evidence", {"external_safety_hold"}),
    ]:
        mask &= ~failures.map(lambda value, rejected=rejected: bool(value & rejected))
        funnel.append({"stage": label, "remaining": int(mask.sum())})
    assert mask.equals(current.retained)
    funnel = pd.DataFrame(funnel)
    for filename, frame in {
        "candidate-decisions.csv.gz": current,
        "clinical-target-audit.csv": clinical,
        "intake-decisions.csv.gz": intake,
        "clinical-target-evidence.csv": pd.read_csv(clinical_path),
        "external-safety-reviews.csv": safety,
        "funnel.csv": funnel,
    }.items():
        write_csv(frame, out / filename)
    # Preserve scientific registry fields, excluding investigator contact/location lists.
    if (clinical_inputs / "trial-registry-snapshot.json.gz").exists():
        studies = json.loads(
            gzip.decompress((clinical_inputs / "trial-registry-snapshot.json.gz").read_bytes())
        )
    else:
        records = json.loads((clinical_inputs / "studies.json").read_text())
        records.extend(json.loads(p.read_text()) for p in sorted(clinical_inputs.glob("NCT*.json")))
        selected = {}
        for record in records:
            p = record["protocolSection"]
            selected[p["identificationModule"]["nctId"]] = {
                k: p[k]
                for k in [
                    "identificationModule",
                    "statusModule",
                    "descriptionModule",
                    "armsInterventionsModule",
                    "conditionsModule",
                    "designModule",
                ]
                if k in p
            }
        studies = [selected[k] for k in sorted(selected)]
    (out / "trial-registry-snapshot.json.gz").write_bytes(
        gzip.compress(json.dumps(studies, sort_keys=True).encode(), mtime=0)
    )
    query = json.loads((clinical_inputs / "query.json").read_text())
    (out / "query.json").write_text(json.dumps(query, indent=2) + "\n")
    for symbol in ["BAGE", "GAGE2B"]:
        (out / f"{symbol}-hgnc.json").write_bytes(
            (clinical_inputs / f"{symbol}-hgnc.json").read_bytes()
        )
    registry = pd.DataFrame(
        [
            {
                "nct_id": p["identificationModule"]["nctId"],
                "title": p["identificationModule"]["briefTitle"],
                "status": p["statusModule"]["overallStatus"],
                "source_url": "https://clinicaltrials.gov/study/"
                + p["identificationModule"]["nctId"],
            }
            for p in studies
        ]
    )
    write_csv(registry, out / "registry-search-results.csv")
    render(out, current, clinical, intake, funnel, genes, len(registry))
    inputs = [
        clinical_path,
        safety_path,
        baseline / "candidate-comparison.csv",
        baseline / "intake-decisions.csv",
        categories / "gene-categories.csv.gz",
        categories / "run-manifest.json",
    ]
    manifest = {
        "schema": "oncoref.cta-clinical-curation.v1",
        "reviewed_date": "2026-10-21",
        "inputs": {str(p.resolve().relative_to(ROOT)): sha256(p) for p in inputs},
        "implementation_sha256": {
            "scripts/cta_clinical_curation_report.py": sha256(Path(__file__))
        },
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(funnel.to_string(index=False))
    print(clinical[["Symbol", "clinical_disposition"]].to_string(index=False))


def render(out, current, clinical, intake, funnel, genes, n_registry):
    figure_style.apply()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.barh(funnel.stage, funnel.remaining, color="#377eb8")
    ax.invert_yaxis()
    for i, value in enumerate(funnel.remaining):
        ax.text(value + 35, i, str(value), va="center")
    ax.set_xlim(0, funnel.remaining.max() * 1.12)
    ax.set_xlabel("Distinct identities / candidates")
    ax.set_title("CTA curation: publication and clinical intake through final review")
    fig.subplots_adjust(left=0.40, bottom=0.13, right=0.98, top=0.90)
    for suffix in ["png", "svg"]:
        figure_style.save(fig, out / f"curation-funnel.{suffix}", keep_titles=True)
    plt.close(fig)
    clinical_table = clinical[["Symbol", "evidence_kind", "clinical_disposition"]].sort_values(
        "Symbol"
    )
    categories = genes.groupby(["category", "other_subgroup"]).size().reset_index(name="genes")
    text = f"""# CTA curation, clinical coverage and expression categories — 21 October 2026

**{len(intake):,} intake identities → {len(current):,} mapped coding candidates → {int(current.retained.sum()):,} retained after external safety review.** The source union contains the 13 imported publication sources plus the explicitly reviewed clinical-target evidence. No clinical benchmark target is silently omitted.

## General rules

1. Preserve every source nomination and clinical target upstream, including unresolved official identities. Trial or therapeutic use supplies nomination evidence; it does not override normal-tissue counterevidence.
2. Require reproductive RNA >0 or positive reproductive IHC, complete routine RNA, the approved per-tissue RNA ceilings and selectivity minima, and no positive IHC outside Reproductive/Disposable. [Exact grid](../../normal-tissue-scopes.md).
3. Keep family/nomination reviews explicit. An external normal-tissue injury signal supported by target-expression evidence places a clinical-prioritization hold. Receptor-specific cross-reactivity is recorded separately from gene-wide normal expression.
4. Distinguish **normal protein counterevidence**, **RNA-threshold review**, **selectivity/reference review**, **coverage gaps**, and **external safety holds**. A failed research threshold is not proof of clinical harm. Every clinical target stays visible for review even when it fails automatic retention.

The prior [recalibration comparison](../cta-recalibration-20261021/index.md) retains 649 under its HPA and nomination rules. This report adds KIF20A and CLDN6 from clinical evidence and applies the external safety review. It leaves the approved numerical grid unchanged. Frozen public core/extended panels remain separate.

## Current funnel

{markdown_table(funnel)}

![Curation funnel](curation-funnel.png)

## Clinical target coverage

{markdown_table(clinical_table)}

The audit covers **{len(clinical)} reviewed loci**, including related oncofetal/TAA comparators and one cross-recognition safety context; these are not {len(clinical)} independent approved therapies. Historical antigen labels and shared paralogs remain explicit. The broad registry search returned {n_registry} records after targeted record supplementation; search hits include unrelated/background mentions and are not automatically treated as target interventions. This is a dated, reproducible coverage audit, not a guarantee that every trial worldwide has been indexed.

Examples of primary anchors: [MAGE-A4/Tecelra approval](https://www.cancer.gov/about-cancer/treatment/drugs/afamitresgene-autoleucel), [PRAME/IMA203](https://www.nature.com/articles/s41591-025-03650-6), [KK-LC-1](https://clinicaltrials.gov/study/NCT05035407), [AKAP4/SP17/ropporin/Span-xb trial](https://clinicaltrials.gov/study/NCT02224599), [CLDN6 phase 1](https://www.nature.com/articles/s41591-023-02612-0). Every row has its own source and resolution limits in the downloadable audit.

### Why some clinically studied targets are not automatically retained

- **SPAG9:** substantial CNS RNA and positive normal neuronal IHC, including High staining in cerebral-cortex neurons. Its [institutional DC-vaccine trial report](https://cancerinstitutewia.in/wp-content/uploads/2024/08/Annual-Report-2022-2023-Final-Draft-for-Print-proofread_organized-for-website.pdf) establishes clinical use, not safety for every targeting modality.
- **SP17/SPA17:** positive bronchial ciliated-cell and other normal IHC, plus CNS RNA. These are measured normal-tissue signals under the current gate.
- **XAGE1B and LY6K:** RNA-only threshold failures require review; they are not labelled clinically unsafe. XAGE1B has liver 1.5/lung 1.2 nTPM without usable local IHC relaxation. Its historical antigen/isoform label also requires sequence-level paralog review.
- **CLDN6:** oncofetal clinical evidence is included even though the gene was missing from the classic CTA intake. The adult HPA panel does not measure its fetal reference. A low adult reproductive maximum can make pairwise selectivity fail; this is specifically a developmental-reference review, not a demonstrated clinical-toxicity claim.
- **MAGE-A12:** normal-brain RNA and cross-recognition in a trial with severe neurological injury support an explicit precautionary hold, even though its HPA screen passes. The primary report cannot uniquely attribute injury to this paralog. [Morgan et al.](https://pubmed.ncbi.nlm.nih.gov/23377668/)
- **MAGE-A3:** the titin cross-reactivity tragedy requires receptor-specific assessment. It does not by itself establish cardiac MAGE-A3 expression or a gene-wide exclusion. [Linette et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3743463/)

### BAGE, GAGE2B and SPAG9

All three are present in intake. SPAG9 is mapped and screened. HGNC identifies **BAGE (HGNC:942; NCBI:574; RefSeq NM_001187)** and **GAGE2B (HGNC:31957; NCBI:645037; RefSeq NM_001098411)** as not on the reference assembly and supplies no Ensembl gene ID in the retrieved records. Neither has a row under its current/legacy target label or the historical GAGE2B Ensembl alias in the pinned HPA v23 RNA/IHC tables. They remain unassessed; no paralog is substituted and no absent measurement becomes zero. [HGNC BAGE](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:942), [HGNC GAGE2B](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:31957).

## Reproductive, Normal and Other

{markdown_table(categories)}

These expression groups span all {len(genes):,} canonical coding genes, including non-CTAs. **Other = Disposable ∪ Silent**, under the approved absolute category/IHC limits and with no reproductive-positive requirement or ratio. Silent means below-threshold expression; missing routine RNA remains Unassessed. Missing IHC uses baseline limits and stays explicitly unassessed for protein absence.

For future subsequence analysis, a reproductive screening label must not erase a real normal-expression observation or a shared normal paralog. The [genome-wide export](../normal-expression-categories-20261021/index.md) therefore includes a separate conservative normal-reference flag that can overlap the disjoint screening categories. Exact sequence, transcript and protein occurrences still need to be compared.

## Downloads and reproducibility

- [Complete updated intake](intake-decisions.csv.gz) · [All candidate decisions](candidate-decisions.csv.gz)
- [Clinical target decisions and exact exclusions](clinical-target-audit.csv) · [Primary clinical evidence](clinical-target-evidence.csv)
- [External safety reviews](external-safety-reviews.csv) · [Funnel](funnel.csv)
- [Registry search results](registry-search-results.csv) · [Scientific registry snapshot](trial-registry-snapshot.json.gz) · [Search query](query.json)
- [BAGE HGNC snapshot](BAGE-hgnc.json) · [GAGE2B HGNC snapshot](GAGE2B-hgnc.json) · [Checksums](run-manifest.json)
- [Source publications and their per-source counts](../cta-recalibration-20261021/index.md) · [Full tissue/species mappings](../normal-tissue-recalibration-20261021/index.md)

Reproduce with `python scripts/cta_clinical_curation_report.py --out <empty-directory> --categories <normal-expression-report-directory> --clinical-inputs <this-report-directory>`. The data and numerical policies are research screening evidence, not validated clinical safety limits.
"""
    (out / "index.md").write_text(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--categories", type=Path, required=True)
    parser.add_argument("--clinical-inputs", type=Path, required=True)
    args = parser.parse_args()
    build(args.out, args.categories, args.clinical_inputs)
