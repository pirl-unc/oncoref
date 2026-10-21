#!/usr/bin/env python3
"""Rebuild the complete published intake and the recalibrated human CTA screen."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_permitted_report import (  # noqa: E402
    EXPANDED,
    TARGETS,
    annotate_archive,
    evaluate_gene,
    sha256,
    write_csv,
)
from cta_scope_report import annotate_measurements  # noqa: E402
from cta_tissue_tier_report import checked_input  # noqa: E402

from oncoref import cta_regen, cta_sources, figure_style, hpa, reference_data  # noqa: E402
from oncoref.cta_normal_screen import (  # noqa: E402
    adaptive_normal_rna_limits,
    screening_nomination_review,
)
from oncoref.normal_tissues import (  # noqa: E402
    normal_tissue_screening_policy,
    normal_tissue_taxonomy,
)

PREVIOUS = ROOT / "docs/audits/cta-selectivity-20261020"
OBSERVATIONS = ROOT / "docs/audits/cta-permitted-tissues-20261019"
TIERS = ("vital", "protected", "manageable")
SCENARIOS = (
    "admission_only",
    "selectivity_revised",
    "grid_revised",
    "current_adaptive",
    "current_fixed",
)
NAMED = set(TARGETS) | {
    "BAGE",
    "GAGE2B",
    "SPAG9",
    "SAGE1",
    "GAGE1",
    "SSX4",
    "SSX5",
    "MAGEA6",
    "MAGEA12",
}


def with_grid(audit, rna_policy, selectivity_policy, *, fixed=False):
    result = audit.copy()
    for index, row in result.iterrows():
        if row.screening_tier not in TIERS:
            continue
        tier, reliability = row.screening_tier, row.ihc_reliability_used
        result.at[index, "rna_ceiling_ntpm"] = (
            rna_policy[tier]["maximum_rna"]
            if fixed
            else rna_policy["negative_ihc_rna_limits"][reliability][tier]
        )
        result.at[index, "minimum_reproductive_selectivity"] = (
            selectivity_policy[tier]["minimum_reproductive_selectivity"]
            if fixed
            else selectivity_policy["negative_ihc_minimum_reproductive_selectivity"][reliability][
                tier
            ]
        )
    result["selectivity_margin"] = (
        result.reproductive_selectivity - result.minimum_reproductive_selectivity
    )
    if fixed:
        result["ihc_reliability_used"] = "Baseline (adaptation disabled)"
    return result


def build(out):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory; historical reports are immutable")
    out.mkdir(parents=True, exist_ok=True)
    previous_path = checked_input(PREVIOUS, "candidate-comparison.csv")
    observation_path = checked_input(OBSERVATIONS, "normal-observations.csv.gz")
    previous = pd.read_csv(previous_path, keep_default_na=False).set_index("Ensembl_Gene_ID")
    old_manifest = json.loads((PREVIOUS / "run-manifest.json").read_text())
    old_policy = old_manifest["policy"]
    policy = normal_tissue_screening_policy("ihc_adaptive_cta")
    expected = {
        name: set(map(tuple, pairs)) for name, pairs in old_manifest["routine_contexts"].items()
    }
    membership = cta_sources.reviewed_publication_membership().merge(
        cta_sources.publication_sources()[
            ["source_tag", "doi", "citation", "source_url", "source_table", "input_sha256"]
        ],
        on="source_tag",
        how="left",
        validate="many_to_one",
    )
    assert membership.doi.notna().all()
    coding = membership.loc[membership.protein_candidate_eligible]
    ids = set(coding.Ensembl_Gene_ID)
    added = ids - set(previous.index)
    assert added == {"ENSG00000008294"}, "Review new intake additions before changing this snapshot"
    observations = pd.read_csv(observation_path, keep_default_na=False, low_memory=False)
    receipts, extras = {}, []
    for modality, name in [("bulk_rna", "hpa_rna_consensus"), ("ihc", "hpa_normal_tissue")]:
        receipt = reference_data.provenance(name, "v23")
        assert sha256(receipt["path"]) == receipt["sha256"]
        raw = hpa._read_hpa(name, "v23")
        extras.append(annotate_measurements(raw, modality, added))
        receipts[name] = {k: receipt[k] for k in ["url", "sha256", "bytes", "version"]}
    extra = pd.concat(extras, ignore_index=True)
    observations = annotate_archive(pd.concat([observations, extra], ignore_index=True))
    groups = dict(iter(observations.groupby("Ensembl_Gene_ID")))
    decisions, audits = [], []
    for number, gid in enumerate(sorted(ids), 1):
        if number == 1 or number % 250 == 0:
            print(f"Evaluating {number}/{len(ids)}", flush=True)
        sources = coding.loc[coding.Ensembl_Gene_ID.eq(gid)]
        symbol = sources.Symbol.iloc[0]
        old = (
            previous.loc[gid]
            if gid in previous.index
            else pd.Series(
                {
                    "Symbol": symbol,
                    "family_eligible": True,
                    "nomination_review_eligible": True,
                    "nomination_review_reason": "Primary cancer-associated SPAG9 nomination; current expression gates apply",
                    "ihc_coverage": "unavailable",
                    "retained": False,
                }
            )
        )
        reviewed = screening_nomination_review(
            gid, old.nomination_review_eligible, old.nomination_review_reason
        )
        local = groups.get(gid, observations.iloc[:0])
        unreliable = old.ihc_coverage == "unreliable" or gid in cta_regen.cross_reactive_ihc()
        audit = adaptive_normal_rna_limits(
            local, routine_ihc_contexts=expected["ihc"], reviewed_unreliable=unreliable
        )
        audit["Symbol"] = symbol
        audits.append(audit)
        grids = {
            "admission_only": with_grid(audit, old_policy, old_policy),
            "selectivity_revised": with_grid(audit, old_policy, policy),
            "grid_revised": audit,
            "current_adaptive": audit,
            "current_fixed": with_grid(audit, policy, policy, fixed=True),
        }
        for scenario, active in grids.items():
            current_review = scenario in {"current_adaptive", "current_fixed"}
            eligible = (
                reviewed["nomination_review_eligible"]
                if current_review
                else bool(old.nomination_review_eligible)
            )
            result = evaluate_gene(
                local,
                scope=EXPANDED,
                vital_limit=policy["vital"]["maximum_rna"],
                protected_limit=policy["protected"]["maximum_rna"],
                manageable_limit=policy["manageable"]["maximum_rna"],
                expected=expected,
                unreliable=unreliable,
                rna_limits=active.rna_ceiling_ntpm,
                selectivity_audit=active,
            )
            result.update(
                {
                    **reviewed,
                    "Ensembl_Gene_ID": gid,
                    "Symbol": symbol,
                    "scenario": scenario,
                    "family_eligible": bool(old.family_eligible),
                    "nomination_review_eligible": eligible,
                    "nomination_review_reason": reviewed["nomination_review_reason"]
                    if current_review
                    else old.nomination_review_reason,
                    "screen_review_applied": current_review
                    and reviewed["screen_review_superseded_holdback"],
                    "retained": bool(result["evidence_pass"] and old.family_eligible and eligible),
                    "previous_retained": bool(old.retained),
                    "source_dois": ";".join(sorted(set(sources.doi))),
                    "source_publications": ";".join(sorted(set(sources.citation))),
                }
            )
            for allowed, reason in [
                (old.family_eligible, "family_exclusion"),
                (eligible, "existing_nomination_or_manual_review"),
            ]:
                if not allowed:
                    result["failure_reasons"] = ";".join(
                        filter(None, [result["failure_reasons"], reason])
                    )
            decisions.append(result)
    decisions = pd.DataFrame(decisions)
    frames = {
        name: frame.set_index("Ensembl_Gene_ID") for name, frame in decisions.groupby("scenario")
    }
    chain = [previous, *[frames[s] for s in SCENARIOS[:-1]]]
    for small, large in itertools.pairwise(chain):
        assert set(small.index[small.retained]) <= set(large.index[large.retained])
    assert set(frames["current_fixed"].index[frames["current_fixed"].retained]) <= set(
        frames["current_adaptive"].index[frames["current_adaptive"].retained]
    )
    current = frames["current_adaptive"].copy()
    for name in SCENARIOS:
        current[name + "_retained"] = frames[name].retained
    current["change"] = [
        "gained" if now and not before else "unchanged"
        for now, before in zip(current.retained, current.previous_retained)
    ]
    current["gain_stage"] = [
        next((name for name in SCENARIOS[:-1] if frames[name].loc[gid, "retained"]), "")
        if row.retained and not row.previous_retained
        else ""
        for gid, row in current.iterrows()
    ]
    summary = []
    for name, frame in [("previous_released", previous), *[(s, frames[s]) for s in SCENARIOS]]:
        kept = frame.loc[frame.retained]
        summary.append(
            {
                "scenario": name,
                "candidates": len(frame),
                "evidence_pass": int(frame.evidence_pass.sum()),
                "retained": len(kept),
                "without_ihc": int(kept.ihc_coverage.eq("unavailable").sum()),
                "partial_ihc": int(kept.ihc_coverage.eq("partial").sum()),
                "complete_routine_ihc": int(kept.ihc_coverage.eq("routine_complete").sum()),
            }
        )
    summary = pd.DataFrame(summary)
    intake = []
    for identity, group in membership.groupby("intake_identity", sort=True):
        first = group.iloc[0]
        gid = first.Ensembl_Gene_ID
        screened = gid in current.index
        row = current.loc[gid] if screened else None
        intake.append(
            {
                "intake_identity": identity,
                "Symbol": first.Symbol or first.authority_symbol or first.source_symbol,
                "Ensembl_Gene_ID": gid,
                "identity_status": first.identity_status,
                "protein_candidate_eligible": bool(group.protein_candidate_eligible.any()),
                "screened": screened,
                "retained": bool(row.retained) if screened else False,
                "screen_status": ("retained" if row.retained else "excluded")
                if screened
                else "unassessed_no_canonical_mapping"
                if not gid
                else "outside_protein_coding_screen",
                "failure_reasons": row.failure_reasons
                if screened
                else "missing_canonical_mapping"
                if not gid
                else "not_protein_coding",
                "source_dois": ";".join(sorted(set(group.doi))),
                "authority_url": first.authority_url,
                "identity_rationale": first.rationale,
                "assay_caveat": first.assay_caveat,
            }
        )
    intake = pd.DataFrame(intake)
    funnel = [{"stage": "Published source identities", "remaining": len(intake)}]
    funnel.append({"stage": "Mapped protein-coding candidates", "remaining": len(current)})
    failures = current.failure_reasons.str.split(";").map(set)
    mask = pd.Series(True, index=current.index)
    for label, rejected in [
        ("Positive reproductive evidence", {"insufficient_positive_reproductive_evidence"}),
        (
            "Usable resolved evidence and routine RNA coverage",
            {
                "incomplete_routine_rna",
                "unusable_measurement",
                "unresolved_positive_anatomy",
                "reviewed_unreliable_ihc",
                "unresolved_reproductive_selectivity",
            },
        ),
        ("No positive IHC outside permitted tissues", {"positive_ihc_outside_permitted"}),
        ("Absolute RNA ceilings", {"rna_above_tissue_ceiling"}),
        ("Pairwise selectivity", {"reproductive_selectivity_below_threshold"}),
        ("Family eligibility", {"family_exclusion"}),
        ("Nomination and current review eligibility", {"existing_nomination_or_manual_review"}),
    ]:
        mask &= ~failures.map(lambda reasons, rejected=rejected: bool(reasons & rejected))
        funnel.append({"stage": label, "remaining": int(mask.sum())})
    assert mask.equals(current.retained)
    funnel = pd.DataFrame(funnel)
    source_counts = []
    for doi, group in membership.groupby("doi"):
        subset = current.loc[current.index.isin(group.Ensembl_Gene_ID)]
        source_counts.append(
            {
                "doi": doi,
                "citation": group.citation.iloc[0],
                "intake_identities": group.intake_identity.nunique(),
                "canonical_mapped": group.loc[
                    group.Ensembl_Gene_ID.ne(""), "Ensembl_Gene_ID"
                ].nunique(),
                "coding_candidates": len(subset),
                "evidence_pass": int(subset.evidence_pass.sum()),
                "retained": int(subset.retained.sum()),
            }
        )
    source_counts = pd.DataFrame(source_counts)
    outputs = {
        "candidate-decisions.csv": decisions,
        "candidate-comparison.csv": current.reset_index(),
        "membership-changes.csv": current.loc[current.change.eq("gained")].reset_index(),
        "notable-targets.csv": current.loc[current.Symbol.isin(NAMED)].reset_index(),
        "intake-decisions.csv": intake,
        "intake-target-audit.csv": intake.loc[intake.Symbol.isin(["BAGE", "GAGE2B", "SPAG9"])],
        "source-membership.csv.gz": membership,
        "publication-sources.csv": cta_sources.publication_sources(),
        "source-counts.csv": source_counts,
        "rna-selectivity-audit.csv.gz": pd.concat(audits).reset_index(),
        "new-source-observations.csv.gz": extra,
        "scope-summary.csv": summary,
        "funnel.csv": funnel,
    }
    for filename, frame in outputs.items():
        write_csv(frame, out / filename)
    (out / "tissue-taxonomy.json").write_text(json.dumps(normal_tissue_taxonomy(), indent=2) + "\n")
    render(out, current, intake, summary, funnel, policy, source_counts)
    implementation = [
        "scripts/cta_recalibration_report.py",
        "scripts/cta_permitted_report.py",
        "scripts/cta_scope_report.py",
        "scripts/cta_tissue_tier_report.py",
        "oncoref/cta_normal_screen.py",
        "oncoref/normal_tissues.py",
        "oncoref/cta_sources.py",
    ]
    inputs = [
        previous_path,
        observation_path,
        PREVIOUS / "run-manifest.json",
        OBSERVATIONS / "run-manifest.json",
    ]
    inputs += [
        ROOT / "oncoref/data" / name
        for name in [
            "cta-publication-membership.csv",
            "cta-publication-sources.csv",
            "cta-intake-identity-reviews.csv",
            "cta-intake-identity-snapshot.json",
            "cta-screening-review-overrides.csv",
            "cta-specificity-audit.csv",
            "cta-ihc-unreliable.csv",
            "normal-tissue-taxonomy.yaml",
            "normal-tissue-mappings.csv",
        ]
    ]
    manifest = {
        "report_schema": "oncoref.cta-recalibration-report.v1",
        "policy": policy,
        "routine_contexts": old_manifest["routine_contexts"],
        "new_measurement_sources": receipts,
        "checks": [
            "Pinned inputs verified",
            "Historical retained set is nested through each relaxed policy stage",
            "Fixed retained set is a subset of adaptive",
            "All intake identities retained including unassessed NCBI identities",
            "Final funnel reconciles exactly to current decisions",
        ],
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "implementation_sha256": {name: sha256(ROOT / name) for name in implementation},
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(funnel.to_string(index=False))
    print(intake.loc[intake.Symbol.isin(["BAGE", "GAGE2B", "SPAG9"])].to_string(index=False))


def markdown_table(frame):
    return "\n".join(
        "| " + " | ".join(str(v).replace("|", "/") for v in row) + " |"
        for row in [
            frame.columns,
            ["---"] * len(frame.columns),
            *frame.itertuples(index=False, name=None),
        ]
    )


def render(out, current, intake, summary, funnel, policy, source_counts):
    figure_style.apply()
    import matplotlib.pyplot as plt

    labels = [
        "v1.11.0",
        "Reproductive RNA >0",
        "Revised selectivity",
        "Protected baseline 3",
        "Current adaptive + reviews",
        "Current fixed + reviews",
    ]
    for filename, names, values, title in [
        (
            "policy-comparison",
            labels,
            summary.retained,
            "Reproductive admission, selectivity and review changes",
        ),
        (
            "intake-funnel",
            funnel.stage,
            funnel.remaining,
            "Published CTA intake through the current screen",
        ),
    ]:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.barh(list(names), list(values), color="#377eb8")
        ax.invert_yaxis()
        for i, value in enumerate(values):
            ax.text(value + max(values) * 0.01, i, str(value), va="center")
        ax.set_xlim(0, max(values) * 1.15)
        ax.set_xlabel(
            "Candidates" if filename == "policy-comparison" else "Distinct identities / candidates"
        )
        ax.set_title(title)
        fig.subplots_adjust(left=0.43, bottom=0.13, right=0.98, top=0.90)
        for suffix in ["png", "svg"]:
            figure_style.save(fig, out / f"{filename}.{suffix}", keep_titles=True)
        plt.close(fig)
    grid = pd.DataFrame(
        [
            {
                "IHC": confidence,
                **{
                    tier.title(): f"{policy['negative_ihc_rna_limits'][confidence][tier]:g} nTPM / {policy['negative_ihc_minimum_reproductive_selectivity'][confidence][tier]:.0%}"
                    for tier in TIERS
                },
            }
            for confidence in ["Missing", "Uncertain", "Approved", "Supported", "Enhanced"]
        ]
    )
    named = current.loc[
        current.Symbol.isin(
            [
                "CT83",
                "CTAG2",
                "SAGE1",
                "CTCFL",
                "GAGE1",
                "SSX4",
                "SSX5",
                "MAGEA6",
                "MAGEA10",
                "MAGEA12",
                "MAGEA11",
                "SPAG9",
            ]
        ),
        ["Symbol", "retained", "failure_reasons", "screen_review_applied"],
    ].sort_values("Symbol")
    short_intake = intake.loc[
        intake.Symbol.isin(["BAGE", "GAGE2B", "SPAG9"]),
        ["Symbol", "intake_identity", "screened", "screen_status"],
    ]
    text = f"""# CTA recalibration and complete intake — 21 October 2026

Policy `{policy["policy_id"]}`. {len(source_counts)} publication sources, {len(intake):,} distinct intake identities and {len(current):,} mapped protein-coding candidates. Human HPA v23 measurements remain pinned.

## Current rules

- Reproductive admission requires **valid RNA >0 nTPM OR positive reproductive IHC**. A measured zero does not qualify; unavailable observations remain missing.
- Disposable expression is exempt from exclusion and cannot supply the reproductive reference or admission evidence.
- `Selectivity = Rmax / (Rmax + tissue RNA)`, with Rmax the highest valid raw reproductive nTPM for the same gene/source/version/species.
- Each Manageable, Protected and Vital tissue must pass **both** its absolute ceiling and selectivity minimum. Any positive IHC, including Low, still excludes.
- The same complete routine local negative IHC selects both thresholds, using the weakest reliability. Missing, uncertain, incomplete or unusable IHC receives baseline limits. Zero/missing reproductive RNA leaves selectivity unresolved.

### Maximum nTPM / minimum selectivity

{markdown_table(grid)}

The user approved this target-informed research recalibration. The previous 99% minimum required a 99-fold reproductive/tissue ratio and could exclude 0.1 nTPM observations despite a 1 nTPM absolute ceiling. The new percentages were selected for the pairwise metric while retaining absolute and IHC gates; they are not validated safety limits. Protected's baseline changes generally from 2 to 3 nTPM, which admits CT83's salivary 2.6 nTPM. There are no gene-specific numerical exemptions. Named-target recovery is a calibration check, not independent validation.

## Policy comparison

{markdown_table(summary)}

The first row is the archived v1.11.0 cohort. Subsequent rows include the new SPAG9 candidate. `admission_only` changes reproductive admission; `selectivity_revised` then changes the percentage grid; `grid_revised` also changes Protected's baseline; `current_adaptive` additionally applies explicit holdback reviews. `current_fixed` uses baseline RNA/selectivity values without IHC relaxation.

![Policy comparison](policy-comparison.png)

## Current funnel

{markdown_table(funnel)}

This sequential ordering attributes each removal to its first failing stage. Failure reasons overlap; every reason remains in the candidate decisions. An intake identity lacking a canonical mapping is visible upstream and is unassessed, never assigned zero RNA.

![Source-to-screen funnel](intake-funnel.png)

## Requested targets and review provenance

{markdown_table(named)}

CTAG2's old categorical RNA holdback is explicitly superseded **for this screen**, conditional on all current expression and coverage gates. Its pinned heart 5.1 nTPM and Enhanced Not detected cardiomyocyte staining remain visible. The [heart evidence review](../ctag2-heart-evidence.md) retains unresolved cell origin, antigen-presentation questions and confounded clinical cardiac-event reports. MAGEA11 receives the same review of its historical trace-RNA holdback. Neither adjudication bypasses positive IHC or a numerical failure. The frozen core/extended panels retain their original reviews.

### BAGE, GAGE2B and SPAG9 intake

{markdown_table(short_intake)}

BAGE and GAGE2B have official human NCBI identities, confirmed in a checksum-pinned snapshot, but lack a usable mapping in the pinned Ensembl 115 canonical space. Their memberships, citations and identity reviews are preserved in the funnel. No paralog is substituted, and their HPA normal-tissue expression is unassessed here. Jamin's GAGE2B nomination uses shared microarray probes; it does not establish locus-specific expression. SPAG9 enters the coding screen from its primary ovarian-cancer study and receives the same HPA gates as every other gene.

New primary sources: [Boël et al. 1995, BAGE](https://pubmed.ncbi.nlm.nih.gov/7895173/) and [Garg et al. 2007, SPAG9](https://pubmed.ncbi.nlm.nih.gov/17332284/). Identity authorities: [NCBI BAGE](https://www.ncbi.nlm.nih.gov/gene/574), [NCBI GAGE2B](https://www.ncbi.nlm.nih.gov/gene/645037). Source assays and limitations are preserved separately from identity resolution.

## Coverage

Complete routine RNA is required. IHC-unavailable and partial-IHC candidates remain explicitly distinguished from routine-complete observations; none proves exhaustive protein or peptide-HLA absence. The routine panel contains 50 RNA and 80 IHC tissue/cell contexts. Missing developmental and anatomical compartments remain unassessed. Human assay limits do not transfer to dog or mouse. See [category definitions and mappings](../../normal-tissue-scopes.md).

## Downloads

- [All source identities and their stage/status](intake-decisions.csv)
- [BAGE/GAGE2B/SPAG9 intake audit](intake-target-audit.csv)
- [Per-publication funnel counts](source-counts.csv)
- [Exact source memberships and assay resolution](source-membership.csv.gz)
- [Publication registry](publication-sources.csv)
- [Candidate decisions, all scenarios](candidate-decisions.csv)
- [Current comparison](candidate-comparison.csv)
- [Membership changes](membership-changes.csv)
- [Named target decisions](notable-targets.csv)
- [Every RNA value, reproductive reference, selectivity and IHC support](rna-selectivity-audit.csv.gz)
- [Added SPAG9 observations](new-source-observations.csv.gz)
- [Policy and taxonomy](tissue-taxonomy.json)
- [Checksums and reproducibility](run-manifest.json)

Reproduce with `python scripts/cta_recalibration_report.py --out <empty-directory>` at this implementation revision. Existing candidate measurements are reused from the verified archive; the added SPAG9 observations come from the same checksum-pinned HPA v23 inputs.
"""
    (out / "index.md").write_text(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
