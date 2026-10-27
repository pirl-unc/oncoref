#!/usr/bin/env python3
"""Reproduce descriptive breadth from pinned HPA v23, without changing membership."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import cta_policy, hpa, reference_data  # noqa: E402
from oncoref.normal_expression_breadth import (  # noqa: E402
    POLICY_ID,
    anatomical_panel,
    burden,
    donor_summaries,
)


def digest(path):
    """Return the content SHA-256 used to pin an audit input or output."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(frame: pd.DataFrame, path):
    """Write a data table with deterministic compression when requested."""
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.name.endswith(".gz") else None,
    )


def build(out: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory; preserve historical reports")
    out.mkdir(parents=True, exist_ok=True)
    current = ROOT / "docs/audits/cta-policy-20261025"
    screening = ROOT / "docs/audits/cta-screen-20261023"
    audit_path = screening / "tissue-evidence.csv.gz"
    screen_manifest = json.loads((screening / "run-manifest.json").read_text())
    if digest(audit_path) != screen_manifest["outputs"][audit_path.name]:
        raise ValueError("Pinned tissue evidence checksum mismatch")
    candidates = pd.read_csv(current / "candidate-results.csv.gz")
    if len(candidates) != 2540 or not candidates["Gene ID"].is_unique:
        raise ValueError("Expected 2,540 unique mapped nominations")
    candidate_ids = set(candidates["Gene ID"])
    # Prespecified ordinary-normal controls, independent of score values.
    controls = {"ENSG00000075624": "ACTB", "ENSG00000111640": "GAPDH", "ENSG00000196262": "PPIA"}
    wanted = candidate_ids | controls.keys()
    receipt = reference_data.provenance("hpa_rna_consensus", "v23")
    if digest(receipt["path"]) != receipt["sha256"]:
        raise ValueError("HPA cache checksum mismatch")
    native = hpa._read_hpa("hpa_rna_consensus", "v23")
    panel = anatomical_panel()
    native = native.loc[native.Gene.isin(wanted) & native.Tissue.isin(panel.compartment)]
    source = native.rename(columns={"Gene": "gene_id", "Tissue": "compartment", "nTPM": "value"})
    baseline = pd.MultiIndex.from_product(
        [sorted(wanted), panel.compartment], names=["gene_id", "compartment"]
    ).to_frame(index=False)
    baseline = baseline.merge(
        source[["gene_id", "compartment", "value"]],
        how="left",
        on=["gene_id", "compartment"],
        validate="one_to_one",
    )
    categories = panel.set_index("compartment").category
    manifest = cta_policy.manifest()
    ceilings = {k: v[0] for k, v in manifest["rna_ceiling_ntpm"].items()}
    baseline["ceiling"] = baseline.compartment.map(categories).map(ceilings)
    baseline["source"], baseline["version"], baseline["taxonomy_id"] = (
        "hpa_rna_consensus",
        "v23",
        9606,
    )
    baseline["unit"], baseline["cohort"], baseline["summary"] = "nTPM", "HPA consensus", "consensus"
    audit = pd.read_csv(audit_path)
    audit = audit.loc[audit.tissue.isin(panel.compartment)]
    if not audit.rules_policy_id.eq(manifest["measurement_policy"]).all():
        raise ValueError("Tissue evidence does not match the active measurement policy")
    effective = baseline.loc[baseline.gene_id.isin(candidate_ids)].drop(
        columns=["value", "ceiling"]
    )
    effective = effective.merge(
        audit[["Ensembl_Gene_ID", "tissue", "rna_ntpm", "rna_ceiling_ntpm"]].rename(
            columns={
                "Ensembl_Gene_ID": "gene_id",
                "tissue": "compartment",
                "rna_ntpm": "value",
                "rna_ceiling_ntpm": "ceiling",
            }
        ),
        on=["gene_id", "compartment"],
        how="left",
        validate="one_to_one",
    )
    effective["ceiling"] = effective.ceiling.fillna(
        effective.compartment.map(categories).map(ceilings)
    )
    results = []
    for scenario, profiles in [("fixed_baseline", baseline), ("applicable_ihc", effective)]:
        for weight in [0.5, 1.0, 2.0]:
            parts, scores = burden(profiles, anatomical_panel(brain_weight=weight))
            scores["scenario"], scores["brain_weight"] = scenario, weight
            results.append(scores)
            if weight == 1:
                write(parts, out / f"{scenario}-profiles.csv.gz")
    scores = pd.concat(results, ignore_index=True)
    names = candidates.set_index("Gene ID").Gene.to_dict() | controls
    scores["Gene"] = scores.gene_id.map(names)
    scores["Category"] = scores.gene_id.map(candidates.set_index("Gene ID").Category).fillna(
        "Normal control"
    )
    write(scores, out / "breadth.csv.gz")
    write(panel, out / "anatomical-panel.csv")
    donors = pd.read_csv(ROOT / "oncoref/data/cta-donor-rna.csv").rename(
        columns={
            "Ensembl_Gene_ID": "gene_id",
            "tissue": "compartment",
            "source_version": "version",
            "rna_ntpm": "value",
        }
    )
    if not donors.unit.eq("HPA-v23-nTPM").all():
        raise ValueError("Unexpected donor units")
    donors["source"], donors["unit"] = "hpa_tissue_donors", "nTPM"
    donor_summary = donor_summaries(donors)
    write(donor_summary, out / "donor-summaries.csv")
    donor_summary["ceiling"] = donor_summary.compartment.map(categories).map(ceilings)
    donor_profiles, donor_scores = burden(donor_summary, panel)
    donor_profiles = donor_profiles.merge(
        donor_summary.drop(columns=["value", "ceiling"]),
        on=[
            "gene_id",
            "source",
            "version",
            "taxonomy_id",
            "unit",
            "cohort",
            "summary",
            "compartment",
        ],
        how="left",
        validate="one_to_one",
    )
    write(donor_profiles, out / "donor-profiles.csv.gz")
    write(donor_scores, out / "donor-breadth.csv")
    spotlight = scores.loc[
        scores.Gene.isin(
            [
                "PRAME",
                "MAGEA4",
                "CTAG1A",
                "CTAG1B",
                "ZNF730",
                "MAGEA12",
                "MAGEA3",
                "CTAG2",
                "ACTB",
                "GAPDH",
                "PPIA",
            ]
        )
        & scores.brain_weight.eq(1)
    ]
    write(spotlight, out / "benchmarks.csv")
    (out / "index.md").write_text(f"""# Continuous normal-expression breadth

Descriptive policy `{POLICY_ID}`; active eligibility remains `{manifest["policy"]}`.
All 2,540 mapped nominees and three prespecified ordinary-normal controls are
represented. Scores do not admit or exclude any gene. Clinical benchmark genes
were chosen before scoring; this is a descriptive comparison, not validation.

For each category, `B = sum(weight * RNA / (RNA + ceiling))`. Vital, Protected
and Manageable are separate. Reproductive abundance never enters this formula.
The [complete profiles](applicable_ihc-profiles.csv.gz) use the current archived
per-gene IHC-adjusted RNA ceilings. A [fixed-baseline comparison](fixed_baseline-profiles.csv.gz)
uses 2/6/10 nTPM for every gene, including ACTB, GAPDH and PPIA controls.
[All scores and weighting sensitivity](breadth.csv.gz) · [Benchmarks](benchmarks.csv).

## Fixed anatomy and missing coverage

The [panel](anatomical-panel.csv) is specified independently of each gene's
observed tissues. Each organ has equal weight within its category; eight brain
regions divide one brain weight. Retina and unmeasured non-retinal eye divide
one eye weight. Broad small intestine is used once, without adding overlapping
duodenum. Bronchus and embedded vascular/immune populations are not counted as
additional independent whole organs. Cellular blind spots still require the
separate cell-expression audit. This panel is not whole-body cell coverage.

Blood, peripheral nerve, inner ear, non-retinal eye and several protected
compartments are explicit missing coverage where this bulk source has no sample.
Observed weighted contribution is the lower bound; adding missing weight gives
the upper bound. These are coverage bounds, **not confidence intervals**. An
all-missing gene gets [0,1], never a negative finding. Six nominees absent from
the archived tissue screen retain their missing profiles in that scenario.
Sensitivity halves/doubles the prespecified brain organ weight, conserving
within-brain weights and renormalizing the fixed panel, never just observed rows.

## Donors, sources and units

HPA consensus is a published cross-source summary, not a donor median or an
independent cohort. Individual-donor counts and sampling intervals cannot be
inferred from it. [Available reviewed donor summaries](donor-summaries.csv)
currently cover CTAG2 heart only: four independent donors permit a median;
the prespecified ten-donor minimum leaves q90 unavailable. Separate
[donor burdens](donor-breadth.csv) and [profiles](donor-profiles.csv.gz) use the
fixed baseline ceilings, because these RNA donors are not paired to IHC. The
unavailable q90 retains entirely missing coverage, never a zero measurement.
A marginal 95%
percentile-bootstrap interval uses 2,000 donor resamples with seed 0; it is
reported separately from anatomical coverage and is unstable with four donors.
This tiny selected
sample is not an estimate of population risk, nor matched RNA/IHC specimens.
No whole-profile sampling interval is claimed. Native TPM/CP10K, cell atlases,
other species and HPA versions are excluded pending their own calibration.

## Interpretation and reproduction

No classification cutoff has been fitted. Accepted focal protein findings,
hard RNA exclusions and current holds remain in the eligibility policy; a low
burden cannot erase them. This score is a review-priority descriptor, not injury
probability or a body-wide expression fraction.

Run `.venv/bin/python scripts/build_normal_expression_breadth.py --out <empty-directory>`.
The manifest pins inputs, implementation, source receipt and all outputs.
""")
    inputs = [
        audit_path,
        current / "candidate-results.csv.gz",
        ROOT / "oncoref/data/cta-policy.json",
        ROOT / "oncoref/data/cta-donor-rna.csv",
        Path(__file__),
        ROOT / "oncoref/normal_expression_breadth.py",
    ]
    report = {
        "policy": POLICY_ID,
        "eligibility_changes": 0,
        "source": {k: receipt[k] for k in ["url", "sha256", "bytes", "version"]},
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "outputs": {p.name: digest(p) for p in sorted(out.iterdir())},
    }
    (out / "run-manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {len(scores)} category/scenario records for {len(wanted)} genes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
