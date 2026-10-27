#!/usr/bin/env python3
"""Export the measured contexts behind PRAME, CTAG2 and TBXT decisions.

This audit does not tune cutoffs or confer a clinical safety classification.
Clinical observations are separately curated with denominators and source URLs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_quantitative_threshold_report import positive_profile_controls  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import cell_rna_policy as policy  # noqa: E402
from oncoref import hpa, normal_cell_data, normal_cell_screen, reference_data  # noqa: E402

TARGETS = {
    "ENSG00000185686": "PRAME",
    "ENSG00000126890": "CTAG2",
    "ENSG00000164458": "TBXT",
}
BASELINE = ROOT / "docs/audits/cta-expanded-funnel-20261027"
CURRENT = ROOT / "docs/audits/cta-harmonized-20261027"
OUTPUT = ROOT / "docs/audits/cta-target-evidence-20261027"


def main() -> None:
    """Replay the existing cell policy and publish complete target context rows."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(BASELINE / "gene-classifications.csv.gz", low_memory=False)
    conflicts = set(facts.loc[facts["Antibody concern"], "gene_id"])
    ihc = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    controls = positive_profile_controls(ihc)
    raw = normal_cell_data.table("hpa_single_cell_tissue", genes=list(TARGETS))
    write_csv(raw, OUTPUT / "hpa-cluster-observations.csv")
    write_csv(ihc[ihc.Gene.isin(TARGETS)], OUTPUT / "hpa-ihc-observations.csv")
    d = raw.rename(
        columns={"Gene": "gene_id", "Tissue": "tissue", "Cell type": "cell_type", "nTPM": "value"}
    ).copy()
    d["source"], d["version"], d["taxonomy_id"] = "hpa_single_cell_tissue", "v23", 9606
    d["state"], d["stage"], d["unit"] = "unspecified", "adult", "nTPM"
    d["value"] = d.value.where(d["Read count"].ge(20))
    d = normal_cell_screen.contexts(d)
    d = policy.reproductive_immune_contexts(d)
    d = policy.protected_immune_contexts(d)
    d = policy.match_negative_ihc(
        d, ihc[ihc.Gene.isin(TARGETS)], positive_controls=controls, conflicting_genes=conflicts
    )
    d = policy.assess(
        d,
        policy.linked_thresholds(mode="two_axis_doubling"),
        ihc_factors=dict.fromkeys(["Missing/uncertain", *policy.GRADES], 1),
        inclusive_ceiling=True,
    )
    write_csv(d, OUTPUT / "assessed-cluster-contexts.csv")
    bulk = pd.read_csv(BASELINE / "normal/bulk-contexts.csv.gz", low_memory=False)
    write_csv(bulk[bulk.Ensembl_Gene_ID.isin(TARGETS)], OUTPUT / "assessed-bulk-contexts.csv")
    peaks = pd.read_csv(BASELINE / "normal/normal-peak-evidence.csv.gz")
    write_csv(peaks[peaks.gene_id.isin(TARGETS)], OUTPUT / "all-assay-peaks.csv")
    decisions = pd.read_csv(CURRENT / "named-antigen-decisions.csv")
    write_csv(
        decisions[decisions.representative_gene_id.isin(TARGETS)], OUTPUT / "current-decisions.csv"
    )
    inputs = [
        Path(__file__),
        ROOT / "oncoref/cell_rna_policy.py",
        ROOT / "oncoref/normal_cell_screen.py",
        ROOT / "scripts/cta_quantitative_threshold_report.py",
        BASELINE / "normal/manifest.json",
        BASELINE / "gene-classifications.csv.gz",
        CURRENT / "manifest.json",
    ]
    sources = []
    for source in ["hpa_single_cell_tissue", "hpa_normal_tissue"]:
        receipt = reference_data.provenance(source, verify_content=True)
        sources.append({k: receipt[k] for k in ["name", "version", "url", "sha256"]})
    manifest = {
        "review_date": "2026-10-27",
        "clinical_evidence_scope": "targeted review; unassessed targets are not negative evidence",
        "classification_override": False,
        "inputs": input_hashes(inputs),
        "sources": sources,
        "outputs": {
            p.name: sha256(p) for p in sorted([*OUTPUT.glob("*.csv"), OUTPUT / "index.md"])
        },
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        decisions.loc[
            decisions.representative_gene_id.isin(TARGETS), ["protein_group_label", "category"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
