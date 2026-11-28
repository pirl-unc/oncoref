"""Tests for the OncoRef antigen preprint and surface/ectodomain modality audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from oncoref import target_selection

ROOT = Path(__file__).resolve().parents[1]
PREPRINT_DIR = ROOT / "docs/preprint"


def test_preprint_artifacts_exist_and_match_policy():
    meta = json.loads((PREPRINT_DIR / "data/uniprot-candidate-topology.json").read_text())
    gz_bytes = (PREPRINT_DIR / "data/uniprot-candidate-topology.json.gz").read_bytes()
    assert hashlib.sha256(gz_bytes).hexdigest() == meta["sha256"]
    assert meta["policy_id"] == target_selection.POLICY_ID

    md_path = PREPRINT_DIR / "index.md"
    tex_path = PREPRINT_DIR / "oncoref-antigen-preprint.tex"
    pdf_path = PREPRINT_DIR / "oncoref-antigen-preprint.pdf"
    assert md_path.is_file() and md_path.stat().st_size > 20_000
    assert tex_path.is_file() and tex_path.stat().st_size > 15_000
    assert pdf_path.is_file() and pdf_path.read_bytes().startswith(b"%PDF-")
    assert pdf_path.stat().st_size > 100_000

    for stem in (
        "fig1-conceptual-architecture-and-attrition",
        "fig2-multispecies-parameterization-and-orthology",
        "fig3-human-deep-dive-and-ca-only-biology",
        "fig4-clinical-tcr-gold-standard",
        "fig5-surface-ectodomain-modality-routing",
    ):
        assert (PREPRINT_DIR / "figures" / f"{stem}.png").is_file()
        assert (PREPRINT_DIR / "figures" / f"{stem}.pdf").is_file()


def test_surface_topology_audit_covers_all_human_candidates_and_benchmarks():
    topo = pd.read_csv(PREPRINT_DIR / "tables/candidate-surface-topology.csv")
    ts = target_selection.results()
    cand_ids = set(ts.loc[ts["is_target_candidate"].astype(bool), "entity_id"])
    assert cand_ids <= set(topo["entity_id"])

    ca_topo = topo[topo["category"].isin(["CTA", "Cancer Antigen only"])]
    assert len(ca_topo) == 375
    surf_ca = ca_topo[ca_topo["cell_surface_exposed"].astype(bool)]
    assert len(surf_ca) == 24

    by_lbl = topo.set_index("protein_group_label")
    # Intracellular TCR-T anchors have 0 surface ECD exposure
    for intracell in ("PRAME", "NY-ESO-1", "MAGEA4", "CTAG2", "DUX4", "AKAP4"):
        assert not bool(by_lbl.loc[intracell, "cell_surface_exposed"])

    # Actionable >=50 aa surface ectodomain CTAs & CA-only targets
    for surf_50 in (
        "BTN1A1",
        "ERVV-1",
        "ERVV-2",
        "ADAM2",
        "ADAM18",
        "FMR1NB",
        "CT83",
        "ZP4",
        "IZUMO1R",
        "LYPD4",
        "CIROP",
        "CRIPTO3",
        "PKDREJ",
    ):
        assert bool(by_lbl.loc[surf_50, "cell_surface_exposed"])
        assert int(by_lbl.loc[surf_50, "ecd_max_continuous_aa"]) >= 50

    # All 24 clinical ADC/CAR-T/RLT surface benchmarks are Unsafe under TCR-grade screening
    bench = topo[topo["clinical_alias"].fillna("").ne("")]
    assert len(bench) == 24
    assert bench["category"].eq("Unsafe").all()
