"""The geometric comparison must preserve coverage, controls and count accounting."""

from pathlib import Path

import pandas as pd


def test_geometric_audit_replays_prior_membership_and_accounts_for_every_loss():
    root = Path(__file__).resolve().parents[1] / "docs/audits/cta-geometric-thresholds-20261027"
    genes = pd.read_csv(root / "gene-comparison.csv.gz")
    assert len(genes) == 20097 and genes.gene_id.is_unique
    assert genes.legacy_maximum10_selected.sum() == 337
    assert genes.ihc_doubling_maximum10_selected.sum() == 337
    assert genes.two_axis_doubling_maximum10_selected.sum() == 319
    assert genes.membership_change.eq("Lost").sum() == 18
    assert not genes.membership_change.eq("Gained").any()
    named = genes.set_index("symbol")
    for name in [
        "SMC1B",
        "PAGE1",
        "MAGEA4",
        "CTAG1A",
        "CTAG1B",
        "XAGE1A",
        "XAGE1B",
        "SSX2",
        "CT83",
        "TBXT",
    ]:
        assert named.loc[name, "two_axis_doubling_maximum10_selected"]
    for name in ["PRAME", "CTAG2", "PSMA8"]:
        assert not named.loc[name, "two_axis_doubling_maximum10_selected"]
    steps = pd.read_csv(root / "filter-cascade.csv")
    for _, d in steps.groupby("mode", sort=False):
        assert d.before.sub(d.removed).equals(d.remaining)
        assert d.before.iloc[1:].tolist() == d.remaining.iloc[:-1].tolist()
        assert d.removed.sum() + d.remaining.iloc[-1] == len(genes)
    reasons = pd.read_csv(root / "membership-change-reasons.csv")
    assert set(reasons.gene_id) == set(genes.loc[genes.membership_change.eq("Lost"), "gene_id"])
    assert reasons.previous_ceiling.eq(6).all() and reasons.new_ceiling.eq(4).all()
