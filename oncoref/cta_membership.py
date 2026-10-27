"""Biological CTA/oncofetal membership, independent of normal-expression risk.

Established membership is source-backed, not a gene-name exception. Expression
nominees remain a separate evidence tier; none of these calls clear a peptide,
receptor, patient population, or normal-protein subtraction reference.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.cta-membership.2026-10-26.2"
DATA = Path(__file__).parent / "data"
CLASSES = {"Cancer-germline", "Cancer-placental", "Cancer-reproductive", "Oncofetal"}


@dataclass(frozen=True)
class Thresholds:
    """Positive expression support only; never normal-tissue safety limits."""

    reproductive_ntpm: float = 1.0
    reproductive_fold: float = 4.0
    tumor_prevalence: float = 0.05
    tumor_samples: int = 3
    tumor_ihc_cases: int = 2
    papers: int = 2

    def __post_init__(self):
        for key, value in asdict(self).items():
            if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"Positive finite threshold required: {key}")
        if self.tumor_prevalence > 1 or self.reproductive_fold < 1:
            raise ValueError("Invalid prevalence or enrichment threshold")
        for key in ("tumor_samples", "tumor_ihc_cases", "papers"):
            if int(getattr(self, key)) != getattr(self, key):
                raise ValueError(f"Integer threshold required: {key}")


def evidence() -> pd.DataFrame:
    """Reviewed identity, biological context, primary evidence and limitations."""
    return pd.read_csv(DATA / "cta-membership-evidence.csv", keep_default_na=False)


def controls() -> pd.DataFrame:
    """Benchmarks are evaluation labels, never classifier inputs."""
    return pd.read_csv(DATA / "cta-membership-controls.csv", keep_default_na=False)


def classify(
    facts: pd.DataFrame, reviewed: pd.DataFrame, *, thresholds: Thresholds | None = None
) -> pd.DataFrame:
    """Apply three common evidence routes to every canonical human gene.

    1. Reviewed primary evidence establishes reproductive/developmental and
       tumor association. Modern normal-expression concerns do not revoke it.
    2. At least two distinct class-specific cancer/reproductive nomination
       papers support a publication-consensus call, even if the current atlas
       misses rare tumors or reproductive cell states reported in those papers.
    3. Measured reproductive enrichment plus published or measured tumor
       expression supports a candidate, including previously unlisted genes.

    pTPM tumor values are never compared with normal nTPM. Tumor observations
    are aggregated upstream within each cohort, not across cohorts. Their
    recorded support thresholds must match this call's parameters.
    """
    t = thresholds or Thresholds()
    required = {
        "Gene ID",
        "Gene",
        "Taxonomy ID",
        "Reference",
        "Core RNA (nTPM)",
        "Core comparator maximum (nTPM)",
        "Core tissues measured",
        "Core comparator tissues measured",
        "Cancer-germline papers",
        "Cancer-placental papers",
        "Publication classes",
        "Tumor RNA positive",
        "Tumor IHC positive",
        "Testis RNA (nTPM)",
        "Ovary RNA (nTPM)",
        "Placenta RNA (nTPM)",
        "Tumor minimum samples",
        "Tumor minimum prevalence",
        "Tumor IHC minimum cases",
        "Other reproductive RNA (nTPM)",
        "Other reproductive tissues measured",
        "Nonreproductive maximum (nTPM)",
        "Nonreproductive tissues measured",
    }
    if not required <= set(facts):
        raise ValueError(f"Missing membership facts: {sorted(required - set(facts))}")
    if facts.empty or facts["Gene ID"].duplicated().any() or facts["Gene ID"].isna().any():
        raise ValueError("Canonical gene IDs must be unique and present")
    if not facts["Taxonomy ID"].eq(9606).all() or facts.Reference.nunique() != 1:
        raise ValueError("A single explicit human reference is required")
    required_evidence = {
        "Gene ID",
        "Gene",
        "Taxonomy ID",
        "Reference",
        "Class",
        "Identity scope",
        "Biological evidence",
        "Tumor evidence",
        "Sources",
        "Review depth",
        "Limitations",
        "Review date",
        "Admission role",
        "Admission rationale",
    }
    if not required_evidence <= set(reviewed):
        raise ValueError("Incomplete reviewed membership evidence")
    if reviewed.duplicated(["Gene ID", "Class"]).any():
        raise ValueError("One consolidated review per gene/class is required")
    if not set(reviewed["Gene ID"]) <= set(facts["Gene ID"]):
        raise ValueError("Reviewed identities must exist in the canonical universe")
    if not reviewed.Class.isin(CLASSES).all():
        raise ValueError("Unknown biological class")
    if not reviewed["Admission role"].isin(["Membership support", "Context only"]).all():
        raise ValueError("Explicit reviewed-evidence admission role required")
    if not reviewed["Identity scope"].isin(["Locus", "Shared antigen; locus unresolved"]).all():
        raise ValueError("Unresolved families cannot be assigned to arbitrary loci")
    if (
        not reviewed["Taxonomy ID"].eq(9606).all()
        or not reviewed.Reference.eq(facts.Reference.iloc[0]).all()
    ):
        raise ValueError("Evidence reference/species mismatch")
    identity = facts.set_index("Gene ID").Gene
    if not reviewed.Gene.eq(reviewed["Gene ID"].map(identity)).all():
        raise ValueError("Reviewed gene symbol/reference identity mismatch")
    for col in required_evidence - {"Taxonomy ID"}:
        if reviewed[col].isna().any() or reviewed[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Explicit evidence required: {col}")
    if not reviewed.Sources.str.startswith("https://").all():
        raise ValueError("Primary source URLs required")
    pd.to_datetime(reviewed["Review date"], format="%Y-%m-%d", errors="raise")
    d = facts.copy()
    for col in [
        "Core RNA (nTPM)",
        "Testis RNA (nTPM)",
        "Ovary RNA (nTPM)",
        "Placenta RNA (nTPM)",
        "Core comparator maximum (nTPM)",
        "Other reproductive RNA (nTPM)",
        "Nonreproductive maximum (nTPM)",
    ]:
        values = pd.to_numeric(d[col], errors="raise")
        if np.isinf(values).any() or values.lt(0).any():
            raise ValueError("Expression must be nonnegative finite values or explicit missing")
    for col in ["Cancer-germline papers", "Cancer-placental papers"]:
        values = pd.to_numeric(d[col], errors="raise")
        if not np.isfinite(values).all() or values.lt(0).any() or values.mod(1).ne(0).any():
            raise ValueError("Publication counts must be nonnegative integers")
    for col in ["Tumor RNA positive", "Tumor IHC positive"]:
        if not d[col].map(lambda value: isinstance(value, (bool, np.bool_))).all():
            raise ValueError(f"Explicit booleans required: {col}")
    for col, value in {
        "Tumor minimum samples": t.tumor_samples,
        "Tumor minimum prevalence": t.tumor_prevalence,
        "Tumor IHC minimum cases": t.tumor_ihc_cases,
    }.items():
        if not d[col].eq(value).all():
            raise ValueError("Reaggregate tumor support at the requested thresholds")
    measured = d["Core RNA (nTPM)"].ge(t.reproductive_ntpm)
    tumor = d["Tumor RNA positive"] | d["Tumor IHC positive"]
    consensus_classes = {
        "Cancer-germline": d["Cancer-germline papers"].ge(t.papers),
        "Cancer-placental": d["Cancer-placental papers"].ge(t.papers),
    }
    consensus = consensus_classes["Cancer-germline"] | consensus_classes["Cancer-placental"]
    core_enriched = (
        measured
        & (tumor | d[["Cancer-germline papers", "Cancer-placental papers"]].max(axis=1).ge(1))
        & d["Core tissues measured"].eq(3)
        & d["Core comparator tissues measured"].eq(47)
        & d["Core RNA (nTPM)"].ge(t.reproductive_fold * d["Core comparator maximum (nTPM)"])
    )
    # Extended reproductive organs require an explicit cancer/reproductive
    # publication. Normal breast/prostate lineage expression alone must not turn
    # a generic tumor-associated differentiation antigen into a CTA nominee.
    other_enriched = (
        d[["Cancer-germline papers", "Cancer-placental papers"]].max(axis=1).ge(1)
        & d["Other reproductive RNA (nTPM)"].ge(t.reproductive_ntpm)
        & d["Other reproductive tissues measured"].eq(8)
        & d["Nonreproductive tissues measured"].eq(39)
        & d["Other reproductive RNA (nTPM)"].ge(
            t.reproductive_fold * d["Nonreproductive maximum (nTPM)"]
        )
    )
    enriched = core_enriched | other_enriched
    admitted = reviewed.loc[reviewed["Admission role"].eq("Membership support")]
    admission_groups = admitted.groupby("Gene ID", sort=True)
    grouped = reviewed.groupby("Gene ID", sort=True)
    d["Reviewed classes"] = (
        d["Gene ID"].map(admission_groups.Class.agg(lambda x: "|".join(sorted(set(x))))).fillna("")
    )
    d["Reviewed admission role"] = (
        d["Gene ID"]
        .map(grouped["Admission role"].agg(lambda x: "|".join(sorted(set(x)))))
        .fillna("No gene-specific review")
    )
    d["Reviewed admission rationale"] = (
        d["Gene ID"].map(grouped["Admission rationale"].agg(lambda x: " | ".join(x))).fillna("")
    )
    d["Membership sources"] = (
        d["Gene ID"].map(grouped.Sources.agg(lambda x: "|".join(sorted(set(x))))).fillna("")
    )
    d["Membership limitations"] = (
        d["Gene ID"].map(grouped.Limitations.agg(lambda x: " | ".join(x))).fillna("")
    )
    d["Identity scope"] = (
        d["Gene ID"]
        .map(grouped["Identity scope"].first())
        .fillna("RNA/protein mapping; paralog specificity unverified")
    )
    established = d["Reviewed classes"].ne("")
    d["Member"] = established | consensus | enriched
    d["Established route"] = established
    d["Consensus route"] = consensus
    d["Expression route"] = enriched
    d["Core expression route"] = core_enriched
    d["Other reproductive route"] = other_enriched
    d["Evidence tier"] = np.select(
        [established, consensus, enriched],
        ["Established", "Publication consensus", "Expression-supported"],
        default="Insufficient evidence",
    )
    classes = []
    for index, row in d.iterrows():
        values = set(row["Reviewed classes"].split("|")) - {""}
        if consensus.loc[index]:
            for category in ["Cancer-germline", "Cancer-placental"]:
                if consensus_classes[category].loc[index]:
                    values.add(category)
        if core_enriched.loc[index]:
            for tissue, category in [
                ("Testis", "Cancer-germline"),
                ("Ovary", "Cancer-reproductive"),
                ("Placenta", "Cancer-placental"),
            ]:
                value = row[tissue + " RNA (nTPM)"]
                if (
                    value >= t.reproductive_ntpm
                    and value >= t.reproductive_fold * row["Core comparator maximum (nTPM)"]
                ):
                    values.add(category)
        if other_enriched.loc[index]:
            values.add("Cancer-reproductive")
        classes.append("|".join(sorted(values - {""})))
    d["Classes"] = classes
    d["Membership reason"] = np.select(
        [established, consensus, enriched],
        [
            "Reviewed reproductive/developmental and tumor association",
            "At least two distinct class-specific cancer/reproductive publications",
            "Measured reproductive enrichment plus published or measured tumor support",
        ],
        default="No qualifying reviewed, publication-consensus or expression route",
    )
    d["Membership policy"] = POLICY_ID
    return d


def results() -> pd.DataFrame:
    """Genome-wide biological membership and retained normal-expression flags."""
    d = pd.read_csv(DATA / "cta-membership-gene-sets.csv.gz", keep_default_na=False)
    for col in [
        "Gene",
        "Reproductive RNA (nTPM)",
        "Core RNA (nTPM)",
        "Testis RNA (nTPM)",
        "Ovary RNA (nTPM)",
        "Placenta RNA (nTPM)",
        "Core comparator maximum (nTPM)",
        "Core / comparator maximum",
        "Other reproductive RNA (nTPM)",
        "Nonreproductive maximum (nTPM)",
    ]:
        d[col] = d[col].replace("", np.nan)
        if col != "Gene":
            d[col] = pd.to_numeric(d[col])
    return d


def gene_ids(*, established_only: bool = False, category: str | None = None) -> set[str]:
    """Holistic union by default; optionally restrict evidence tier or class."""
    if category is not None and category not in CLASSES:
        raise ValueError(f"Unknown biological class: {category}")
    d = results()
    mask = d.Member
    if established_only:
        mask &= d["Evidence tier"].eq("Established")
    if category is not None:
        mask &= d.Classes.str.split("|").map(lambda values: category in values)
    return set(d.loc[mask, "Gene ID"])
