# Licensed under the Apache License, Version 2.0.
"""Ordered CTA research rules and an excluded, reconsiderable Variable category.

Donor measurements are never pooled across cohorts/releases or imputed from
missing data. Variable is an exclusion from the current screen, not a claim
that an RNA-zero donor is eligible for treatment.
"""

from __future__ import annotations

import json
from functools import lru_cache

import numpy as np
import pandas as pd

from .cta_continuous_screen import ContinuousPolicy, continuous_tissue_screen
from .load_dataset import get_data
from .normal_tissues import SCREENED_NORMAL_TISSUE_TIERS, resolve_normal_tissue

POLICY_ID = "oncoref.cta-screen.2026-10-23.1"
BASELINES = {"vital": (2.0, 0.80), "protected": (6.0, 0.75), "manageable": (10.0, 0.70)}
RELIABLE = {"Approved", "Supported", "Enhanced"}
POSITIVE = {"Low", "Medium", "High", "Ascending", "Descending"}
NUMERICAL = {"rna_above_tissue_ceiling", "reproductive_selectivity_below_threshold"}
# Independent HPA v23 routine RNA inventory from the published 21 October
# expression-category manifest, restricted to screened normal compartments.
NORMAL_PANEL = frozenset(
    {
        "adipose tissue",
        "adrenal gland",
        "amygdala",
        "basal ganglia",
        "bone marrow",
        "cerebellum",
        "cerebral cortex",
        "choroid plexus",
        "colon",
        "duodenum",
        "esophagus",
        "heart muscle",
        "hippocampal formation",
        "hypothalamus",
        "kidney",
        "liver",
        "lung",
        "lymph node",
        "midbrain",
        "pancreas",
        "parathyroid gland",
        "pituitary gland",
        "rectum",
        "retina",
        "salivary gland",
        "skeletal muscle",
        "skin",
        "small intestine",
        "smooth muscle",
        "spinal cord",
        "spleen",
        "stomach",
        "thymus",
        "thyroid gland",
        "tongue",
        "urinary bladder",
    }
)


def reviewed_donor_rna():
    """Reviewed individual RNA measurements; presently CTAG2 heart only.

    No row means unreviewed, never homogeneous/negative. Units are pinned HPA
    v23 nTPM; raw GTEx TPM is deliberately not compared with these ceilings.
    """
    return get_data("cta-donor-rna").copy()


def results():
    """Published ordered-rule results; Retained is not clinical clearance."""
    return get_data("cta-screen-results").copy()


def retained_gene_names():
    """Known symbols in the selection; use retained_gene_ids for every locus."""
    frame = results()
    return set(frame.loc[frame.Category.eq("Retained"), "Gene"].dropna()) - {""}


def retained_gene_ids():
    """All selected Ensembl gene IDs, including loci without a current symbol."""
    frame = results()
    return set(frame.loc[frame.Category.eq("Retained"), "Gene ID"])


def variable_gene_names():
    """Excluded candidates with the reviewed donor-variable vital RNA pattern."""
    frame = results()
    return set(frame.loc[frame.Category.eq("Variable"), "Gene"].dropna()) - {""}


def variable_gene_ids():
    """All Variable Ensembl gene IDs; these are excluded from selection."""
    frame = results()
    return set(frame.loc[frame.Category.eq("Variable"), "Gene ID"])


@lru_cache(maxsize=512)
def _annotation(tissue, cell=None):
    return resolve_normal_tissue(
        tissue, cell_type=cell, species="human", source_id="hpa_normal_tissue", source_version="v23"
    )


def reproductive_controls(audit):
    """Reliable positive reproductive compartments in the curated gene profile.

    The consensus export does not identify antibody/run/donor pairs. This is
    profile-level support, not a verified matched-specimen positive control or
    an assay-sensitivity calibration. Cell anatomy, not just organ, must map.
    """
    found = {}
    for row in audit.itertuples():
        for cell in json.loads(row.supporting_ihc):
            if cell["ihc_level"] not in POSITIVE or cell["ihc_reliability"] not in RELIABLE:
                continue
            mapped = _annotation(cell["tissue"], cell["cell_type"])
            if "reproductive" in mapped.scopes:
                found.setdefault(row.Ensembl_Gene_ID, set()).add(
                    f"{cell['tissue']}: {cell['cell_type']} ({cell['ihc_level']})"
                )
    return {gid: "; ".join(sorted(values)) for gid, values in found.items()}


def screen_tissues(audit):
    """Apply the tighter selectivity ladder and reproductive-control requirement.

    All selectivity minima increase by five percentage points. Existing absolute
    ceilings remain. Negative-IHC relaxation requires a reliable reproductive
    positive in the same curated HPA gene profile; otherwise baseline limits
    apply. The continuous allowance remains 0.1 nTPM and R must exceed N.
    """
    if "rules_policy_id" in audit:
        raise ValueError("Supply the original HPA audit, not already revised tissue rules")
    controls = reproductive_controls(audit)
    revised = audit.copy()
    for key in ["rna_ceiling_ntpm", "minimum_reproductive_selectivity"]:
        revised["previous_" + key] = revised[key]
        revised[key] = pd.to_numeric(revised[key], errors="coerce")
    revised["reproductive_ihc_control"] = revised.Ensembl_Gene_ID.map(controls).fillna("")
    tested = revised.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)
    revised.loc[tested, "minimum_reproductive_selectivity"] = (
        revised.loc[tested, "minimum_reproductive_selectivity"] + 0.05
    ).round(6)
    usable = revised.reproductive_ihc_control.ne("") & revised.ihc_evidence_status.isin(
        ["negative_approved", "negative_supported", "negative_enhanced"]
    )
    for tier, (ceiling, selectivity) in BASELINES.items():
        baseline = revised.screening_tier.eq(tier) & ~usable
        revised.loc[baseline, "rna_ceiling_ntpm"] = ceiling
        revised.loc[baseline, "minimum_reproductive_selectivity"] = selectivity
    revised["negative_ihc_relaxation_used"] = tested & usable
    result = continuous_tissue_screen(revised, policy=ContinuousPolicy())
    result["rules_policy_id"] = POLICY_ID
    return result


def screen_donors(audit, donors):
    """Compare every reviewed donor with exactly the main screen's tissue rules.

    Audit must be the original pinned HPA v23 input, not an already revised
    screen. A reported zero and a failing donor must be different donors in the
    same gene/tissue/cohort/version. Technical duplicates cannot supply either.
    """
    columns = [
        "Gene ID",
        "Tissue",
        "Cohort",
        "Donors",
        "Reported zero",
        "Failing donors",
        "Highest RNA (nTPM)",
        "Absolute ceiling (nTPM)",
        "Selectivity ceiling (nTPM)",
        "Reproductive maximum (nTPM)",
        "Negative local IHC",
        "Variable pattern",
        "Donor findings",
        "RNA/IHC pairing",
        "Source",
    ]
    if donors.empty:
        return pd.DataFrame(columns=columns)
    required = {
        "Ensembl_Gene_ID",
        "tissue",
        "cohort",
        "source_version",
        "donor_id",
        "sample_id",
        "rna_ntpm",
        "unit",
        "taxonomy_id",
        "source_url",
        "rna_ihc_pairing",
    }
    if missing := required - set(donors):
        raise ValueError(f"Missing donor fields: {sorted(missing)}")
    if not (
        donors.source_version.eq("v23").all()
        and donors.unit.eq("HPA-v23-nTPM").all()
        and pd.to_numeric(donors.taxonomy_id).eq(9606).all()
    ):
        raise ValueError("Donor rules require human HPA v23 nTPM; do not substitute raw TPM")
    keys = ["Ensembl_Gene_ID", "tissue", "cohort", "source_version"]
    for key in [*keys, "donor_id", "sample_id", "source_url", "rna_ihc_pairing"]:
        if donors[key].isna().any() or donors[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing donor provenance: {key}")
    if donors.duplicated([*keys, "donor_id"]).any():
        raise ValueError("Duplicate donors within a cohort/gene/tissue")
    if donors.duplicated([*keys, "sample_id"]).any():
        raise ValueError("Duplicate RNA specimens within a cohort/gene/tissue")
    donors = donors.copy()
    if not donors.index.is_unique:
        raise ValueError("Donor observation indices must be unique")
    donors["rna_ntpm"] = pd.to_numeric(donors.rna_ntpm, errors="raise")
    if not (np.isfinite(donors.rna_ntpm).all() and donors.rna_ntpm.ge(0).all()):
        raise ValueError("Invalid donor RNA; missing measurements are not zeros")
    screened = screen_tissues(audit).set_index(["Ensembl_Gene_ID", "tissue"])
    records = []
    for (gid, tissue, cohort, _), group in donors.groupby(keys, sort=True):
        if (gid, tissue) not in screened.index:
            raise ValueError(f"Donor tissue has no matched screening context: {gid}, {tissue}")
        context = screened.loc[(gid, tissue)]
        if context.screening_tier != "vital":
            continue
        n = group.rna_ntpm
        r = float(context.reproductive_rna_max_ntpm)
        c = float(context.rna_ceiling_ntpm)
        limit = float(context.continuous_selectivity_ceiling_ntpm)
        if not np.isfinite([r, c, limit]).all():
            raise ValueError("Cannot evaluate donor RNA without usable tissue rules")
        failed = (n > np.nextafter(c, np.inf)) | (n > np.nextafter(limit, np.inf)) | (n >= r)
        local_negative = context.ihc_evidence_status in {
            "negative_approved",
            "negative_supported",
            "negative_enhanced",
        }
        zero = int(n.eq(0).sum())
        records.append(
            {
                "Gene ID": gid,
                "Tissue": tissue,
                "Cohort": cohort,
                "Donors": len(group),
                "Reported zero": zero,
                "Failing donors": int(failed.sum()),
                "Highest RNA (nTPM)": n.max(),
                "Absolute ceiling (nTPM)": c,
                "Selectivity ceiling (nTPM)": limit,
                "Reproductive maximum (nTPM)": r,
                "Negative local IHC": local_negative,
                "Variable pattern": bool(zero and failed.any() and local_negative),
                "Donor findings": "; ".join(
                    f"{row.sample_id}: {row.rna_ntpm:g} nTPM ({'fails' if failed.loc[i] else 'passes'})"
                    for i, row in group.iterrows()
                ),
                "RNA/IHC pairing": "; ".join(sorted(set(group.rna_ihc_pairing))),
                "Source": "; ".join(sorted(set(group.source_url))),
            }
        )
    return pd.DataFrame(records, columns=columns)


def screen_candidates(candidates, audit, *, donors=None):
    """Return gene decisions, tissue evidence and donor evidence.

    Retained means a research point-screen pass. Variable and Excluded are both
    excluded from selection. Variable never overrides positive protein staining,
    unresolved quality, a family exclusion, or a separate safety hold.
    """
    if donors is None:
        donors = reviewed_donor_rna()
        donors = donors.loc[donors.Ensembl_Gene_ID.isin(candidates.Ensembl_Gene_ID)]
    tissues = screen_tissues(audit)
    donor_results = screen_donors(audit, donors)
    groups = dict(iter(tissues.groupby("Ensembl_Gene_ID", sort=False)))
    donor_groups = dict(iter(donor_results.groupby("Gene ID", sort=False)))
    records = []
    clinical_evidence = candidates.set_index("Ensembl_Gene_ID")["Clinical evidence"]
    for gene in candidates.itertuples():
        rows = groups.get(gene.Ensembl_Gene_ID, tissues.iloc[:0])
        normal = rows.loc[rows.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)]
        reviewed = donor_groups.get(gene.Ensembl_Gene_ID, donor_results.iloc[:0])
        failure_codes = set(filter(None, gene.failure_reasons.split(";"))) - NUMERICAL
        if gene.external_safety_hold:
            failure_codes.add("external safety hold")
        if not gene.family_eligible:
            failure_codes.add("family exclusion")
        if not gene.nomination_review_eligible:
            failure_codes.add("nomination unresolved")
        if not gene.reproductive_positive:
            failure_codes.add("reproductive expression not demonstrated")
        if len(normal) and normal.disallowed_positive_ihc.any():
            failure_codes.add("positive normal IHC")
        quality = list(filter(None, gene.quality_flags.split(";")))
        missing_tissues = NORMAL_PANEL - set(normal.tissue)
        if missing_tissues:
            quality.append("missing normal RNA: " + ", ".join(sorted(missing_tissues)))
        if normal.empty or normal.point_result.eq("Unassessed").any():
            quality.append("incomplete tissue evidence")
        failures = normal.loc[normal.point_result.eq("Fail")]
        donor_failure = bool(len(reviewed) and reviewed["Failing donors"].gt(0).any())
        variable = bool(len(reviewed) and reviewed["Variable pattern"].any())
        # A negative annotation must exist; absence of IHC is never all-negative.
        measured = normal.loc[normal.ihc_evidence_status.str.startswith("negative_")]
        can_reconsider = bool(variable and len(measured) and not failure_codes and not quality)
        if failure_codes:
            status = "Excluded"
        elif can_reconsider:
            status = "Variable"
        elif donor_failure or len(failures):
            status = "Excluded"
        elif quality:
            status = "Unassessed"
        else:
            status = "Retained"
        findings = [x.replace("_", " ") for x in [*sorted(failure_codes), *quality]]
        findings += [
            f"{r.tissue}: {r.rna_ntpm:g} nTPM; {r.point_reasons.replace('_', ' ')}"
            for r in failures.itertuples()
        ]
        if donor_failure:
            findings.append("individual vital-tissue RNA donors fail the same inclusion limits")
        controls = sorted(set(rows.reproductive_ihc_control) - {""})
        records.append(
            {
                "Gene": gene.Symbol,
                "Gene ID": gene.Ensembl_Gene_ID,
                "Category": status,
                "Selected": status == "Retained",
                "Reason": "; ".join(findings) or "All available point-screen rules pass",
                "Reproductive IHC": "; ".join(controls) or "Not demonstrated",
                "Donor review": "Reviewed" if len(reviewed) else "Not assessed",
                "Donor findings": " | ".join(reviewed["Donor findings"]),
                "RNA/IHC pairing": " | ".join(reviewed["RNA/IHC pairing"]),
                "Reconsideration": (
                    "Determine why donor RNA differs; obtain paired RNA/protein in RNA-positive "
                    "and RNA-zero vital specimens, establish cell/isoform identity, and validate "
                    "any proposed subgroup before an explicit new review. RNA zero alone does "
                    "not establish eligibility."
                )
                if can_reconsider
                else "",
                "Clinical evidence": clinical_evidence.loc[gene.Ensembl_Gene_ID],
                "Policy": POLICY_ID,
            }
        )
    return pd.DataFrame(records).sort_values("Gene"), tissues, donor_results
