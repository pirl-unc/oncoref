"""Reviewed TCR-T experience and its limited protein-level RNA screening credit.

Clinical evidence belongs to a peptide, HLA, receptor, regimen and population.
The credit is a research-policy inference to a protein's aggregated RNA profile,
not clinical clearance of every isoform, peptide, receptor or patient. Sentinel
events are not an exhaustive adverse-event table. Unreviewed is never negative.
"""

from __future__ import annotations

import json
from copy import deepcopy
from importlib.resources import files

import pandas as pd

EVIDENCE_FILE = "tcr-trial-evidence.json"
GOLD_STANDARD_FILE = "tcr-clinical-gold-standard.csv"
SPECIFICITY_FINDINGS = {"none_identified", "suspected_or_confirmed", "unassessed"}
SPECIFICITY_ATTRIBUTIONS = {
    "suspected",
    "confirmed",
    "not_supported",
    "not_established",
    "unresolved",
    "unassessed",
}
CLINICAL_STRATA = {
    "safe_cta_potency_proven",
    "toxic_on_target",
    "toxic_cross_reactivity",
    "conservative_gene_level_unsafe",
    "excluded_not_potency_proven",
}
TOXICITY_MECHANISMS = {
    "none_attributed",
    "on_target_off_tumor",
    "paralog_or_off_target_cross_reactivity",
    "unproven_potency",
}


def validate(record: dict) -> None:
    """Reject ambiguous identities, denominators and contradictory safety reviews."""
    if (record.get("schema"), record.get("schema_version")) != ("oncoref.tcr-trial-evidence", 1):
        raise ValueError("Unsupported TCR trial evidence schema")
    products = {p["product_id"]: p for p in record["products"]}
    cohorts = {c["cohort_id"]: c for c in record["cohorts"]}
    if len(products) != len(record["products"]) or len(cohorts) != len(record["cohorts"]):
        raise ValueError("Duplicate product or cohort identity")
    if len({e["event_id"] for e in record["events"]}) != len(record["events"]):
        raise ValueError("Duplicate event identity")
    sources = record["sources"]
    for p in products.values():
        if p["modality"] != "TCR-T" or p["taxon_id"] != 9606:
            raise ValueError("This curation is scoped to human TCR-T")
        if not p["hla"] or not p["targets"] or not p["identity_sources"]:
            raise ValueError("Missing peptide/HLA/product identity")
        if not set(p["identity_sources"]) <= sources.keys():
            raise ValueError("Unknown identity source")
        for target in p["targets"]:
            if not all(target.get(k) for k in ["gene_id", "symbol", "peptide", "isoform_scope"]):
                raise ValueError("Incomplete target identity")
    for c in cohorts.values():
        if c["product_id"] not in products or c["source"] not in sources:
            raise ValueError("Unknown cohort product/source")
        n = c["patients_treated"]
        if n is not None and (type(n) is not int or n <= 0):
            raise ValueError("Treated denominator must be positive or unknown")
        if c["grade45_specificity"] not in SPECIFICITY_FINDINGS:
            raise ValueError("Unknown specificity finding")
        if c["review_status"] not in {"reviewed", "pending"}:
            raise ValueError("Unknown review status")
        if c["grade45_specificity"] == "none_identified" and (
            c["review_status"] != "reviewed" or n is None or not c["source_sections"]
        ):
            raise ValueError("Favorable evidence requires reviewed results and denominator")
    for e in record["events"]:
        if e["product_id"] not in products or e["source"] not in sources:
            raise ValueError("Unknown event product/source")
        if e["specificity_attribution"] not in SPECIFICITY_ATTRIBUTIONS:
            raise ValueError("Unknown specificity attribution")
        for finding in e.get("causal_findings", []):
            if finding["source"] not in sources or finding["interpretation_status"] not in {
                "accepted",
                "unresolved",
                "rejected",
            }:
                raise ValueError("Invalid causal finding source or interpretation")
            if not all(
                finding.get(k) for k in ["protein", "peptide", "hla", "evidence", "limitations"]
            ):
                raise ValueError("Causal findings require peptide, HLA and evidence scope")
        if e["grade"] is not None and (type(e["grade"]) is not int or not 1 <= e["grade"] <= 5):
            raise ValueError("Grade must be 1-5 or unreported")
        cohort = e["cohort_id"]
        if cohort is not None:
            if cohort not in cohorts or cohorts[cohort]["product_id"] != e["product_id"]:
                raise ValueError("Event/cohort product mismatch")
            if (
                _specificity_block(e)
                and cohorts[cohort]["grade45_specificity"] == "none_identified"
            ):
                raise ValueError("Favorable cohort contradicts a specificity event")


def _specificity_block(event: dict) -> bool:
    # An ungraded suspected specificity injury cannot masquerade as below grade 4.
    return event["specificity_attribution"] in {"suspected", "confirmed"} and (
        event["grade"] is None or event["grade"] >= 4
    )


def evidence() -> dict:
    """Return a validated, detached snapshot with linked primary sources, offline."""
    record = json.loads(
        files("oncoref").joinpath("data", EVIDENCE_FILE).read_text(encoding="utf-8")
    )
    validate(record)
    return record


def products(*, record: dict | None = None) -> pd.DataFrame:
    """One row per reviewed receptor/product, never one observation per gene alias.

    Favorable credit needs a reviewed treated cohort with no identified grade 4-5
    specificity signal, and no such signal in any reviewed cohort/event for that
    receptor. A pending cohort cannot supply credit; an explicitly unassessed
    specificity event blocks credit until reviewed. Treatment-related CRS,
    conditioning toxicity and events with no established specificity mechanism
    remain visible without being renamed antigen-specific tissue destruction.
    """
    record = evidence() if record is None else deepcopy(record)
    validate(record)
    rows = []
    for p in record["products"]:
        cohorts = [c for c in record["cohorts"] if c["product_id"] == p["product_id"]]
        events = [e for e in record["events"] if e["product_id"] == p["product_id"]]
        positive = any(c["grade45_specificity"] == "suspected_or_confirmed" for c in cohorts)
        blocked = positive or any(_specificity_block(e) for e in events)
        unresolved = any(
            e["specificity_attribution"] in {"unassessed", "unresolved"} for e in events
        )
        favorable = [
            c["cohort_id"]
            for c in cohorts
            if c["review_status"] == "reviewed" and c["grade45_specificity"] == "none_identified"
        ]
        credit = bool(favorable) and not blocked and not unresolved
        rows.append(
            p
            | {
                "rna_credit": credit,
                "credit_cohorts": tuple(sorted(favorable)) if credit else (),
                "clinical_finding": "specificity_signal"
                if blocked
                else "unassessed"
                if unresolved or not favorable
                else "no_identified_severe_specificity_signal",
                "event_ids": tuple(e["event_id"] for e in events),
            }
        )
    return pd.DataFrame(rows)


def cohorts() -> pd.DataFrame:
    """Reported cohort denominators; do not add repeated publications or cases."""
    return pd.DataFrame(evidence()["cohorts"])


def events() -> pd.DataFrame:
    """Sentinel clinical events, including ungraded and unresolved attribution."""
    return pd.DataFrame(evidence()["events"])


def gene_exclusions() -> pd.DataFrame:
    """Explicit gene consequences, separate from product injury and credit.

    An unresolved causal peptide is context; independent normal-cell evidence
    can still exclude its gene through the appropriate normal-evidence gate.
    """
    return pd.DataFrame(
        evidence().get("gene_exclusions", []),
        columns=["gene_id", "symbol", "effect", "event_id", "rationale"],
    )


def protein_credit(gene_ids, *, taxon_id: int = 9606, record: dict | None = None) -> dict:
    """Assess modest protein-profile credit; every member locus must be supported.

    No symbol-prefix/family/orthology transfer. Multiple trials never stack.
    Returning credit does not clear other peptides, untested isoforms or a new
    receptor. Use ``products`` and its exact target mappings for peptide work.
    """
    ids = set(gene_ids)
    if not ids or isinstance(gene_ids, str):
        raise ValueError("Supply a nonempty collection of gene IDs")
    if taxon_id != 9606:
        return {
            "rna_credit": False,
            "product_ids": (),
            "cohort_ids": (),
            "scope": "no_species_transfer",
        }
    p = products(record=record)
    selected = []
    covered = set()
    for row in p.to_dict("records"):
        matches = ids & {t["gene_id"] for t in row["targets"]}
        if row["rna_credit"] and matches:
            selected.append(row)
            covered |= matches
    credit = ids <= covered
    return {
        "rna_credit": credit,
        "product_ids": tuple(sorted(row["product_id"] for row in selected)) if credit else (),
        "cohort_ids": tuple(sorted({c for row in selected for c in row["credit_cohorts"]}))
        if credit
        else (),
        "scope": "protein_profile_only_tested_epitope_limits_apply"
        if credit
        else "unassessed_or_not_qualified",
    }


def annotate_proteins(facts: pd.DataFrame) -> pd.DataFrame:
    """Attach credit and provenance to the replay's reference-bound member sets."""
    record = evidence()
    supported = {t["gene_id"] for p in record["products"] for t in p["targets"]}
    d = facts.copy()
    credits = []
    for members in d.member_gene_ids:
        ids = json.loads(members) if isinstance(members, str) else members
        credits.append(
            protein_credit(ids, record=record)
            if supported.intersection(ids)
            else {"rna_credit": False, "product_ids": (), "cohort_ids": (), "scope": "unreviewed"}
        )
    d["clinical_rna_credit"] = [c["rna_credit"] for c in credits]
    d["clinical_product_ids"] = [";".join(c["product_ids"]) for c in credits]
    d["clinical_cohort_ids"] = [";".join(c["cohort_ids"]) for c in credits]
    d["clinical_credit_scope"] = [c["scope"] for c in credits]
    return d


def _parse_strict_bool_series(series: pd.Series, col_name: str) -> pd.Series:
    parsed = []
    for val in series.tolist():
        if isinstance(val, bool) or type(val).__name__ == "bool_":
            parsed.append(bool(val))
        elif isinstance(val, str) and val.strip() in {"True", "False"}:
            parsed.append(val.strip() == "True")
        else:
            raise ValueError(f"Column {col_name!r} requires strict booleans, got {val!r}")
    return pd.Series(parsed, index=series.index, dtype=bool)


def validate_gold_standard(df: pd.DataFrame) -> None:
    """Validate the potency-proven clinical TCR-T / ImmTAC gold-standard benchmark."""
    required = {
        "benchmark_id",
        "target_symbol",
        "entity_id",
        "ensembl_gene_ids",
        "cross_recognized_gene_ids",
        "product_name",
        "linked_product_id",
        "modality",
        "hla_restriction",
        "peptide",
        "patients_treated",
        "orr_pct",
        "potency_proven",
        "clinical_stratum",
        "toxicity_mechanism",
        "severe_specificity_toxicity",
        "observed_clinical_toxicities",
        "strongest_proxy_assay",
        "strongest_proxy_ratio",
        "citation",
        "source_url",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Gold-standard table missing required columns: {sorted(missing)}")
    if not df["benchmark_id"].is_unique or df["benchmark_id"].isna().any():
        raise ValueError("Duplicate or missing benchmark_id in clinical gold standard")
    if not set(df["clinical_stratum"]) <= CLINICAL_STRATA:
        raise ValueError("Unknown clinical_stratum in clinical gold standard")
    if not set(df["toxicity_mechanism"]) <= TOXICITY_MECHANISMS:
        raise ValueError("Unknown toxicity_mechanism in clinical gold standard")
    for bcol in ("potency_proven", "severe_specificity_toxicity"):
        for val in df[bcol].tolist():
            if not (isinstance(val, bool) or type(val).__name__ == "bool_"):
                raise ValueError(f"Column {bcol!r} must contain boolean values, got {val!r}")

    toxic_strata = {"toxic_on_target", "toxic_cross_reactivity"}
    known_products = {p["product_id"] for p in evidence()["products"]}
    for r in df.itertuples(index=False):
        if (
            pd.notna(r.linked_product_id)
            and str(r.linked_product_id).strip()
            and str(r.linked_product_id) not in known_products
        ):
            raise ValueError(f"Unknown linked_product_id {r.linked_product_id!r}")
        if bool(r.severe_specificity_toxicity) != (r.clinical_stratum in toxic_strata):
            raise ValueError(
                f"Row {r.benchmark_id!r} has inconsistent clinical_stratum={r.clinical_stratum!r} "
                f"and severe_specificity_toxicity={r.severe_specificity_toxicity!r}"
            )
        if bool(r.potency_proven):
            if r.clinical_stratum == "excluded_not_potency_proven":
                raise ValueError(f"Potency-proven row {r.benchmark_id!r} cannot be excluded")
            if pd.isna(r.patients_treated) or int(r.patients_treated) <= 0:
                raise ValueError(
                    f"Potency-proven row {r.benchmark_id!r} needs patients_treated > 0"
                )
            # Potency requires either objective clinical tumor responses (ORR > 0) or
            # rapid lethal in vivo cytolysis that precluded radiologic staging (magea3_a3a).
            has_orr = pd.notna(r.orr_pct) and float(r.orr_pct) > 0.0
            if not has_orr and not bool(r.severe_specificity_toxicity):
                raise ValueError(
                    f"Row {r.benchmark_id!r} marked potency_proven without ORR > 0 or lethal in vivo lysis"
                )
        else:
            if r.clinical_stratum != "excluded_not_potency_proven":
                raise ValueError(
                    f"Non-potency-proven row {r.benchmark_id!r} must use excluded_not_potency_proven"
                )


def gold_standard(*, potency_proven_only: bool = False) -> pd.DataFrame:
    """Return curated endogenous-protein TCR-T / ImmTAC trials with response & toxicity outcomes.

    Human clinical TCR trials with confirmed anti-tumor cytolytic potency (``ORR > 0%``
    or rapid *in vivo* target-cell lysis) serve as the empirical toxicity gold standard
    against which OncoRef's healthy-tissue RNA, single-cell, IHC, rare-cell UMI, and
    1-mismatch sequence screens act as proxies. Trials with ``0% ORR`` (e.g. ``IMA202``)
    or pending clinical readouts (e.g. ``CT83``) remain visible with
    ``potency_proven == False`` and are excluded when ``potency_proven_only=True``.
    """
    from io import StringIO

    raw = files("oncoref").joinpath("data", GOLD_STANDARD_FILE).read_text(encoding="utf-8")
    df = pd.read_csv(
        StringIO(raw), dtype={"potency_proven": object, "severe_specificity_toxicity": object}
    )
    df["potency_proven"] = _parse_strict_bool_series(df["potency_proven"], "potency_proven")
    df["severe_specificity_toxicity"] = _parse_strict_bool_series(
        df["severe_specificity_toxicity"], "severe_specificity_toxicity"
    )
    validate_gold_standard(df)
    if potency_proven_only:
        df = df[df["potency_proven"]].reset_index(drop=True)
    return df


def evaluate_against_gold_standard(
    catalog_df: pd.DataFrame | None = None, *, potency_proven_only: bool = False
) -> pd.DataFrame:
    """Evaluate OncoRef's protein-level safety proxies against the clinical TCR gold standard."""
    from . import target_selection

    cat = target_selection.results() if catalog_df is None else catalog_df
    by_eid = cat.set_index("entity_id", drop=False)
    by_gid: dict[str, pd.Series] = {}
    for row in cat.itertuples(index=False):
        s = pd.Series(row._asdict())
        for gid in json.loads(row.member_gene_ids):
            by_gid[gid] = s

    gs = gold_standard(potency_proven_only=potency_proven_only).copy()
    cats, classes, normal_tpms, unsafe_reasons = [], [], [], []
    primary_unsafe_col, cross_cats, cross_unsafe_col, any_unsafe, alignments = [], [], [], [], []

    for r in gs.itertuples(index=False):
        if r.entity_id not in by_eid.index:
            raise ValueError(f"Benchmark entity_id {r.entity_id!r} missing from catalog")
        prow = by_eid.loc[r.entity_id]
        member_gids = set(json.loads(str(prow["member_gene_ids"])))
        declared_gids = {g.strip() for g in str(r.ensembl_gene_ids).split(";") if g.strip()}
        if not declared_gids or not (declared_gids <= member_gids):
            raise ValueError(
                f"Benchmark {r.benchmark_id!r} ensembl_gene_ids {declared_gids} mismatch {member_gids}"
            )
        cats.append(str(prow["category"]))
        classes.append(str(prow["classification"]))
        normal_tpms.append(
            float(prow["normal_reference_tpm"])
            if pd.notna(prow["normal_reference_tpm"])
            else float("nan")
        )
        unsafe_reasons.append(
            str(prow["Unsafe_reasons"]) if pd.notna(prow["Unsafe_reasons"]) else ""
        )

        x_gids = [
            g.strip()
            for g in (
                str(r.cross_recognized_gene_ids) if pd.notna(r.cross_recognized_gene_ids) else ""
            ).split(";")
            if g.strip()
        ]
        unknown_xgids = [g for g in x_gids if g not in by_gid]
        if unknown_xgids:
            raise ValueError(
                f"Benchmark {r.benchmark_id!r} unknown cross_recognized_gene_ids: {unknown_xgids}"
            )
        x_cat_list = [str(by_gid[g]["category"]) for g in x_gids]
        cross_cats.append(";".join(x_cat_list))

        prim_unsafe = str(prow["category"]) == "Unsafe"
        x_unsafe = "Unsafe" in x_cat_list
        flagged_unsafe = prim_unsafe or x_unsafe
        primary_unsafe_col.append(prim_unsafe)
        cross_unsafe_col.append(x_unsafe)
        any_unsafe.append(flagged_unsafe)

        if not bool(r.potency_proven):
            alignments.append("excluded_potency_unproven")
        elif bool(r.severe_specificity_toxicity):
            alignments.append("true_unsafe_toxic" if flagged_unsafe else "false_safe_miss")
        elif str(prow["category"]) == "CTA":
            alignments.append("true_safe_cta")
        elif flagged_unsafe:
            alignments.append("conservative_unsafe_epitope_tolerated")
        else:
            alignments.append("other_held")

    gs["oncoref_category"] = cats
    gs["oncoref_classification"] = classes
    gs["oncoref_normal_reference_tpm"] = normal_tpms
    gs["oncoref_unsafe_reasons"] = unsafe_reasons
    gs["primary_target_unsafe"] = primary_unsafe_col
    gs["cross_recognized_categories"] = cross_cats
    gs["cross_recognized_unsafe"] = cross_unsafe_col
    gs["effective_oncoref_unsafe"] = any_unsafe
    gs["benchmark_alignment"] = alignments
    return gs


def gold_standard_summary(catalog_df: pd.DataFrame | None = None) -> dict:
    """Return quantitative proxy-calibration metrics against potency-proven human TCR trials."""
    ev = evaluate_against_gold_standard(catalog_df=catalog_df, potency_proven_only=False)
    potent = ev[ev["potency_proven"]]
    toxic = potent[potent["severe_specificity_toxicity"]]
    safe_cta = potent[potent["clinical_stratum"].eq("safe_cta_potency_proven")]
    conservative = potent[potent["benchmark_alignment"].eq("conservative_unsafe_epitope_tolerated")]
    excluded = ev[~ev["potency_proven"]]

    n_toxic = len(toxic)
    n_toxic_primary = int(toxic["primary_target_unsafe"].sum())
    n_toxic_caught = int(toxic["effective_oncoref_unsafe"].sum())
    n_safe_cta = len(safe_cta)
    n_safe_cta_kept = int(safe_cta["oncoref_category"].eq("CTA").sum())

    return {
        "total_entries": len(ev),
        "potency_proven_entries": len(potent),
        "excluded_non_potent_or_pending_entries": len(excluded),
        "severe_specificity_toxic_entries": n_toxic,
        "severe_specificity_toxic_primary_unsafe": n_toxic_primary,
        "severe_specificity_toxic_flagged_unsafe": n_toxic_caught,
        "false_safe_misses": int((toxic["benchmark_alignment"] == "false_safe_miss").sum()),
        "primary_target_sensitivity_to_toxicity": (n_toxic_primary / n_toxic) if n_toxic else 0.0,
        "sensitivity_to_clinical_toxicity": (n_toxic_caught / n_toxic) if n_toxic else 0.0,
        "potency_proven_safe_cta_entries": n_safe_cta,
        "potency_proven_safe_cta_retained": n_safe_cta_kept,
        "recall_on_safe_ctas": (n_safe_cta_kept / n_safe_cta) if n_safe_cta else 0.0,
        "conservative_non_cta_unsafe_entries": len(conservative),
        "conservative_non_cta_targets": sorted(set(conservative["target_symbol"].astype(str))),
        "safe_cta_targets": sorted(set(safe_cta["target_symbol"].astype(str))),
        "toxic_targets": sorted(set(toxic["target_symbol"].astype(str))),
    }
