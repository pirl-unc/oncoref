#!/usr/bin/env python3
"""Build shared species audit tables from the current pinned selection replay."""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import cta_selection, reference_data  # noqa: E402
from oncoref.species_antigens import cta_snapshot  # noqa: E402

RELEASE = ROOT / "docs/audits/cta-selection-20261028"
OUTPUT = ROOT / "docs/audits/cta-source-walkthrough-20261028"
DATA = ROOT / "oncoref/data"
BASE_URL = "https://pirl-unc.github.io/oncoref/audits/cta-source-walkthrough-20261028/"
COMPARATORS = (
    "PRAME",
    "MAGEA1",
    "MAGEA3",
    "MAGEA4",
    "MAGEA12",
    "NY-ESO-1",
    "CTAG2",
    "XAGE1A/B",
    "SSX2/B",
    "CT83",
    "TBXT",
    "SMC1B",
    "PAGE1",
)
TABLE_COLUMNS = {
    "sources": [
        "source_id",
        "citation",
        "source_url",
        "entity_kind",
        "mapped_entities",
        "mapped_gene_loci",
        "CTA",
        "Unsafe",
        "Other",
        "Unassigned",
        "unique_entities",
        "primary_source",
        "covered_by",
        "evidence_scope",
        "status",
    ],
    "memberships": ["source_id", "entity_id", "label", "gene_ids", "category"],
    "overlaps": ["source_a", "source_b", "mapped_overlap", "selected_overlap"],
    "funnel": [
        "stage",
        "status",
        "candidate_removals",
        "eligible",
        "Unsafe",
        "Other",
        "Unassigned",
        "note",
    ],
    "criteria": ["stage", "status", "rule", "source_ids"],
    "explanations": [
        "entity_id",
        "label",
        "gene_ids",
        "category",
        "layer",
        "decision",
        "tissue",
        "cell_type",
        "value",
        "unit",
        "candidate_ceiling",
        "baseline_ceiling",
        "clinical_product_ids",
        "clinical_credit_applied",
        "selection_tier",
        "ihc_coverage",
        "source_id",
        "source_version",
        "source_context",
        "reason",
        "source_url",
    ],
}
HUMAN_RULES = {
    "Reproductive or developmental expression": (
        "Full Reproductive bulk ≥5 nTPM OR supported reproductive-cell RNA ≥5 nTPM (≥20 reads; indirect sites additionally need complete outside-cluster coverage and every outside value <5, with no curated normal-induction conflict) OR reliable mapped positive reproductive IHC for every member without an antibody conflict OR supported developmental stage/lineage mean ≥5 TPM; ≥20 cells, ≥5 positive and ≥10% positive at ≥1 TPM. Presence evidence never overrides normal exclusions.",
        [
            "hpa_rna_consensus",
            "hpa_normal_tissue",
            "hpa_single_cell_tissue",
            "normal-induction-evidence",
            "tyser_E-MTAB-9388",
            "petropoulos_E-MTAB-3929",
        ],
    ),
    "Cancer expression": (
        "At least one measured compatible tumor sample ≥2 in native TPM or HPA pTPM OR reviewed reliable tumor protein; do not add maxima from different samples.",
        ["oncoref_tumor_union", "hpa_cancer"],
    ),
    "HPA normal protein": (
        "Reliable non-Low positives in Vital/Protected supply Unsafe evidence; Manageable positives block CTA only; uncertain evidence remains unresolved. Low is flagged.",
        ["hpa_normal_tissue"],
    ),
    "Reproductive concentration": (
        "HPA v23 average normal/(normal+Reproductive maximum) above 10% blocks the adult/placental route when its reference is positive and the comparison is complete. Use 38 contexts/29 anatomy units, including Disposable and excluding Permissive thymus. Supported embryonic evidence remains independent; its contrast is unassessed. This is not safety evidence.",
        ["hpa_rna_consensus"],
    ),
    "Required RNA coverage": (
        "Complete GTEx v11 bulk core (52 native contexts), required nonreproductive bulk anatomy across compatible sources, and distinct HPA tissue-cell/immune panels; every protein member measured. DICE and cornea coverage is optional.",
        [
            "hpa_rna_consensus",
            "hpa_brain",
            "gtex_tissue_rna",
            "hpa_single_cell_tissue",
            "hpa_immune_cell",
            "hpa_immune_monaco",
        ],
    ),
    "Normal bulk": (
        "One inclusive tissue/IHC ceiling: above it in Vital/Protected supplies Unsafe evidence; above it in Manageable blocks CTA. No second RNA boundary.",
        ["hpa_rna_consensus", "hpa_brain", "gtex_tissue_rna"],
    ),
    "HPA tissue clusters": (
        "RNA requires ≥20 reads and uses the single ceiling. IHC credit requires an exact or declared compartment match; unrelated cells in the same tissue cannot inherit it.",
        ["hpa_single_cell_tissue"],
    ),
    "HPA immune": (
        "Protected-tier native nTPM limits without matched negative-IHC credit; reproductive/placental location remains neutral.",
        ["hpa_immune_cell"],
    ),
    "Monaco immune": (
        "Protected-tier native TPM limits without matched negative-IHC credit.",
        ["hpa_immune_monaco"],
    ),
    "DICE immune": (
        "Third-highest reported TPM value; native protected-tier limits. Unpaired group distributions retain bounds. Missing coverage alone does not exclude.",
        ["dice"],
    ),
    "Adult cornea": (
        "Recurrent native CP10K: ≥5 positive cells and ≥20 counts in ≥2 sample groups. Nominal RNA limits /100; not a TPM conversion. Missing coverage alone does not exclude.",
        ["GSE155683_adult"],
    ),
    "Accepted normal-protein evidence": (
        "Accepted positive protein findings exclude unless extracted Low intensity retains a flag; unquantified does not mean negative.",
        ["cta_policy_findings"],
    ),
    "Clinical safety exclusion": (
        "Retain reviewed product-specific clinical injury counterevidence and causal-attribution limits; no admission shortcut from clinical use.",
        ["cta_policy_findings"],
    ),
    "External validated normal protein": (
        "Use categorical/quantitative staining only within verified reagent, tissue and cell scope; no inherited HPA reliability label.",
        ["external_ihc"],
    ),
}


def records(frame: pd.DataFrame) -> list[dict]:
    """Convert pandas nulls/scalars to portable JSON values."""
    return json.loads(frame.to_json(orient="records"))


def citation(tag: str) -> str:
    """Use one human-readable author/year style; source keys remain separate."""
    match = re.match(r"([A-Za-z]+)(\d{4})", tag)
    if not match:
        raise ValueError(f"No author/year in source tag: {tag}")
    author = {"daSilva": "da Silva", "Boel": "Boël"}.get(match[1], match[1])
    return f"{author} et al. ({match[2]})"


def display(value: object) -> str:
    """Render missing evidence explicitly, and numbers with useful precision."""
    if pd.isna(value):
        return "unavailable"
    return f"{value:.6g}" if isinstance(value, (float, int)) else str(value)


def membership_tables(species: str) -> dict[str, list[dict]]:
    """Union tables within each paper, then count distinct entities and loci."""
    members, info, scopes = [], {}, {}
    if species == "human":
        registry = pd.read_csv(DATA / "cta-publication-sources.csv").set_index("source_tag")
        known = pd.read_csv(RELEASE / "known-cta-decisions.csv.gz", low_memory=False)
        for doi, rows in registry.groupby("doi"):
            info[doi] = (citation(rows.index[0]), "https://doi.org/" + doi)
            scopes[doi] = "; ".join(sorted(set(rows.evidence_scope)))
        for row in known.dropna(subset=["membership_sources"]).itertuples():
            for tag in row.membership_sources.split(";"):
                members.append(
                    {
                        "source_id": registry.loc[tag, "doi"],
                        "entity_id": row.protein_equivalence_group_id,
                        "label": row.canonical_symbol
                        if pd.notna(row.canonical_symbol)
                        else row.gene_id,
                        "gene_id": row.gene_id,
                        "category": row.protein_group_category,
                    }
                )
        # Author-recognized prior CT antigens annotate the paper without being
        # misrepresented as new discoveries or quantitative admission evidence.
        context = pd.read_csv(DATA / "cta-publication-context.csv")
        current = cta_selection.results().set_index("entity_id")
        # New publication sets use the common mapped-intake contract. The
        # frozen known-CTA comparison remains a separate historical cohort.
        registered = {tag for tags in known.membership_sources.dropna() for tag in tags.split(";")}
        intake = pd.read_csv(DATA / "cta-publication-membership.csv").fillna("")
        gene_to_entity = {g: entity for entity, row in current.iterrows() for g in row.gene_ids}
        for row in intake[~intake.source_tag.isin(registered)].itertuples():
            entity = gene_to_entity.get(row.Ensembl_Gene_ID)
            if row.mapping_status != "mapped" or entity is None:
                continue
            target = current.loc[entity]
            members.append(
                {
                    "source_id": registry.loc[row.source_tag, "doi"],
                    "entity_id": entity,
                    "label": target.label,
                    "gene_id": row.Ensembl_Gene_ID,
                    "category": target.category,
                }
            )
        for row in context.itertuples():
            target = current.loc[row.entity_id]
            if set(json.loads(row.member_gene_ids)) != set(target.gene_ids):
                raise ValueError(f"Stale publication protein-group identity: {row.entity_id}")
            for gene_id in target.gene_ids:
                members.append(
                    {
                        "source_id": row.source_id,
                        "entity_id": row.entity_id,
                        "label": row.label,
                        "gene_id": gene_id,
                        "category": target.category,
                    }
                )
        frame = pd.DataFrame(members)
        grouped = (
            frame.groupby(["source_id", "entity_id"], sort=True)
            .agg(
                label=("label", lambda x: "/".join(sorted(set(x)))),
                gene_ids=("gene_id", lambda x: sorted(set(x))),
                category=("category", "first"),
            )
            .reset_index()
        )
        kind, status = "protein_equivalence_group", "annotation"
    else:
        snapshot = cta_snapshot(species)
        native = snapshot.as_dict()
        for row in snapshot.candidates():
            for source in row["derivation"]:
                members.append(
                    {
                        "source_id": source["source_id"],
                        "entity_id": f"nomination:{native['taxon_id']}:{row['id']}",
                        "label": row["label"],
                        "gene_ids": row["gene_ids"],
                        "category": None,
                    }
                )
        grouped = pd.DataFrame(members, columns=TABLE_COLUMNS["memberships"]).drop_duplicates(
            ["source_id", "entity_id"]
        )
        info = (
            {
                "hindriks2025": (
                    "Hindriks et al. (2025)",
                    "https://doi.org/10.1016/j.vetimm.2025.110946",
                ),
                "nemec2019": ("Nemec et al. (2019)", "https://doi.org/10.1111/vco.12475"),
            }
            if species == "dog"
            else {
                source["id"]: (f"PMID:{source['pmid']}", source["url"])
                for source in native["sources"]
            }
        )
        kind, status = "native_nomination", "unassessed"
        scopes = {source["id"]: source["scope"] for source in native["sources"]}
    all_sets = {sid: set(grouped.loc[grouped.source_id.eq(sid), "entity_id"]) for sid in info}
    selected = {
        sid: set(grouped.loc[grouped.source_id.eq(sid) & grouped.category.eq("CTA"), "entity_id"])
        for sid in info
    }
    # Only retire a single-digit list when retained sources cover its entire
    # mapped set. Sequential decisions preserve the union even for mutually
    # covering small lists. Provenance memberships always remain available.
    primary = {sid for sid, values in all_sets.items() if len(values) >= 10}
    for sid in sorted(set(info) - primary, key=lambda x: (-len(all_sets[x]), x)):
        covered = set().union(*(all_sets[x] for x in primary))
        if species != "human" or not all_sets[sid] <= covered:
            primary.add(sid)
    counts = []
    for sid, (label, url) in info.items():
        part = grouped[grouped.source_id.eq(sid)]
        others = set().union(*(values for other, values in all_sets.items() if other != sid))
        counts.append(
            {
                "source_id": sid,
                "citation": label,
                "source_url": url,
                "entity_kind": kind,
                "mapped_entities": len(part),
                "mapped_gene_loci": len({gene for ids in part.gene_ids for gene in ids}),
                **{c: int(part.category.eq(c).sum()) for c in cta_selection.CATEGORIES},
                "Unassigned": int(part.category.isna().sum()),
                "unique_entities": len(all_sets[sid] - others),
                "primary_source": sid in primary,
                "covered_by": sorted(x for x in primary if all_sets[x] & all_sets[sid])
                if sid not in primary
                else [],
                "evidence_scope": scopes[sid],
                "status": status,
            }
        )
    overlaps = [
        {
            "source_a": a,
            "source_b": b,
            "mapped_overlap": len(all_sets[a] & all_sets[b]),
            "selected_overlap": len(selected[a] & selected[b]) if species == "human" else None,
        }
        for a in info
        for b in info
    ]
    return {
        "sources": sorted(counts, key=lambda x: -x["mapped_entities"]),
        "memberships": records(grouped),
        "overlaps": overlaps,
    }


HUMAN_RULES["CTA and Other assignments"] = (
    "Other is the complete-normal-panel, resolved-evidence fallback when every normal gate passes but no candidate subtype qualifies. Tumor or developmental missingness alone does not block Other. Other:Tumor-supported requires low measured normal/developmental expression and two distinct tumor samples ≥10 TPM or reviewed reliable tumor protein; it is not target admission. Manageable holds remain unassigned, never cleared by breadth.",
    [
        "hpa_rna_consensus",
        "hpa_single_cell_tissue",
        "tyser_E-MTAB-9388",
        "oncoref_tumor_union",
        "hpa_normal_tissue",
    ],
)


def observation_source(row: dict | pd.Series) -> dict:
    """Retain the actual assay/version and native context in every RNA reason."""
    source = str(row.get("source", ""))
    version = row.get("version")
    version = str(version) if pd.notna(version) else None
    native = row.get("source_context")
    native = str(native) if pd.notna(native) else None
    if source in reference_data.REFERENCE_SOURCES and version:
        url = reference_data.REFERENCE_SOURCES[source]["urls"].get(version, BASE_URL)
    elif source.startswith("hpa"):
        url = "https://v23.proteinatlas.org/about/download"
    elif source.startswith("dice"):
        url = "https://dice-database.org/"
    elif "cornea" in source.lower() or "GSE155683" in source:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683"
    else:
        url = BASE_URL
    return {
        "source_id": source or None,
        "source_version": version,
        "source_context": native,
        "source_url": url,
    }


def human_explanations(facts: pd.DataFrame) -> pd.DataFrame:
    """Explain every contributing layer with a measured context when available."""
    evidence_fields = {
        "gene_id",
        "assay",
        "source",
        "version",
        "source_context",
        "tissue",
        "cell_type",
        "value",
        "unit",
        "candidate_ceiling",
        "baseline_ceiling",
        "clinical_product_ids",
        "clinical_credit_applied",
        "selection_tier",
        "ihc_match_scope",
        "negative_ihc_grade",
        "Read count",
        "value_upper",
        "ratio",
        "candidate_exclusion",
        "unsafe_evidence",
    }
    peaks = pd.read_csv(
        RELEASE / "normal-peak-evidence.csv.gz",
        low_memory=False,
        usecols=lambda column: column in evidence_fields,
    )
    proteins = pd.read_csv(RELEASE / "hpa-protein-evidence.csv.gz", low_memory=False)
    imported = pd.read_csv(RELEASE / "imported-protein-review-decisions.csv")
    coverage = pd.read_csv(OUTPUT / "normal-panel-coverage-gaps.csv.gz")
    coverage_by_entity = dict(iter(coverage.groupby("entity_id")))
    combined_bulk = pd.read_csv(OUTPUT / "combined-rna-coverage.csv.gz").set_index("entity_id")
    rows = []
    peak_index, protein_index = {}, {}
    for unsafe, column in [(False, "candidate_exclusion"), (True, "unsafe_evidence")]:
        best = (
            peaks[peaks[column]]
            .sort_values("ratio", ascending=False)
            .drop_duplicates(["gene_id", "assay"])
        )
        peak_index.update(
            {(row["gene_id"], row["assay"], unsafe): row for row in best.to_dict("records")}
        )
        best_protein = proteins[proteins[column]].drop_duplicates("representative_gene_id")
        protein_index.update(
            {
                (row["representative_gene_id"], unsafe): row
                for row in best_protein.to_dict("records")
            }
        )
    hold_columns = [column for column in facts if column.startswith("hold:")]
    for r in facts.to_dict("records"):
        base = {
            "entity_id": r["protein_equivalence_group_id"],
            "label": r["protein_group_label"],
            "gene_ids": json.loads(r["member_gene_ids"]),
            "category": None if pd.isna(r["category"]) else r["category"],
        }
        safety_layers = [
            c
            for c in hold_columns
            if c[5:]
            not in {
                "Reproductive / embryonic association",
                "Cancer expression",
                "Required RNA coverage",
                "Reproductive concentration",
            }
        ]
        safety_status = (
            "unsafe"
            if r["vital_protected_exclusion"]
            else "hold"
            if any(r[c] for c in safety_layers)
            else "unknown"
            if not r["required_coverage"]
            else "pass"
        )
        for axis, status, reason in [
            (
                "Safety",
                safety_status,
                "Absolute tissue/IHC RNA limits, qualifying normal protein and clinical safety evidence; required coverage is separate from absence of a detected exclusion.",
            ),
            (
                "Cancer expression",
                "pass" if r["cancer_association"] else "hold" if r["tumor_measured"] else "unknown",
                f"One sample ≥2: OncoRef {display(r['oncoref_maximum_tpm'])} TPM; HPA {display(r['hpa_tumor_group_lower_bound_ptpm'])} pTPM; reviewed reliable tumor protein: {r['reliable_tumor_protein']}. No tumor-versus-normal ratio.",
            ),
            (
                "Reproductive/developmental association",
                "pass" if r["positive_context"] else "hold",
                f"Reproductive RNA {display(r['reproductive_ntpm'])} nTPM; qualifying reproductive IHC {r['reproductive_protein_association']}; supported lineage-cell RNA {display(r['reproductive_cell_ntpm'])} nTPM; indirect cell support {r['reproductive_cell_indirect_association']} (outside clusters complete: {display(r['outside_cluster_complete'])}; outside maximum: {display(r['outside_cluster_max_ntpm'])} nTPM; normal induction counterevidence: {r['normal_induction_evidence']}); supported developmental RNA {display(r['developmental_mean_tpm'])} TPM. RNA routes require 5 in their own assay; supported cells and reliable reproductive IHC are independent presence routes. Normal and cancer gates still apply.",
            ),
        ]:
            row = dict.fromkeys(TABLE_COLUMNS["explanations"])
            row.update(
                base,
                layer="Assessment: " + axis,
                decision=status,
                reason=reason,
                source_url=BASE_URL,
            )
            rows.append(row)
        for column in hold_columns:
            layer = column[5:]
            hold, unsafe = r[column], r["unsafe:" + layer]
            if not hold:
                continue
            row = dict.fromkeys(TABLE_COLUMNS["explanations"])
            row.update(
                base,
                layer=layer.replace(
                    "Reproductive / embryonic association",
                    "Reproductive or developmental expression",
                ),
                decision="unsafe" if unsafe else "hold",
                source_url=BASE_URL,
            )
            candidate = peak_index.get((r["representative_gene_id"], layer, unsafe))
            if candidate is not None:
                p = pd.Series(candidate)
                row.update(
                    tissue=p.tissue,
                    cell_type=p.cell_type,
                    value=p.value,
                    unit=p.unit,
                    candidate_ceiling=p.candidate_ceiling,
                    baseline_ceiling=p.get("baseline_ceiling"),
                    clinical_product_ids=p.get("clinical_product_ids"),
                    clinical_credit_applied=p.get("clinical_credit_applied", False),
                    selection_tier=p.selection_tier,
                    ihc_coverage=f"{display(p.get('ihc_match_scope'))} / {p.negative_ihc_grade}",
                )
                row["reason"] = (
                    f"{p.tissue} / {p.cell_type}: {display(p.value)} {p.unit}; {p.selection_tier} ceiling {display(p.candidate_ceiling)}; {row['ihc_coverage']}."
                )
                if p.get("clinical_credit_applied", False):
                    row["reason"] += (
                        f" Baseline ceiling {display(p.baseline_ceiling)}; one-step clinical RNA credit"
                        f" from {p.clinical_product_ids}. Tested epitope/HLA/product limits apply."
                    )
                if pd.notna(p.get("Read count")):
                    row["reason"] += f" {int(p['Read count'])} reads."
                if pd.notna(p.get("value_upper")):
                    row["reason"] += (
                        f" Unpaired upper bound {display(p.value_upper)}; not a paired measurement."
                    )
                row.update(observation_source(p))
                if row["source_context"]:
                    row["reason"] += f" Native context: {row['source_context']}."
            elif layer == "HPA normal protein":
                candidate = protein_index.get((r["representative_gene_id"], unsafe))
                if candidate is not None:
                    p = pd.Series(candidate)
                    row.update(
                        tissue=p.tissue,
                        cell_type=p.cell_type,
                        value=p.Level,
                        unit="HPA staining category",
                        selection_tier=p.selection_tier,
                        reason=f"{p.Level} HPA staining in {p.tissue} / {p.cell_type}; {p.Reliability} reliability; {p.selection_tier} compartment.",
                        source_url=f"https://v23.proteinatlas.org/{p.gene_id}/tissue",
                    )
            elif layer == "Accepted normal-protein evidence":
                p = imported[
                    imported["Gene ID"].isin(base["gene_ids"]) & imported.retained_protein_exclusion
                ].iloc[0]
                row.update(
                    tissue=p.Context,
                    value=p.protein_level,
                    unit=p.protein_scale,
                    selection_tier=p["Tissue category"],
                    reason=f"Accepted protein evidence: {p.Context}; {p.protein_level} ({p.protein_scale}).",
                    source_url=p.Source,
                )
            elif layer == "Reproductive / embryonic association":
                row["reason"] = (
                    f"Reproductive {display(r['reproductive_ntpm'])} nTPM; supported developmental mean {display(r['developmental_mean_tpm'])} TPM; supported cell RNA {display(r['reproductive_cell_ntpm'])} nTPM; indirect cell support {r['reproductive_cell_indirect_association']}. No RNA/protein/cell route qualifies; normal induction counterevidence can invalidate an indirect-only route."
                )
                if r["normal_induction_evidence"]:
                    row["source_url"] = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/"
                    row["reason"] += (
                        " Published human plasmacytoid dendritic-cell induction contradicts indirect-only reproductive localization."
                    )
            elif layer == "Reproductive concentration":
                row["reason"] = (
                    f"HPA diffuseness {100 * r['hpa_diffuseness']:.2f}% exceeds 10%; reference {r['hpa_diffuseness_reference']:.3g} nTPM. Complete 29-unit adult/placental comparison, excluding thymus; no other qualifying Reproductive route. This is an identity failure, not a safety exclusion."
                )
                row["source_url"] = (
                    "https://pirl-unc.github.io/oncoref/audits/cta-unassigned-20261029/diffuseness-impact/"
                )
            elif layer == "Cancer expression":
                row["reason"] = (
                    f"OncoRef tumor maximum {display(r['oncoref_maximum_tpm'])} TPM; HPA lower bound {display(r['hpa_tumor_group_lower_bound_ptpm'])} pTPM. No qualifying tumor RNA or reviewed protein evidence."
                )
            elif layer == "Required RNA coverage":
                details = []
                combined = combined_bulk.loc[base["entity_id"]]
                if not combined.bulk_complete:
                    details.append(
                        f"Combined normal bulk missing required tissues: {combined.missing_bulk_tissues}"
                    )
                if not combined.gtex_core_complete:
                    details.append(
                        "Required GTEx v11 core panel is incomplete (52 native bulk contexts); "
                        f"missing member identities {combined.missing_gtex_member_ids}; "
                        f"uncovered tissue concepts {combined.missing_gtex_tissues}. "
                        "A missing measurement is not zero expression."
                    )
                for gap in coverage_by_entity.get(base["entity_id"], pd.DataFrame()).itertuples():
                    if gap.source_id == "hpa_rna_consensus":
                        continue
                    contexts = json.loads(gap.missing_contexts)
                    preview = ", ".join(contexts[:6])
                    if len(contexts) > 6:
                        preview += f", plus {len(contexts) - 6} more"
                    details.append(
                        f"{gap.source}: least-covered member has {gap.minimum_observed_contexts}/{gap.expected_contexts_per_member} contexts; incomplete members {gap.incomplete_member_gene_ids}; missing {preview}"
                    )
                row["reason"] = (
                    " | ".join(details)
                    or "Group-member coverage remains unresolved; inspect the source coverage audit."
                )
                row["source_url"] = BASE_URL + "coverage/"
            if not row["reason"]:
                row["reason"] = HUMAN_RULES.get(row["layer"], (layer, []))[0]
            rows.append(row)
        row = dict.fromkeys(TABLE_COLUMNS["explanations"])
        row.update(
            base,
            layer="Final assignment",
            decision="assignment",
            source_url=BASE_URL,
            reason=(
                f"Normal screen passes; no candidate subtype qualifies. Category subtype: {r['subtype']}. Categorical assignment: {r['classification']}. Cancer or reproductive evidence can be missing without implying missing normal coverage."
                if r["category"] == "Other"
                else f"{r['assessment_status']}: "
                + (
                    str(r["category"])
                    if pd.notna(r["category"])
                    else "no biological category assigned"
                )
            ),
        )
        rows.append(row)
        row = dict.fromkeys(TABLE_COLUMNS["explanations"])
        row.update(
            base,
            layer="Assessment: Tumor-supported",
            decision="lead" if r["tumor_restricted_profile"] else "not_qualified",
            reason=f"Distinct ≥10 TPM tumor samples (lower bound): {display(r['tumor_samples_ge10'])}; resolved tumor identity: {r['tumor_identity_resolved']}; reviewed tumor protein: {r['reliable_tumor_protein']}; low measured normal/developmental profile: {r['tumor_restricted_normal_low']}; normal-coverage review accepted: {r['normal_restriction_reviewed']}; lead awaiting review: {r['tumor_restricted_lead']}. All normal safety gates remain independent.",
            source_url=BASE_URL,
        )
        rows.append(row)
        row = dict.fromkeys(TABLE_COLUMNS["explanations"])
        row.update(
            base,
            layer="Assessment: Reproductive concentration",
            decision="hold"
            if r["reproductive_concentration_hold"]
            else "pass"
            if r["adult_profile_status"] == "concentrated"
            else "unknown",
            reason=f"Adult/placental contrast: {r['adult_profile_status']}; HPA score {display(100 * r['hpa_diffuseness'])}%; reference {display(r['hpa_diffuseness_reference'])} nTPM. Developmental contrast: {r['developmental_profile_status']}. Unassessed is not a measured pass.",
            source_url="https://pirl-unc.github.io/oncoref/audits/cta-unassigned-20261029/diffuseness-impact/",
        )
        rows.append(row)
        for flag in [
            "Antibody concern",
            "external_review_flag",
            "low_normal_protein_flag",
            "limited_protein_review_flag",
            "stem_progenitor_expression_flag",
        ]:
            if r[flag]:
                row = dict.fromkeys(TABLE_COLUMNS["explanations"])
                row.update(
                    base,
                    layer="Retained evidence flag",
                    decision="flag",
                    value=flag,
                    reason=flag.replace("_", " "),
                    source_url=BASE_URL,
                )
                rows.append(row)
    # Preserve low-valued comparator observations and matching provenance too.
    contexts = pd.read_csv(RELEASE / "widely-studied-contexts.csv.gz", low_memory=False)
    contexts = contexts[contexts.value.notna() & contexts.candidate_ceiling.notna()]
    wanted = contexts[
        contexts.get("ihc_match_scope").eq("reviewed_compartment") & contexts.value.ge(5)
    ]
    contexts = pd.concat(
        [
            contexts.sort_values("ratio", ascending=False).drop_duplicates(["gene_id", "assay"]),
            wanted,
        ]
    ).drop_duplicates(["gene_id", "assay", "tissue", "cell_type", "value"])
    lookup = facts.set_index("representative_gene_id")
    for p in contexts.to_dict("records"):
        r = lookup.loc[p["gene_id"]]
        row = dict.fromkeys(TABLE_COLUMNS["explanations"])
        row.update(
            entity_id=r.protein_equivalence_group_id,
            label=r.protein_group_label,
            gene_ids=json.loads(r.member_gene_ids),
            category=None if pd.isna(r.category) else r.category,
            layer=p["assay"],
            decision="observation",
            tissue=p["tissue"],
            cell_type=p["cell_type"],
            value=p["value"],
            unit=p["unit"],
            candidate_ceiling=p["candidate_ceiling"],
            baseline_ceiling=p.get("baseline_ceiling"),
            clinical_product_ids=p.get("clinical_product_ids"),
            clinical_credit_applied=p.get("clinical_credit_applied", False),
            selection_tier=p["selection_tier"],
            ihc_coverage=f"{display(p.get('ihc_match_scope'))}; {p['negative_ihc_grade']}; IHC cells: {display(p.get('ihc_cell_type'))}",
            reason=f"{p['tissue']} / {p['cell_type']}: {display(p['value'])} {p['unit']}; {p['selection_tier']} ceiling {display(p['candidate_ceiling'])}.",
            **observation_source(p),
        )
        if row["source_context"]:
            row["reason"] += f" Native context: {row['source_context']}."
        if p.get("clinical_credit_applied", False):
            row["reason"] += (
                f" Baseline ceiling {display(p['baseline_ceiling'])}; one-step clinical RNA credit"
                f" from {p['clinical_product_ids']}. Tested epitope/HLA/product limits apply."
            )
        rows.append(row)
    return pd.DataFrame(rows, columns=TABLE_COLUMNS["explanations"])


def native_assessment(species: str) -> dict:
    """Frozen audit of native evidence availability; no admission inference."""
    taxon = cta_snapshot(species).as_dict()["taxon_id"]
    data = {}
    candidates = cta_snapshot(species).candidates()
    stages = [
        (
            "Curated native nominations",
            "available" if candidates else "unavailable",
            f"{len(candidates)} native nominations; not a genome-wide screen",
        ),
        (
            "Reference and assay identity",
            "partial",
            "13 unique canine mappings and one ambiguous MAGE-family assay",
        ),
        (
            "Normal and reproductive/developmental RNA",
            "not_evaluated",
            "Bgee TPM and BarkBase TMM-CPM require native reference reconciliation",
        ),
        (
            "Tumor RNA",
            "not_evaluated",
            "GSE304066 count/FPKM and source-reference conflicts unresolved",
        ),
        (
            "Normal protein",
            "not_evaluated",
            "Study IHC is limited; Pfizer pooled iBAQ is not cell-resolved IHC",
        ),
        (
            "Quantitative selection",
            "not_evaluated",
            "No native species selection policy; no biological category assignment",
        ),
    ]
    if species == "mouse":
        stages = [
            *stages[:1],
            (
                "Reference and assay identity",
                "partial",
                "Pinned single-gene RefSeq subset; historical BALB/c assays do not establish C57BL/6J expression",
            ),
            (
                "Normal and reproductive/developmental RNA",
                "partial",
                "Five BALB/c normal organs with qualitative P1A RNA; somatic positives retained, donor counts unknown",
            ),
            (
                "Tumor RNA",
                "partial",
                "J558 cell-line RNA only; no tumor cohort or absolute abundance",
            ),
            ("Normal protein", "unavailable", "No protein observations imported in this seed"),
            (
                "Quantitative selection",
                "not_evaluated",
                "No quantitative mouse selection policy; no biological category assignment",
            ),
        ]
    data["funnel"] = [
        {
            "stage": stage,
            "status": status,
            "candidate_removals": None,
            "eligible": None,
            "Unsafe": None,
            "Other": None,
            "Unassigned": len(candidates),
            "note": note,
        }
        for stage, status, note in stages
    ]
    data["explanations"] = [
        dict.fromkeys(TABLE_COLUMNS["explanations"])
        | {
            "entity_id": f"nomination:{taxon}:{c['id']}",
            "label": c["label"],
            "gene_ids": c["gene_ids"],
            "category": None,
            "layer": "Selection coverage",
            "decision": "unassessed",
            "reason": "; ".join(c["evidence_coverage"]["missing"]),
            "source_url": "https://github.com/pirl-unc/oncoref/issues/574",
        }
        for c in candidates
    ]
    return data


def build() -> None:
    """Export the same accounting schema for measured and unassessed species."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    from oncoref.cta_selection_policy import POLICY_ID

    explanations = human_explanations(facts)
    files = {}
    for species, taxon in [("human", 9606), ("dog", 9615), ("mouse", 10090)]:
        data = dict(
            taxon_id=taxon,
            status="screened" if species == "human" else "unassessed",
            policy_id=POLICY_ID if species == "human" else None,
            **membership_tables(species),
        )
        if species == "human":
            funnel = pd.read_csv(RELEASE / "filtering-funnel.csv")
            data["funnel"] = [
                {
                    "stage": r.layer.replace(
                        "Reproductive / embryonic association",
                        "Reproductive or developmental expression",
                    ),
                    "status": "evaluated",
                    "candidate_removals": r.new_candidate_removals,
                    "eligible": r.CTA_eligible_so_far,
                    "Unsafe": r.Unsafe,
                    "Other": r.Other,
                    "Unassigned": r.Unassigned,
                    "note": "CTA eligibility so far; the final row assigns CTA, Unsafe and mutually exclusive Other categories",
                }
                for r in funnel.itertuples()
            ]
            data["explanations"] = records(explanations)
        else:
            data.update(native_assessment(species))
        data["criteria"] = [
            {
                "stage": r["stage"],
                "status": r["status"],
                "rule": HUMAN_RULES[r["stage"]][0] if species == "human" else r["note"],
                "source_ids": HUMAN_RULES[r["stage"]][1] if species == "human" else [],
            }
            for r in data["funnel"]
        ]
        for key, columns in TABLE_COLUMNS.items():
            pd.DataFrame(data[key], columns=columns).to_csv(
                OUTPUT / f"{species}-{key}.csv.gz",
                index=False,
                compression={"method": "gzip", "mtime": 0},
            )
        name = f"cta-selection-audit-{species}.json.gz"
        content = gzip.compress(
            json.dumps(data, allow_nan=False, separators=(",", ":")).encode(), mtime=0
        )
        (DATA / name).write_bytes(content)
        files[name] = hashlib.sha256(content).hexdigest()
    receipt = {
        "schema": "oncoref.cta-selection-audit",
        "version": 1,
        "decision_policy": POLICY_ID,
        "columns": TABLE_COLUMNS,
        "files": files,
        "inputs": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                RELEASE / "manifest.json",
                Path(__file__),
                DATA / "cta-publication-sources.csv",
                DATA / "cta-publication-membership.csv",
                RELEASE / "alternative-evidence/manifest.json",
                DATA / "cta-publication-context.csv",
                OUTPUT / "coverage-audit-receipt.json",
                *sorted((DATA / "species" / "dog").glob("cta-*.json")),
                *sorted((DATA / "species" / "mouse").glob("cta-*.json")),
            ]
        },
    }
    (DATA / "cta-selection-audit.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Exported {len(explanations):,} human explanations and all species tables")


if __name__ == "__main__":
    build()
