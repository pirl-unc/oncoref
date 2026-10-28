#!/usr/bin/env python3
"""Build shared species audit tables from the frozen 1.21.0 selection evidence.

This is an explanation/export pass: it does not change eligibility thresholds,
reclassify genes, substitute antibody matches, or infer canine measurements.
"""

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

from oncoref import cta_selection, hpa, normal_cell_screen  # noqa: E402
from oncoref.species_antigens import cta_snapshot  # noqa: E402

RELEASE = ROOT / "docs/audits/cta-release-1.21.0"
OUTPUT = ROOT / "docs/audits/cta-source-walkthrough-20261028"
DATA = ROOT / "oncoref/data"
BASE_URL = "https://pirl-unc.github.io/oncoref/audits/cta-release-1.21.0/"
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
        "Neither",
        "unique_entities",
        "status",
    ],
    "memberships": ["source_id", "entity_id", "label", "gene_ids", "category"],
    "overlaps": ["source_a", "source_b", "mapped_overlap", "selected_overlap"],
    "funnel": ["stage", "status", "candidate_removals", "eligible", "Unsafe", "Neither", "note"],
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
        "unsafe_threshold",
        "ihc_coverage",
        "reason",
        "source_url",
    ],
}
HUMAN_RULES = {
    "Reproductive or developmental expression": (
        "Reproductive ≥5 nTPM and > adult bulk maximum OR Tyser embryonic lineage mean ≥5 TPM; ≥20 cells, ≥5 positive and ≥10% positive at ≥1 TPM.",
        ["hpa_rna_consensus", "tyser_E-MTAB-9388"],
    ),
    "Tumor sample ≥2 TPM/pTPM": (
        "At least one measured compatible tumor sample ≥2 in native TPM or HPA pTPM; do not add maxima from different samples.",
        ["oncoref_tumor_union", "hpa_cancer"],
    ),
    "HPA normal protein": (
        "Medium/High or other non-Low positives exclude; accepted reliability supplies Unsafe evidence; uncertain positives retain a hold. Low remains flagged.",
        ["hpa_normal_tissue"],
    ),
    "Required RNA coverage": (
        "Required bulk and tissue-cell/immune panels; missing members are not zeros. DICE and cornea are optional coverage.",
        ["hpa_rna_consensus", "hpa_single_cell_tissue", "hpa_immune_cell", "hpa_immune_monaco"],
    ),
    "Normal bulk": (
        "Native nTPM at or below tissue/IHC ceiling passes; >ceiling holds; ≥2× supplies Unsafe evidence.",
        ["hpa_rna_consensus"],
    ),
    "HPA tissue clusters": (
        "Require ≥20 reads; apply exact cell/tissue/IHC ceiling. Same-tissue negative IHC alone does not establish cell-type coverage.",
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
    members, info = [], {}
    if species == "human":
        registry = pd.read_csv(DATA / "cta-publication-sources.csv").set_index("source_tag")
        known = pd.read_csv(RELEASE / "known-cta-decisions.csv.gz", low_memory=False)
        for doi, rows in registry.groupby("doi"):
            info[doi] = (citation(rows.index[0]), "https://doi.org/" + doi)
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
        for row in cta_snapshot(species).candidates():
            for source in row["derivation"]:
                members.append(
                    {
                        "source_id": source["source_id"],
                        "entity_id": f"nomination:9615:{row['id']}",
                        "label": row["label"],
                        "gene_ids": row["gene_ids"],
                        "category": "Neither",
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
            else {}
        )
        kind, status = "native_nomination", "unassessed"
    all_sets = {sid: set(grouped.loc[grouped.source_id.eq(sid), "entity_id"]) for sid in info}
    selected = {
        sid: set(grouped.loc[grouped.source_id.eq(sid) & grouped.category.eq("CTA"), "entity_id"])
        for sid in info
    }
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
                "unique_entities": len(all_sets[sid] - others),
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


def rna_peaks() -> pd.DataFrame:
    """Keep one supported worst-margin context per entity/assay for explanations.

    Complete context exports remain linked. Bounds for unpaired multi-locus
    DICE/cornea observations stay explicit; upper bounds are not measurements.
    """
    peaks = []
    keep = {
        "gene_id",
        "tissue",
        "cell_type",
        "value",
        "value_upper",
        "upper_supported_value",
        "unit",
        "negative_ihc_grade",
        "candidate_ceiling",
        "strong_boundary",
        "Read count",
        "candidate_exclusion",
        "unsafe_evidence",
        "assay",
        "category",
        "selection_tier",
    }
    for filename, assay in [
        ("bulk-contexts.csv.gz", "Normal bulk"),
        ("cluster-peak-evidence.csv.gz", "HPA tissue clusters"),
        ("immune-cornea-contexts.csv.gz", None),
        ("dice-group-contexts.csv.gz", "DICE immune"),
        ("cornea-group-contexts.csv.gz", "Adult cornea"),
    ]:
        for frame in pd.read_csv(RELEASE / filename, usecols=lambda x: x in keep, chunksize=150000):
            if assay:
                frame["assay"] = assay
            frame["ratio"] = frame.value / frame.candidate_ceiling
            for bound in ("value_upper", "upper_supported_value"):
                if bound in frame:
                    frame["ratio"] = pd.concat(
                        [frame.ratio, frame[bound] / frame.candidate_ceiling], axis=1
                    ).max(axis=1)
            frame["context_file"] = filename
            peaks.append(
                frame.sort_values(
                    ["unsafe_evidence", "candidate_exclusion", "ratio"], ascending=False
                ).drop_duplicates(["gene_id", "assay"])
            )
    return (
        pd.concat(peaks, ignore_index=True)
        .sort_values(["unsafe_evidence", "candidate_exclusion", "ratio"], ascending=False)
        .drop_duplicates(["gene_id", "assay"])
    )


def human_explanations(facts: pd.DataFrame, peaks: pd.DataFrame) -> pd.DataFrame:
    """Attach actual observations to each active layer, without inventing evidence."""
    index = peaks.set_index(["gene_id", "assay"])
    raw = hpa._read_hpa("hpa_normal_tissue", "v23").copy()
    raw = raw.rename(columns={"Gene": "gene_id", "Tissue": "tissue", "Cell type": "cell_type"})
    raw["source"], raw["version"], raw["taxonomy_id"] = "hpa_normal_tissue", "v23", 9606
    normal = normal_cell_screen.contexts(raw)
    normal = normal[normal.category.isin(["vital", "protected", "manageable"])]
    protein = normal[normal.Level.isin(["Medium", "High", "Ascending", "Descending"])].copy()
    protein["reliable"] = protein.Reliability.isin(["Approved", "Supported", "Enhanced"])
    protein["level_rank"] = protein.Level.map(
        {"Medium": 1, "High": 2, "Ascending": 3, "Descending": 3}
    )
    protein = (
        protein.sort_values(["reliable", "level_rank"], ascending=False)
        .drop_duplicates("gene_id")
        .set_index("gene_id")
    )
    negatives = raw[
        raw.Level.eq("Not detected") & raw.Reliability.isin(["Approved", "Supported", "Enhanced"])
    ].copy()
    negatives["description"] = negatives.cell_type + " (" + negatives.Reliability + ")"
    same_tissue = negatives.groupby(["gene_id", "tissue"]).description.agg("; ".join).to_dict()
    imported = pd.read_csv(RELEASE / "imported-protein-review-decisions.csv").set_index("Gene ID")
    findings = pd.read_csv(DATA / "cta-external-safety-reviews.csv")
    gaps = pd.read_csv(
        ROOT / "docs/audits/cta-expanded-funnel-20261027/mandatory-coverage-gaps.csv.gz",
        low_memory=False,
    ).set_index("gene_id")
    compartments = pd.read_csv(RELEASE / "changed-compartment-contexts.csv.gz")
    rows = []
    for fact in facts.to_dict("records"):
        gene = fact["representative_gene_id"]
        base = {
            "entity_id": fact["protein_equivalence_group_id"],
            "label": fact["protein_group_label"],
            "gene_ids": json.loads(fact["member_gene_ids"]),
            "category": fact["category"],
        }

        def add(layer: str, decision: str, reason: str, *, base=base, **values: object) -> None:
            rows.append(
                {
                    **dict.fromkeys(TABLE_COLUMNS["explanations"]),
                    **base,
                    "layer": layer,
                    "decision": decision,
                    "reason": reason,
                    **values,
                }
            )

        repro, adult, embryo = (
            fact["reproductive_ntpm"],
            fact["adult_bulk_max_ntpm"],
            fact["developmental_mean_tpm"],
        )
        embryo_text = (
            f"{embryo:.3g} TPM" if pd.notna(embryo) else "no qualifying supported measurement"
        )
        add(
            "Reproductive or developmental expression",
            "pass" if fact["positive_context"] else "hold",
            f"Reproductive maximum {display(repro)} nTPM; requires ≥5 and > adult bulk maximum {display(adult)} nTPM. "
            f"Supported embryonic lineage: {embryo_text}; requires mean ≥5 TPM with cell-support guards.",
            source_url=BASE_URL,
        )
        add(
            "Tumor sample ≥2 TPM/pTPM",
            "pass" if fact["tumor_union_ge2"] else "hold",
            f"OncoRef maximum {display(fact['oncoref_maximum_tpm'])} TPM, sample {display(fact['oncoref_maximum_sample'])} "
            f"({display(fact['oncoref_cancer_code'])}); HPA group lower bound {display(fact['hpa_tumor_group_lower_bound_ptpm'])} pTPM. Requires ≥2 in either compatible source.",
            source_url="https://pirl-unc.github.io/oncoref/audits/cta-expanded-funnel-20261027/tumor-union/tumor-maxima.csv.gz",
        )
        for layer in [col[5:] for col in facts if col.startswith("hold:")]:
            if layer in ("Reproductive / embryonic association", "Tumor sample ≥2 TPM/pTPM"):
                continue
            excluded = fact["hold:" + layer]
            unsafe = fact["unsafe:" + layer]
            if not (excluded or unsafe) and (gene, layer) not in index.index:
                continue
            decision = "unsafe" if unsafe else "hold" if excluded else "pass"
            if (gene, layer) in index.index:
                obs = index.loc[(gene, layer)]
                unit, value, ceiling, threshold = (
                    obs.unit,
                    obs.value,
                    obs.candidate_ceiling,
                    obs.strong_boundary,
                )
                cell = obs.get("cell_type")
                cell = cell if pd.notna(cell) else "bulk tissue"
                ihc = str(obs.negative_ihc_grade)
                if ihc in ("Missing", "Missing/uncertain", "Uncertain"):
                    negative = same_tissue.get((gene, obs.tissue))
                    ihc = (
                        ("Same-tissue negative IHC; cell-type coverage unresolved: " + negative)
                        if negative
                        else "No usable matched negative IHC"
                    )
                detail = f"{obs.tissue} / {cell}: {display(value)} {unit}; candidate ceiling {display(ceiling)}; Unsafe threshold {display(threshold)}. {ihc}."
                if pd.isna(value):
                    detail += " No supported quantitative summary; absence of a triggered exclusion is not measured zero."
                if excluded and not obs.candidate_exclusion:
                    detail += " Layer hold includes unresolved coverage/anatomy beyond this representative measured context."
                if pd.notna(obs.get("Read count")):
                    detail += f" {int(obs['Read count'])} reads."
                for bound in ("value_upper", "upper_supported_value"):
                    if pd.notna(obs.get(bound)) and obs[bound] > value:
                        detail += f" Unpaired-group upper bound {obs[bound]:g} {unit}; not a paired measurement."
                add(
                    layer,
                    decision,
                    detail,
                    tissue=obs.tissue,
                    cell_type=cell,
                    value=value,
                    unit=unit,
                    candidate_ceiling=ceiling,
                    unsafe_threshold=threshold,
                    ihc_coverage=ihc,
                    source_url=BASE_URL + obs.context_file,
                )
            elif layer == "HPA normal protein":
                options = protein[protein.index.isin(base["gene_ids"])]
                if not options.empty:
                    obs = options.iloc[0]
                    add(
                        layer,
                        decision,
                        f"{obs.tissue} / {obs.cell_type}: {obs.Level} HPA staining, {obs.Reliability} reliability.",
                        tissue=obs.tissue,
                        cell_type=obs.cell_type,
                        value=obs.Level,
                        unit="HPA staining category",
                        source_url=f"https://v23.proteinatlas.org/{options.index[0]}/tissue",
                    )
                else:
                    add(
                        layer,
                        decision,
                        "Archived positive HPA hold; no matching representative context in this export; inspect full member evidence.",
                        source_url=BASE_URL + "hpa-protein-level-summary.csv.gz",
                    )
            elif layer == "Accepted normal-protein evidence":
                options = imported[imported.index.isin(base["gene_ids"])]
                for _, obs in options.iterrows():
                    add(
                        layer,
                        decision,
                        f"{obs['Context']}: {obs['protein_level']}. {obs['Rationale']}",
                        tissue=obs["Context"],
                        value=obs["protein_level"],
                        unit=obs["protein_scale"],
                        source_url=obs["Source"],
                    )
            elif layer == "Clinical safety exclusion":
                options = findings[findings.Ensembl_Gene_ID.isin(base["gene_ids"])]
                for _, obs in options.iterrows():
                    add(
                        layer,
                        decision,
                        f"{obs.rationale} {obs.limitations}",
                        source_url=obs.source_url,
                    )
                if options.empty:
                    add(
                        layer,
                        decision,
                        "Reviewed clinical exclusion; see product-specific finding and attribution limits.",
                        source_url=BASE_URL + "problematic-target-checks.csv",
                    )
            elif layer == "Required RNA coverage":
                details = []
                for member in base["gene_ids"]:
                    if member not in gaps.index:
                        continue
                    evidence = gaps.loc[member]
                    for field, required, panel in [
                        ("normal_bulk_tissues", 50, "bulk tissues"),
                        ("cell_records", 426, "tissue-cell records"),
                        ("hpa_immune_cell_records", 19, "HPA immune records"),
                        ("hpa_immune_monaco_records", 30, "Monaco records"),
                    ]:
                        if evidence[field] != required:
                            details.append(
                                f"{member}: {display(evidence[field])}/{required} {panel}"
                            )
                reason = (
                    "; ".join(details)
                    or "Required measurement or complete-member aggregation unavailable; see member coverage"
                )
                add(
                    layer,
                    decision,
                    reason + ". Optional DICE/cornea absence does not cause this hold.",
                    source_url="https://pirl-unc.github.io/oncoref/audits/cta-expanded-funnel-20261027/mandatory-coverage-gaps.csv.gz",
                )
            else:
                add(
                    layer,
                    decision,
                    "Unresolved source/anatomical or protein evidence; no quantitative value is imputed.",
                    source_url=BASE_URL + "protein-group-classifications.csv.gz",
                )
        if fact["low_normal_protein_flag"]:
            options = imported[
                imported.index.isin(base["gene_ids"]) & imported.low_normal_protein_flag
            ]
            for _, obs in options.iterrows():
                add(
                    "Normal protein observation",
                    "flag",
                    f"{obs['Context']}: {obs['protein_level']} ({obs['protein_scale']}); retained positive finding, no negative-IHC credit.",
                    tissue=obs["Context"],
                    value=obs["protein_level"],
                    unit=obs["protein_scale"],
                    source_url=obs["Source"],
                )
        if fact["protein_group_label"] in COMPARATORS:
            for obs in compartments[compartments.gene_id.eq(gene)].itertuples():
                ihc = same_tissue.get((gene, obs.tissue))
                ihc = (
                    "Same-tissue negative IHC; cell-type coverage unresolved: " + ihc
                    if ihc
                    else "No usable matched negative IHC"
                )
                add(
                    "Additional tissue-cluster context",
                    "observation",
                    f"{obs.tissue} / {obs.cell_type}: {display(obs.value)} {obs.unit}; candidate ceiling {display(obs.candidate_ceiling)}; Unsafe threshold {display(obs.strong_boundary)}. {ihc}.",
                    tissue=obs.tissue,
                    cell_type=obs.cell_type,
                    value=obs.value,
                    unit=obs.unit,
                    candidate_ceiling=obs.candidate_ceiling,
                    unsafe_threshold=obs.strong_boundary,
                    ihc_coverage=ihc,
                    source_url=BASE_URL + "changed-compartment-contexts.csv.gz",
                )
    return pd.DataFrame(rows, columns=TABLE_COLUMNS["explanations"])


def build() -> None:
    """Export source accounting and detailed decisions in one schema for all species."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    peaks = rna_peaks()
    peaks.to_csv(
        OUTPUT / "normal-rna-evidence.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    explanations = human_explanations(facts, peaks)
    explanations.to_csv(
        OUTPUT / "human-explanations.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    files = {}
    for species, taxon in [("human", 9606), ("dog", 9615), ("mouse", 10090)]:
        data = {
            "taxon_id": taxon,
            "status": "screened" if species == "human" else "unassessed",
            "policy_id": "oncoref.cta-selection.1.21.0" if species == "human" else None,
            **membership_tables(species),
        }
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
                    "Neither": r.Neither,
                    "note": "Eligible so far until the final stage; every layer evaluates every group",
                }
                for r in funnel.itertuples()
            ]
            data["explanations"] = records(explanations)
        else:
            candidates = cta_snapshot(species).candidates()
            n = len(candidates)
            data["funnel"] = [
                {
                    "stage": "Curated native nominations",
                    "status": "available" if n else "unavailable",
                    "candidate_removals": None,
                    "eligible": None,
                    "Unsafe": 0,
                    "Neither": n,
                    "note": f"{n} nominations; not a genome-wide screen",
                },
                {
                    "stage": "Reference and assay identity",
                    "status": "partial" if n else "unavailable",
                    "candidate_removals": None,
                    "eligible": None,
                    "Unsafe": 0,
                    "Neither": n,
                    "note": "Dog: 13 unique gene nominations, 1 ambiguous MAGE family assay; unique mapping is not measured protein equivalence"
                    if n
                    else "No integrated native nomination identities",
                },
                {
                    "stage": "Normal and reproductive/developmental RNA",
                    "status": "not_evaluated",
                    "candidate_removals": None,
                    "eligible": None,
                    "Unsafe": 0,
                    "Neither": n,
                    "note": "Dog native matrices retain TPM/TMM-CPM and unresolved quantification-reference equivalence; no selection counts",
                },
                {
                    "stage": "Tumor RNA",
                    "status": "not_evaluated",
                    "candidate_removals": None,
                    "eligible": None,
                    "Unsafe": 0,
                    "Neither": n,
                    "note": "Dog GSE304066 has count/FPKM and source-reference conflicts; study PCR is not a TPM cohort",
                },
                {
                    "stage": "Normal IHC and mass spectrometry",
                    "status": "not_evaluated",
                    "candidate_removals": None,
                    "eligible": None,
                    "Unsafe": 0,
                    "Neither": n,
                    "note": "Dog study IHC is limited; Pfizer iBAQ is pooled bulk protein, not cell IHC or DLA presentation",
                },
                {
                    "stage": "Quantitative expression and protein selection",
                    "status": "not_evaluated",
                    "candidate_removals": None,
                    "eligible": None,
                    "Unsafe": 0,
                    "Neither": n,
                    "note": "Native reference/units, required panels and selection policy unresolved",
                },
            ]
            data["explanations"] = [
                {
                    **dict.fromkeys(TABLE_COLUMNS["explanations"]),
                    "entity_id": f"nomination:{taxon}:{c['id']}",
                    "label": c["label"],
                    "gene_ids": c["gene_ids"],
                    "category": "Neither",
                    "layer": "Selection coverage",
                    "decision": "unassessed",
                    "reason": "; ".join(c["evidence_coverage"]["missing"]),
                    "source_url": "https://github.com/pirl-unc/oncoref/issues/574",
                }
                for c in candidates
            ]
            if species == "mouse":
                data["funnel"] = data["funnel"][:1]
        data["criteria"] = [
            {
                "stage": row["stage"],
                "status": row["status"],
                "rule": HUMAN_RULES[row["stage"]][0] if species == "human" else row["note"],
                "source_ids": HUMAN_RULES[row["stage"]][1] if species == "human" else [],
            }
            for row in data["funnel"]
        ]
        for key in TABLE_COLUMNS:
            frame = pd.DataFrame(data[key], columns=TABLE_COLUMNS[key])
            frame.to_csv(
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
        "decision_policy": "oncoref.cta-selection.1.21.0",
        "columns": TABLE_COLUMNS,
        "files": files,
        "native_table_receipts": {
            "hpa_normal_tissue_v23_csv_sha256": hashlib.sha256(
                hpa._read_hpa("hpa_normal_tissue", "v23").to_csv(index=False).encode()
            ).hexdigest()
        },
        "inputs": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                RELEASE / "manifest.json",
                RELEASE / "protein-group-classifications.csv.gz",
                RELEASE / "known-cta-decisions.csv.gz",
                RELEASE / "filtering-funnel.csv",
                RELEASE / "imported-protein-review-decisions.csv",
                RELEASE / "bulk-contexts.csv.gz",
                RELEASE / "cluster-peak-evidence.csv.gz",
                RELEASE / "immune-cornea-contexts.csv.gz",
                RELEASE / "dice-group-contexts.csv.gz",
                RELEASE / "cornea-group-contexts.csv.gz",
                RELEASE / "changed-compartment-contexts.csv.gz",
                ROOT / "docs/audits/cta-expanded-funnel-20261027/mandatory-coverage-gaps.csv.gz",
                DATA / "cta-publication-sources.csv",
                DATA / "cta-external-safety-reviews.csv",
                DATA / "normal-tissue-taxonomy.yaml",
                ROOT / "oncoref/normal_cell_screen.py",
                ROOT / "oncoref/species_antigens.py",
                *sorted((DATA / "species").glob("*/*cta*")),
                Path(__file__),
            ]
        },
    }
    (DATA / "cta-selection-audit.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {"explanation_rows": len(explanations), "output": str(OUTPUT), "artifacts": files},
            indent=2,
        )
    )


if __name__ == "__main__":
    build()
