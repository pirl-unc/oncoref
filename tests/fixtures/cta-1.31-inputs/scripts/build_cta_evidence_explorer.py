#!/usr/bin/env python3
"""Explain released CTA calls without changing the selection policy or data.

Normal context replay reuses the release's identity resolution and same-context
protein sums. Top-five observations are ranked within a source and native unit;
unsupported measurements and reproductive contexts remain distinguishable.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from oncoref import cta_selection, gene_knowledge, tissue_hazards
from oncoref import cta_selection_policy as policy

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-selection-20261028"
OUT = ROOT / "docs/audits/cta-evidence-guide-20261106"
ASSETS = ROOT / "scripts/cta_evidence_assets"
KEYS = ["gene_id", "source", "unit", "scope"]
CONTEXT = [*KEYS, "tissue", "cell_type"]
FIELDS = [
    *CONTEXT,
    "value",
    "raw_value",
    "selection_tier",
    "candidate_ceiling",
    "candidate_exclusion",
    "unsafe_evidence",
    "unresolved_context",
    "negative_ihc_grade",
    "source_context",
    "Read count",
    "value_upper",
    "bound_reason",
    "statistic",
    "version",
]
REVIEW_LABELS = {
    "incomplete_evidence": "Missing measurements",
    "normal_evidence_hold": "Conflicting or uncertain normal evidence",
    "translation_review": "Uncertain translation",
    "unassessed": "Not assessed",
}
ROLE_LABELS = {"allow": "No sequence exclusion", "review": "Unresolved", "exclude": "Exclude"}


def clean(value):
    """Serialize missing observations as null, never as measured zero."""
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if value is None or value is pd.NA:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not pd.notna(value):
        return None
    if isinstance(value, float) and abs(value) == float("inf"):
        return "unbounded"
    return value


def write_json(path, value, compressed=False):
    content = (
        json.dumps(clean(value), ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n"
    ).encode()
    path.write_bytes(gzip.compress(content, mtime=0) if compressed else content)


def write_csv(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if str(path).endswith(".gz") else None,
    )


def source_key(source, assay):
    if assay == "DICE immune":
        return "dice"
    if assay == "Adult cornea":
        return "GSE155683_adult"
    return source


def top_contexts(frame):
    """Highest supported values, not the largest limit ratios, per compartment."""
    d = frame.sort_values(
        ["value", "raw_value", "tissue", "cell_type"],
        ascending=[False, False, True, True],
        na_position="last",
    )
    return d.drop_duplicates(CONTEXT).groupby(KEYS, dropna=False, sort=False).head(5)


class ExplorerCollector:
    def __init__(self, facts):
        self.facts = facts.set_index("representative_gene_id")
        self.parts = []

    def collect(self, frame, assay):
        d = frame[frame.gene_id.isin(self.facts.index)].copy()
        for name in [
            "clinical_rna_credit",
            "clinical_product_ids",
            "clinical_cohort_ids",
            "clinical_credit_scope",
        ]:
            d[name] = d.gene_id.map(self.facts[name])
        d = policy.assess_contexts(d)
        d["source"] = d.source.map(lambda source: source_key(source, assay))
        d["raw_value"] = d.get("raw_maximum", d.value).fillna(d.value)
        d["scope"] = d.selection_tier.map(
            lambda tier: (
                "Reproductive/developmental"
                if tier in {"reproductive", "developmental"}
                else "Other normal"
            )
        )
        if "value_upper" in d:
            complete = d.measured_members.eq(d.member_count)
            d["candidate_exclusion"] |= complete & d.value_upper.gt(d.candidate_ceiling)
            d["unresolved_context"] |= (
                complete & d.value_upper.gt(d.candidate_ceiling) & ~d.value.gt(d.candidate_ceiling)
            )
        if "upper_supported_value" in d:
            possible = d.upper_supported_value.gt(d.candidate_ceiling)
            d["candidate_exclusion"] |= possible
            d["unresolved_context"] |= possible & ~d.value.gt(d.candidate_ceiling)
        d = d.reindex(columns=FIELDS)
        self.parts.append(top_contexts(d))

    def finish(self):
        d = top_contexts(pd.concat(self.parts, ignore_index=True))
        write_csv(d, OUT / "normal-top-contexts.csv.gz")
        return d


def replay(facts):
    from replay_cta_selection import replay_contexts

    # The release's replay helper also emits IHC and source receipts. Isolate
    # those writes; no frozen release artifact is refreshed by this report.
    with tempfile.TemporaryDirectory(prefix="oncoref-explorer-") as tmp:
        scratch = Path(tmp)
        for path in AUDIT.glob("*-accepted-identities.csv.gz"):
            shutil.copyfile(path, scratch / path.name)
        replay_contexts(facts, scratch, collector=ExplorerCollector(facts))


class BiologyCollector:
    """Retain all supported category observations, before the display's top-five cut."""

    def __init__(self, facts):
        self.genes = set(facts.representative_gene_id)
        self.parts = []

    def collect(self, frame, assay):
        if assay not in {"Normal bulk", "HPA tissue clusters"}:
            return
        selected = frame[
            frame.gene_id.isin(self.genes)
            & (frame.tissue.isin({"testis", "placenta"}) | frame.cell_type.eq("oocytes"))
            & frame.value.ge(5)
            & frame.unit.isin({"TPM", "nTPM"})
        ].copy()
        selected["scope"] = "Reproductive/developmental"
        self.parts.append(selected.reindex(columns=[*CONTEXT, "value", "Read count", "version"]))

    def finish(self):
        result = pd.concat(self.parts, ignore_index=True).drop_duplicates(CONTEXT)
        result = result.sort_values(CONTEXT)
        write_csv(result, OUT / "biology-rna-contexts.csv.gz")
        return result


def replay_biology(facts):
    from replay_cta_selection import replay_contexts

    ctas = set(cta_selection.results(category="CTA").entity_id)
    selected = facts[facts.protein_equivalence_group_id.isin(ctas)]
    with tempfile.TemporaryDirectory(prefix="oncoref-biology-") as tmp:
        scratch = Path(tmp)
        for path in AUDIT.glob("*-accepted-identities.csv.gz"):
            shutil.copyfile(path, scratch / path.name)
        replay_contexts(selected, scratch, collector=BiologyCollector(selected))
        shutil.copyfile(scratch / "cluster-source.json", OUT / "biology-rna-source.json")


def route(row):
    adult, dev = bool(row.reproductive_route), bool(row.developmental_route)
    return (
        "Adult reproductive + developmental"
        if adult and dev
        else "Developmental only"
        if dev
        else "Adult reproductive only"
        if adult
        else "No qualifying association"
    )


def review_categories(row):
    categories = []
    if row.assessment_status in REVIEW_LABELS:
        categories.append(REVIEW_LABELS[row.assessment_status])
    if row.identity_status != "resolved":
        categories.append("Uncertain identity")
    reasons = str(row.get("translation_qualification_reasons", ""))
    if "ambiguous_or_interrupted_sequence" in reasons or "noncanonical_start_residue" in reasons:
        categories.append("Ambiguous sequence")
    if row.shared_sequence == "review" and not categories:
        categories.append("Not assessed")
    return categories


def biological_categories(protein, fact, observations, developmental, registry):
    """Attach biological context, retaining the exact observation behind each tag."""
    if protein.category != "CTA":
        return [], []
    evidence = []
    for assignment in registry["curated_assignments"]:
        if set(protein.gene_ids) & set(assignment["gene_ids"]):
            evidence.append({**assignment, "kind": "Publication"})
    for observation in observations:
        value = observation["value"]
        if (
            observation["scope"] != "Reproductive/developmental"
            or observation["unit"] not in {"TPM", "nTPM"}
            or not isinstance(value, (int, float))
            or value < registry["rules"]["rna_minimum"]
        ):
            continue
        tissue, cell = observation["tissue"], observation["cell_type"]
        category = (
            "oocyte_maternal"
            if cell == "oocytes"
            else "testis"
            if tissue == "testis"
            else "placenta"
            if tissue == "placenta"
            else None
        )
        if category:
            evidence.append(
                {
                    "category_id": category,
                    "kind": "RNA measurement",
                    "source_id": observation["source"],
                    "gene_ids": list(protein.gene_ids),
                    "model": "Human reproductive tissue",
                    "basis": f"{tissue} / {cell}: {value:g} {observation['unit']}",
                }
            )
    for observation in developmental:
        if observation["Mean TPM"] < registry["rules"]["rna_minimum"]:
            continue
        category = {
            "petropoulos_E-MTAB-3929": "early_embryo",
            "tyser_E-MTAB-9388": "gastrulation",
        }.get(observation["source"])
        if category:
            evidence.append(
                {
                    "category_id": category,
                    "kind": "RNA measurement",
                    "source_id": observation["source"],
                    "gene_ids": [observation["Gene ID"]],
                    "model": "Human " + observation["stage"],
                    "basis": f"{observation['Lineage']}: {observation['Mean TPM']:g} TPM",
                }
            )
    assigned = {e["category_id"] for e in evidence}
    for needed, category, specific in [
        (
            fact.reproductive_route,
            "other_reproductive",
            {"oocyte_maternal", "reproductive_secretion", "testis", "placenta"},
        ),
        (fact.developmental_route, "other_developmental", {"early_embryo", "gastrulation"}),
    ]:
        if needed and not assigned.intersection(specific):
            assigned.add(category)
            evidence.append(
                {
                    "category_id": category,
                    "kind": "Released association",
                    "source_id": "current_decisions",
                    "gene_ids": list(protein.gene_ids),
                    "model": "Released admission record",
                    "basis": "Reproductive association"
                    if category == "other_reproductive"
                    else "Developmental association",
                }
            )
    if not assigned:
        assigned.add("uncategorized")
    return [c["id"] for c in registry["categories"] if c["id"] in assigned], evidence


def source_guide(audit):
    sources = cta_selection.sources()
    sources = sources[sources.taxon_id.eq(9606)]
    items = []
    for _sid, part in sources.groupby("source_id", sort=False):
        row = part.iloc[0].to_dict()
        row["evidence_category"] = ", ".join(part.evidence_category.unique())
        items.append(row)
    byid = {s["source_id"]: s for s in items}
    for sid, name, unit, url, limits in [
        (
            "CELLxGENE",
            "CELLxGENE rare-cell screen",
            "% cells with detected RNA",
            "https://cellxgene.cziscience.com/",
            "Existing brain/spinal-cord, heart and lung panel only; >=3,000 cells and >=10 positive cells. >=1% detection is a policy trigger, not an abundance measurement or proof of receptor recognition.",
        ),
        (
            "Ensembl-115",
            "Ensembl 115 coding translations",
            "annotation / amino-acid sequence",
            "https://sep2025.archive.ensembl.org/Homo_sapiens/Info/Index",
            "Qualifying annotated CDS is accepted translation support; it does not prove protein detection in a tumor. Exact selected sequence and every member locus must qualify.",
        ),
        (
            "unsafe_sequences",
            "Sequence overlap with released Unsafe proteins",
            "% shared 9-mers",
            "../sequence-bundle-20261104/",
            "Computed from Ensembl 115 coding translations and current Unsafe assignments. >=95% containment blocks CTA admission; it does not itself establish tissue expression.",
        ),
        (
            "current_decisions",
            "All released decision gates",
            "classification / reasons",
            "../cta-definition-20261104/report.html",
            "Derived result, not an independent experimental source. Unsafe takes precedence over CTA admission; Other retains its evidence status.",
        ),
    ]:
        row = {
            "source_id": sid,
            "name": name,
            "native_unit": unit,
            "source_url": url,
            "limitations": limits,
            "stage": "source-specific",
            "evidence_category": "computed screen"
            if sid != "Ensembl-115"
            else "translation annotation",
            "status": "used_in_selection",
        }
        byid[sid] = row
    for r in audit.sources.to_dict("records"):
        byid.setdefault(
            r["source_id"],
            {
                "source_id": r["source_id"],
                "name": r["citation"],
                "native_unit": "publication membership",
                "source_url": r["source_url"],
                "limitations": r.get("evidence_scope", "Nomination is not a safety assessment."),
                "stage": "publication-specific",
                "evidence_category": "nomination",
                "status": "annotation",
            },
        )
    return list(byid.values())


def build(replay_contexts=False):
    OUT.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(AUDIT / "protein-group-classifications.csv.gz", low_memory=False)
    if replay_contexts or not (OUT / "normal-top-contexts.csv.gz").exists():
        replay(facts)
    if replay_contexts or not (OUT / "biology-rna-contexts.csv.gz").exists():
        replay_biology(facts)
    profiles = cta_selection.profiles()
    roles = cta_selection.sequence_dispositions().set_index("entity_id")
    audit = cta_selection.audit()
    members = (
        audit.memberships.groupby("entity_id").source_id.agg(lambda x: sorted(set(x))).to_dict()
    )
    facts = facts.set_index("protein_equivalence_group_id")
    normal = pd.read_csv(OUT / "normal-top-contexts.csv.gz", low_memory=False)
    top = {
        k: clean(v.drop(columns="gene_id").to_dict("records"))
        for k, v in normal.groupby("gene_id", sort=False)
    }
    biology_rna = {
        gene: clean(part.drop(columns="gene_id").to_dict("records"))
        for gene, part in pd.read_csv(OUT / "biology-rna-contexts.csv.gz").groupby("gene_id")
    }
    knowledge = gene_knowledge.findings()
    knowledge_by_gene = {}
    for finding in knowledge.to_dict("records"):
        for gene in finding["gene_ids"]:
            knowledge_by_gene.setdefault(gene, {})[finding["id"]] = finding
    ledgers = {}
    for sid, path, gene_column in [
        ("cta_policy_findings", ROOT / "oncoref/data/cta-policy-findings.csv", "Gene ID"),
        (
            "normal-induction-evidence",
            ROOT / "oncoref/data/normal-induction-evidence.csv",
            "gene_id",
        ),
        ("external_ihc", AUDIT / "external-protein-evidence.csv.gz", "gene_id"),
    ]:
        ledger = pd.read_csv(path, keep_default_na=False)
        ledgers[sid] = {g: part.to_dict("records") for g, part in ledger.groupby(gene_column)}
    literature = json.loads((ASSETS / "literature-review.json").read_text())
    biology = json.loads((ASSETS / "biology-categories.json").read_text())
    lit_by_label = {}
    for item in literature["findings"]:
        for label in item["proteins"]:
            lit_by_label.setdefault(label, []).append(item)
    # Native developmental lineage means remain separated by study and stage.
    dev = pd.read_csv(AUDIT / "developmental-lineages.csv.gz", low_memory=False)
    dev["supported"] = (
        dev.Cells.ge(20) & dev["Positive cells TPM>=1"].ge(5) & dev["Fraction TPM>=1"].ge(0.1)
    )
    dev = dev[dev.supported].sort_values("Mean TPM", ascending=False)
    devtop = {
        g: d.groupby("source", sort=False).head(5).to_dict("records")
        for g, d in dev.groupby("Gene ID", sort=False)
    }
    ihc = pd.read_csv(AUDIT / "hpa-protein-evidence.csv.gz", low_memory=False)
    ihc["rank"] = ihc.Level.map({"High": 3, "Medium": 2, "Low": 1, "Not detected": 0}).fillna(-1)
    ihc = ihc.sort_values(["rank", "tissue", "cell_type"], ascending=[False, True, True])
    ihctop = {
        g: d.head(5).drop(columns="rank").to_dict("records")
        for g, d in ihc.groupby("gene_id", sort=False)
    }
    sources = source_guide(audit)
    normal_ids = [
        "hpa_rna_consensus",
        "hpa_brain",
        "gtex_tissue_rna",
        "hpa_single_cell_tissue",
        "hpa_immune_cell",
        "hpa_immune_monaco",
        "dice",
        "GSE155683_adult",
    ]
    numeric_columns = [
        *normal_ids,
        "oncoref_tumor_union",
        "hpa_cancer",
        "tyser_E-MTAB-9388",
        "petropoulos_E-MTAB-3929",
        "CELLxGENE",
        "unsafe_sequences",
    ]
    columns = [
        *numeric_columns,
        "hpa_normal_tissue",
        "Ensembl-115",
        "gene-knowledge",
        "cta_policy_findings",
        "external_ihc",
        "normal-induction-evidence",
        *audit.sources.source_id,
    ]
    explanations = {
        k: d.to_dict("records") for k, d in audit.explanations.groupby("entity_id", sort=False)
    }
    rows, subset, details = [], [], {}
    for p in profiles.itertuples(index=False):
        f = facts.loc[p.entity_id]
        role = roles.loc[p.entity_id]
        current = {key: clean(getattr(p, key)) for key in profiles.columns}
        role_for_review = role.copy()
        role_for_review["translation_qualification_reasons"] = p.translation_qualification_reasons
        current.update(
            papers=members.get(p.entity_id, []),
            route=route(f),
            sequence_role=role.shared_sequence,
            sequence_label=ROLE_LABELS[role.shared_sequence],
            review_categories=review_categories(role_for_review),
        )
        obs = top.get(f.representative_gene_id, [])
        values = dict.fromkeys(columns)
        for source in normal_ids:
            measured = [
                r["value"]
                for r in obs
                if r["source"] == source and isinstance(r["value"], (int, float))
            ]
            values[source] = max(measured) if measured else None
        values.update(
            oncoref_tumor_union=clean(p.tumor_max_tpm),
            hpa_cancer=clean(p.tumor_max_ptpm),
            CELLxGENE=clean(p.rare_cell_max_pct) if p.rare_cell_screened else None,
            unsafe_sequences=clean(p.unsafe_sequence_overlap_pct),
        )
        # Do not sum maxima across disjoint developmental contexts for shared
        # proteins. Group qualification is reported separately in release facts.
        drows = []
        for gene in p.gene_ids:
            drows.extend(devtop.get(gene, []))
        current["biology_categories"], current["biology_evidence"] = biological_categories(
            p, f, biology_rna.get(f.representative_gene_id, []), drows, biology
        )
        for source in ["tyser_E-MTAB-9388", "petropoulos_E-MTAB-3929"]:
            measured = [r["Mean TPM"] for r in drows if r["source"] == source]
            values[source] = max(measured) if measured and len(p.gene_ids) == 1 else None
        prows = [r for g in p.gene_ids for r in ihctop.get(g, [])]
        ranks = {"High": 3, "Medium": 2, "Low": 1, "Not detected": 0}
        levels = [r["Level"] for r in prows if r["Level"] in ranks]
        values["hpa_normal_tissue"] = max(levels, key=ranks.get) if levels else None
        values["Ensembl-115"] = "Qualified" if p.translation_qualified else "Review"
        known = list(
            {
                k: v for gene in p.gene_ids for k, v in knowledge_by_gene.get(gene, {}).items()
            }.values()
        )
        values["gene-knowledge"] = len(known) or None
        exp = explanations.get(p.entity_id, [])
        ledger_rows = {}
        for source in ["cta_policy_findings", "external_ihc", "normal-induction-evidence"]:
            linked = [e for gene in p.gene_ids for e in ledgers[source].get(gene, [])]
            values[source] = len(linked) or None
            ledger_rows[source] = linked
        for source in current["papers"]:
            values[source] = "Listed"
        current["values"] = values
        current["literature_review"] = lit_by_label.get(p.label, [])
        current["review_status"] = (
            "Focused literature review"
            if current["literature_review"]
            else "Not individually literature-reviewed in this audit"
        )
        current["route_flags"] = {
            k: clean(f[k])
            for k in [
                "reproductive_rna_association",
                "reproductive_cell_association",
                "reproductive_protein_association",
                "developmental_association",
                "adult_profile_status",
                "developmental_profile_status",
                "Antibody concern",
                "reproductive_cell_ntpm",
                "reproductive_cell_indirect_ntpm",
                "oncoref_maximum_sample",
                "oncoref_cancer_code",
            ]
        }
        detail = {
            "normal": obs,
            "developmental": clean(drows),
            "normal_protein": clean(prows),
            "knowledge": clean(known),
            "ledgers": clean(ledger_rows),
            "explanations": clean(exp),
        }
        details[p.entity_id] = detail
        rows.append(current)
        if p.category == "CTA" and not current["papers"]:
            subset.append(current)
    data = {
        "release": "1.31.0",
        "date": "2026-11-06",
        "rows": rows,
        "sources": sources,
        "columns": columns,
        "numeric_columns": numeric_columns,
        "paper_sources": audit.sources.to_dict("records"),
        "literature": literature,
        "biology": biology,
        "catalog_review": json.loads((ASSETS / "catalog-review.json").read_text()),
    }
    data["biology_counts"] = [
        {
            **category,
            **{
                scope: sum(
                    category["id"] in r["biology_categories"]
                    for r in rows
                    if r["category"] == "CTA"
                    and (scope == "all" or bool(r["papers"]) == (scope == "registered"))
                )
                for scope in ["all", "registered", "unlisted"]
            },
        }
        for category in biology["categories"]
    ]
    data["comparison"] = comparison(rows)
    comparison_figure(data["comparison"])
    biology_figure(data["biology_counts"])
    data["tier_definitions"] = tissue_hazards.registry()["tier_definitions"]
    data["hazard_references"] = tissue_hazards.references().to_dict("records")
    data["hazard_assessments"] = tissue_hazards.assessments().to_dict("records")
    write_json(OUT / "data.json.gz", data, True)
    write_json(OUT / "unlisted-ctas.json", subset)
    write_json(OUT / "comparison.json", data["comparison"])
    write_json(OUT / "catalog-review.json", data["catalog_review"])
    write_json(OUT / "biology-categories.json", {**biology, "counts": data["biology_counts"]})
    flat = []
    for row in rows:
        flat.append(
            {
                **{
                    k: v
                    for k, v in row.items()
                    if k not in {"values", "literature_review", "route_flags"}
                },
                **row["values"],
                **row["route_flags"],
            }
        )
    for item in flat:
        for key, value in item.items():
            if isinstance(value, (list, dict, tuple)):
                item[key] = json.dumps(value, ensure_ascii=False)
    table = pd.DataFrame(flat)
    write_csv(table, OUT / "all-candidate-values.csv.gz")
    write_csv(table[table.category.eq("CTA")], OUT / "cta-biology-categories.csv")
    write_csv(
        table[table.entity_id.isin([r["entity_id"] for r in subset])], OUT / "unlisted-ctas.csv"
    )
    # Details are fetched only for opened proteins. Stable shards keep the
    # initial download small and avoid one request per measurement source.
    shards = {}
    for key, value in details.items():
        shard = hashlib.sha256(key.encode()).hexdigest()[:2]
        shards.setdefault(shard, {})[key] = value
    (OUT / "details").mkdir(exist_ok=True)
    for key, value in shards.items():
        write_json(OUT / "details" / (key + ".json.gz"), value, True)
    for name in ["report.html", "explorer.js", "style.css"]:
        shutil.copyfile(ASSETS / name, OUT / name)
    (OUT / "index.md").write_text(
        "# CTA evidence guide\n\n[Open the filterable candidate table and analysis of the 122 CTAs](report.html).\n\n[All candidate values](all-candidate-values.csv.gz) · [602 CTAs with biological categories](cta-biology-categories.csv) · [122-protein table](unlisted-ctas.csv) · [Comparison](comparison.json)\n"
    )
    write_json(
        OUT / "manifest.json",
        {
            "release": "1.31.0",
            "policy_id": cta_selection.coverage()["policy_id"],
            "counts": {
                "total": len(rows),
                "CTA": sum(r["category"] == "CTA" for r in rows),
                "unlisted_CTA": len(subset),
            },
            "inputs": {
                str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [
                    Path(__file__),
                    ROOT / "scripts/replay_cta_selection.py",
                    *sorted(ASSETS.glob("*")),
                    AUDIT / "protein-group-classifications.csv.gz",
                    AUDIT / "developmental-lineages.csv.gz",
                    AUDIT / "context-replay.json",
                    AUDIT / "hpa-protein-evidence.csv.gz",
                    AUDIT / "external-protein-evidence.csv.gz",
                    ROOT / "oncoref/data/cta-selection-human.csv.gz",
                    ROOT / "oncoref/data/cta-selection-profiles-human.csv.gz",
                    ROOT / "oncoref/data/gene-knowledge.json",
                    ROOT / "oncoref/data/normal-tissue-hazards.yaml",
                    ROOT / "oncoref/data/cta-policy-findings.csv",
                    ROOT / "oncoref/data/normal-induction-evidence.csv",
                    ROOT / "oncoref/data/cta-selection-sources.csv",
                    ROOT / "oncoref/data/cta-selection-audit.json",
                ]
            },
            "outputs": {
                str(p.relative_to(OUT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(OUT.rglob("*"))
                if p.is_file() and p.name != "manifest.json"
            },
        },
    )
    print(json.dumps(data["comparison"], indent=2))


def comparison(rows):
    result = {}
    for name, selected in [
        ("registered", [r for r in rows if r["category"] == "CTA" and r["papers"]]),
        ("unlisted", [r for r in rows if r["category"] == "CTA" and not r["papers"]]),
    ]:
        d = pd.DataFrame(selected)
        result[name] = {
            "proteins": len(d),
            "loci": sum(d.gene_ids.map(len)),
            "routes": d.route.value_counts().to_dict(),
            "tumor_max_median_tpm": d.tumor_max_tpm.median(),
            "reproductive_max_median_ntpm": d.reproductive_max_ntpm.median(),
            "developmental_max_median_tpm": d.developmental_lineage_max_tpm.median(),
            "normal_max_median_ntpm": d.normal_max_ntpm.median(),
            "unconfirmed_normal_signal": int(d.unconfirmed_normal_signal.sum()),
            "tumor_max_below_1": int(d.tumor_max_tpm.lt(1).sum()),
            "tumor_max_below_10": int(d.tumor_max_tpm.lt(10).sum()),
            "reviewed_tumor_protein": int(d.reliable_tumor_protein.sum()),
        }
    return clean(result)


def comparison_figure(comparison):
    """Export the descriptive route comparison, with an explicit denominator."""
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.hashsalt"] = "oncoref-cta-evidence-guide"
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 3.9))
    colors = ["#246994", "#96558a", "#138075"]
    routes = ["Adult reproductive only", "Developmental only", "Adult reproductive + developmental"]
    for y, key in enumerate(["registered", "unlisted"]):
        row = comparison[key]
        left = 0
        for route_name, color in zip(routes, colors, strict=True):
            count = row["routes"][route_name]
            width = 100 * count / row["proteins"]
            ax.barh(
                y, width, left=left, height=0.55, color=color, label=route_name if y == 0 else None
            )
            ax.text(
                left + width / 2,
                y,
                str(count),
                ha="center",
                va="center",
                color="white",
                fontsize=12,
            )
            left += width
    ax.set_yticks([0, 1], ["In nomination papers (480)", "Outside collection (122)"])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percent of proteins within each set; labels show protein counts")
    ax.set_title("How the two CTA sets qualify", loc="left", pad=22)
    ax.legend(loc="upper center", bbox_to_anchor=(0.4, -0.36), frameon=False, ncol=2)
    ax.spines[["right", "top", "left"]].set_visible(False)
    fig.subplots_adjust(left=0.26, bottom=0.34, top=0.8)
    for extension in ["svg", "png"]:
        fig.savefig(
            OUT / f"route-comparison.{extension}",
            dpi=160,
            metadata={"Date": None} if extension == "svg" else None,
        )
    plt.close(fig)


def biology_figure(counts):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matplotlib.rcParams["svg.hashsalt"] = "oncoref-cta-biology"
    active = [c for c in counts if c["all"]]
    fig, ax = plt.subplots(figsize=(9.5, 4.7))
    y = np.arange(len(active))
    for key, total, offset, color, label in [
        ("registered", 480, -0.18, "#315875", "In nomination papers (480)"),
        ("unlisted", 122, 0.18, "#087e8b", "Outside collection (122)"),
    ]:
        values = [100 * c[key] / total for c in active]
        ax.barh(y + offset, values, height=0.32, color=color, label=label)
        for row, value, category in zip(y, values, active, strict=True):
            ax.text(value + 0.6, row + offset, str(category[key]), va="center", fontsize=9)
    ax.set_yticks(y, [c["label"] for c in active])
    ax.invert_yaxis()
    ax.set_xlim(0, 108)
    ax.set_xlabel("Percent of each set; labels show protein counts. Categories overlap.")
    ax.set_title("Biological categories of the 602 CTAs", loc="left", pad=18)
    ax.spines[["right", "top", "left"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=2)
    fig.subplots_adjust(left=0.27, bottom=0.24, top=0.87)
    for extension in ["svg", "png"]:
        fig.savefig(
            OUT / f"biology-categories.{extension}",
            dpi=160,
            metadata={"Date": None} if extension == "svg" else None,
        )
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-contexts", action="store_true")
    build(**vars(parser.parse_args()))
