#!/usr/bin/env python3
"""Build the centralized literature registry and the 96-protein review offline.

This is an evidence/provenance release. It does not rewrite the frozen 1.31.0
selection or turn a publication mention into an admission or exclusion decision.
"""

from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import html
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "scripts/cta_evidence_assets"
OUT = ROOT / "docs/audits/cta-literature-20261107"
GUIDE = ROOT / "docs/audits/cta-evidence-guide-20261106"
REGISTRY = ROOT / "oncoref/data/cta-literature.json.gz"
LABELS = {
    "paper_evidence": "Gene-specific studies",
    "normal_review": "Normal-expression review",
    "identity_review": "Identity or translation review",
    "association_only": "Association or gene-list evidence",
    "related_only": "Related-species or family evidence",
    "no_specific_paper": "No specific study verified",
}
CATALOG_ROLES = {
    "10.1038/s41467-021-22695-y": ["reproductive_developmental"],
    "10.1371/journal.pone.0144155": ["reproductive_developmental"],
    "10.4161/cbt.12.3.16260": ["tumor_expression"],
    "10.1158/0008-5472.CAN-16-0225": ["tumor_expression", "function"],
    "10.1016/s1074-7613(95)80053-0": ["catalog_nomination", "immune_recognition"],
    "10.1038/s41467-020-19141-w": ["catalog_nomination", "immune_recognition"],
}
RETAIN = {
    "10.1038/s41467-020-19141-w": "Placenta/cancer candidates plus VGLL1-specific HLA-peptide and T-cell experiments; the functional result is not inherited by every listed gene.",
    "10.1016/s1074-7613(95)80053-0": "Historical BAGE discovery and T-cell recognition. The mapped MAGE context mentions do not replace the unresolved BAGE target mapping.",
    "10.1158/0008-5472.CAN-16-0225": "EGFL6 tumor/vascular mechanism; it supplies no reproductive association by itself.",
}


def read(path):
    content = path.read_bytes()
    return json.loads(gzip.decompress(content) if path.suffix == ".gz" else content)


def write_json(path, value):
    payload = (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode()
    path.write_bytes(gzip.compress(payload, mtime=0) if path.suffix == ".gz" else payload)


def compact_cover(sets):
    """Deterministic, inclusion-minimal set cover; never claim global optimality."""
    sets = {key: set(values) for key, values in sets.items() if values}
    universe = set().union(*sets.values()) if sets else set()
    owners = defaultdict(list)
    for key, values in sets.items():
        for value in values:
            owners[value].append(key)
    selected = {keys[0] for keys in owners.values() if len(keys) == 1}
    covered = set().union(*(sets[k] for k in selected)) if selected else set()
    while covered != universe:
        key = min(
            (k for k in sets if k not in selected), key=lambda k: (-len(sets[k] - covered), k)
        )
        selected.add(key)
        covered.update(sets[key])
    for key in sorted(selected, key=lambda k: (len(sets[k]), k)):
        others = set().union(*(sets[k] for k in selected - {key})) if len(selected) > 1 else set()
        if sets[key] <= others:
            selected.remove(key)
    return sorted(selected)


def identity_maps(rows):
    ids, symbols, aliases = defaultdict(set), defaultdict(set), defaultdict(set)
    for row in rows:
        for gene in row["gene_identities"]:
            gid = gene["gene_id"]
            ids[gid].add(row["entity_id"])
            # The released identity field uses `symbol`, not the combined protein label.
            symbol = gene.get("gene_name", gene.get("symbol", ""))
            if symbol:
                symbols[symbol].add((row["entity_id"], gid))
            for alias in gene.get("synonyms") or []:
                aliases[alias].add((row["entity_id"], gid))
    return ids, symbols, aliases


def resolve_program(row, ids, symbols, aliases):
    if row["type"] != "gene":
        return [], "non_gene_element"
    gid, name = row["gene_id"], row["symbol"]
    if gid.startswith("ENSG"):
        matches = [(entity, gid) for entity in ids.get(gid.split(".")[0], [])]
        return matches, "exact_ensembl" if matches else "unmapped_ensembl"
    # The historical H3.Y locus is explicitly resolved, rather than guessed by punctuation.
    if name == "H3.Y":
        return [
            (entity, "ENSG00000269466") for entity in ids.get("ENSG00000269466", [])
        ], "reviewed_historical_identity"
    matches = symbols.get(name, set())
    if len({gid for _, gid in matches}) == 1:
        return sorted(matches), "exact_symbol"
    if matches:
        return [], "ambiguous_symbol"
    matches = aliases.get(name, set())
    if len({gid for _, gid in matches}) == 1:
        return sorted(matches), "unique_alias"
    return [], "ambiguous_alias" if matches else "unmapped_symbol"


def build_registry():
    guide = read(GUIDE / "data.json.gz")
    notes = read(ASSETS / "individual-review-notes.json")
    bibliography = read(ASSETS / "review-bibliography.json")
    policy = read(ASSETS / "literature-source-policy.json")
    programs = read(ASSETS / "program-source-rows.json.gz")
    old_review = read(ASSETS / "literature-review.json")
    population = read(OUT / "review-population.json")
    baseline = {r["label"]: r for r in population}
    rows = {r["entity_id"]: r for r in guide["rows"]}
    by_label = {r["label"]: r for r in guide["rows"]}
    sources, edges, mapped_programs = {}, [], []
    with (ROOT / "oncoref/data/cta-publication-sources.csv").open() as handle:
        titles = {r["doi"].lower(): r["title"] for r in csv.DictReader(handle)}
    for source in guide["paper_sources"]:
        key = source["source_id"]
        sources[key] = {
            "source_id": key,
            "citation": source["citation"],
            "title": titles.get(key.lower(), source["citation"]),
            "url": source["source_url"],
            "kind": "historical_catalog",
            "roles": CATALOG_ROLES.get(key, ["catalog_nomination"]),
            "scope": source["evidence_scope"],
            "inspection": "released_source_rows",
            "retention_reason": RETAIN.get(
                key, "Retain the original list, assay scope and independent source provenance."
            ),
        }
    with gzip.open(
        ROOT / "docs/audits/cta-source-walkthrough-20261028/human-memberships.csv.gz", "rt"
    ) as handle:
        for row in csv.DictReader(handle):
            edges.append(
                {
                    "source_id": row["source_id"],
                    "entity_id": row["entity_id"],
                    "protein": row["label"],
                    "kind": "historical_catalog",
                    "gene_ids": ast.literal_eval(row["gene_ids"]),
                    "eligible_for_coverage": True,
                    "scope": "Released mapped membership; see original source-table roles.",
                }
            )
    for catalog in read(ASSETS / "catalog-review.json")["catalogs"]:
        key = (
            "10.1016/j.devcel.2019.06.011"
            if catalog["id"] == "chew2019"
            else "10.1016/j.celrep.2023.113114"
        )
        sources[key] = {
            "source_id": key,
            "citation": catalog["citation"],
            "title": catalog["title"],
            "url": catalog["url"],
            "kind": "program_catalog",
            "roles": ["cancer_embryo_program"],
            "scope": catalog["model"],
            "inspection": "complete_supplementary_tables",
            "source_files": catalog["sources"],
            "retention_reason": "Defined cancer–embryo program membership, with the source rows and upstream gene-set origin retained.",
        }
    upstream = "10.1016/j.stem.2022.01.014"
    sources[upstream] = {
        "source_id": upstream,
        "citation": "Taubenschmid-Stowers et al. (2022)",
        "title": "8C-like cells capture the human zygotic genome activation program in vitro",
        "url": "https://pubmed.ncbi.nlm.nih.gov/35216671/",
        "roles": ["upstream_program"],
        "kind": "upstream_program",
        "inspection": "provenance_in_reused_table",
        "scope": "Original ZGA/8C gene sets reused by Smith 2023 Table S3A.",
        "retention_reason": "Credit the original program definition; Smith's reuse is not an independent replication of gene-set membership.",
    }
    maps = identity_maps(guide["rows"])
    for row in programs["rows"]:
        matches, method = resolve_program(row, *maps)
        record = {
            **row,
            "match_method": method,
            "matches": [{"entity_id": e, "gene_id": g} for e, g in matches],
        }
        if method == "reviewed_historical_identity":
            record["identity_reference"] = "https://www.ncbi.nlm.nih.gov/research/histonedb/human/"
        mapped_programs.append(record)
        for entity, gid in matches:
            edges.append(
                {
                    "source_id": row["source_id"],
                    "entity_id": entity,
                    "protein": rows[entity]["label"],
                    "gene_ids": [gid],
                    "kind": "program_membership",
                    "eligible_for_coverage": True,
                    "scope": row["program"],
                    "table": row["table"],
                    "row": row["row"],
                    "match_method": method,
                }
            )
    # Previously reviewed papers remain in this one registry, with their original scope.
    for finding in old_review["findings"]:
        key = finding["url"]
        if "proteinatlas.org" in key:
            kind, roles = "measurement_annotation", ["normal_tissue"]
        else:
            kind = "previous_individual_review"
            category = finding["category"].lower()
            roles = (
                ["normal_tissue"]
                if any(s in category for s in ["normal", "sensory"])
                else ["tumor_expression"]
                if "fusion" in category
                else ["catalog_nomination", "tumor_expression"]
                if "missing" in category
                else ["reproductive_developmental"]
            )
        sources[key] = {
            "source_id": key,
            "citation": finding["title"],
            "title": finding["title"],
            "url": key,
            "kind": kind,
            "roles": roles,
            "scope": finding["model"] + "; " + finding["limitation"],
            "inspection": "previous_focused_review",
            "retention_reason": finding["finding"],
        }
        for label in finding["proteins"]:
            row = by_label[label]
            # These old findings explicitly require gene/species resolution before exact attribution.
            unresolved = (
                label in {"PATE1", "SP7", "CALHM3"} or "fusion" in finding["category"].lower()
            )
            edges.append(
                {
                    "source_id": key,
                    "entity_id": row["entity_id"],
                    "protein": label,
                    "gene_ids": row["gene_ids"],
                    "kind": kind,
                    "eligible_for_coverage": not unresolved and kind != "measurement_annotation",
                    "scope": finding["finding"],
                    "identity_limit": finding["limitation"],
                }
            )
    pmid_sources = {}
    for pmid, paper in bibliography.items():
        key = paper.get("doi") or "PMID:" + pmid
        pmid_sources[pmid] = key
        assessment = policy["paper_assessments"][pmid]
        author = paper["authors"].split(",")[0]
        sources[key] = {
            **paper,
            **assessment,
            "source_id": key,
            "citation": f"{author} et al. ({paper['year']})",
            "kind": "individual_review",
            "scope": assessment.get(
                "model_limit", "See each linked protein's assay and identity assessment."
            ),
            "retention_reason": assessment["contribution"],
        }
    broad_ids = {r["gene_id"] for r in programs["screen_only"]["rows"]}
    reviews = []
    for note in notes["reviews"]:
        row = baseline[note["protein"]]
        reviewed = {
            **note,
            "entity_id": row["entity_id"],
            "gene_ids": row["gene_ids"],
            "gene_identities": row["gene_identities"],
            "biology_categories": row["biology_categories"],
            "admission_route": row["route"],
            "released_category": "CTA",
            "outcome_label": LABELS[note["outcome"]],
            "baseline_gap": "broad_table_only" if set(row["gene_ids"]) & broad_ids else "no_match",
            "citations": [],
            "review_disposition": "Record biological evidence",
        }
        if note["outcome"] == "normal_review":
            reviewed["review_disposition"] = (
                "Adjudicate normal-expression evidence before the next classification release"
            )
        elif note["outcome"] == "identity_review":
            reviewed["review_disposition"] = "Resolve locus, proteoform or translation evidence"
        elif note["outcome"] in {"related_only", "no_specific_paper"}:
            reviewed["review_disposition"] = (
                "Keep measured-expression support; no exact human study is credited"
            )
        for pmid in note["pmids"]:
            key = pmid_sources[pmid]
            assessment = policy["paper_assessments"][pmid]
            limit = policy["edge_limits"].get(f"{note['protein']}:{pmid}", "exact_human_gene")
            if assessment.get("model_limit"):
                limit = assessment["model_limit"]
            # A named human proteomic/association study supplies paper coverage
            # even when reagent/peptide specificity is still under review.
            eligible = (
                limit
                in {
                    "exact_human_gene",
                    "reagent_identity_unresolved",
                    "peptide_identity_unresolved",
                    "native_expression_annotation_link",
                }
                and note["outcome"] != "related_only"
            )
            link = {
                "source_id": key,
                "pmid": pmid,
                "url": bibliography[pmid]["url"],
                "title": bibliography[pmid]["title"],
                "roles": assessment["roles"],
                "inspection": assessment["inspection"],
                "identity_scope": limit,
                "eligible_for_coverage": eligible,
            }
            reviewed["citations"].append(link)
            edges.append(
                {
                    "source_id": key,
                    "entity_id": row["entity_id"],
                    "protein": note["protein"],
                    "gene_ids": row["gene_ids"],
                    "kind": "individual_review",
                    "scope": limit,
                    "eligible_for_coverage": eligible,
                }
            )
        reviewed["exact_human_paper_coverage"] = any(
            c["eligible_for_coverage"] for c in reviewed["citations"]
        )
        reviews.append(reviewed)
    assert Counter(r["baseline_gap"] for r in reviews) == {"no_match": 41, "broad_table_only": 55}
    sets = defaultdict(set)
    links = defaultdict(set)
    for edge in edges:
        links[edge["source_id"]].add(edge["entity_id"])
        if edge["eligible_for_coverage"]:
            sets[edge["source_id"]].add(edge["entity_id"])
    core = set(compact_cover(sets))
    all_covered = set().union(*sets.values())
    old_keys = {s["source_id"] for s in guide["paper_sources"]}
    old_sets = {k: sets[k] for k in old_keys}
    cta_ids = {r["entity_id"] for r in guide["rows"] if r["category"] == "CTA"}
    unlisted_ids = {
        r["entity_id"] for r in guide["rows"] if r["category"] == "CTA" and not r["papers"]
    }
    for key, source in sources.items():
        others = set().union(*(values for other, values in sets.items() if other != key))
        core_others = {other: sets[other] & sets[key] for other in core if other != key}
        source.update(
            collection="core" if key in core else "supplementary",
            linked_proteins=len(links[key]),
            covered_proteins=len(sets[key]),
            covered_ctas=len(sets[key] & cta_ids),
            covered_unlisted_ctas=len(sets[key] & unlisted_ids),
            unique_proteins=len(sets[key] - others),
            unique_ctas=len((sets[key] - others) & cta_ids),
            covered_by=compact_cover(core_others) if key not in core and sets[key] else [],
            proteins=sorted(rows[e]["label"] for e in links[key]),
            coverage_proteins=sorted(rows[e]["label"] for e in sets[key]),
        )
        if key in core:
            remaining = set().union(*(sets[k] for k in core - {key}))
            source["core_witnesses"] = sorted(rows[e]["label"] for e in sets[key] - remaining)
        else:
            source["core_witnesses"] = []
    coverage = []
    for entity in sorted(cta_ids, key=lambda e: rows[e]["label"]):
        row = rows[entity]
        matches = [e for e in edges if e["entity_id"] == entity]
        coverage.append(
            {
                "entity_id": entity,
                "protein": row["label"],
                "gene_ids": row["gene_ids"],
                "historical_sources": row["papers"],
                "exact_paper_sources": sorted(
                    {e["source_id"] for e in matches if e["eligible_for_coverage"]}
                ),
                "other_sources": sorted(
                    {e["source_id"] for e in matches if not e["eligible_for_coverage"]}
                ),
                "biology_categories": row["biology_categories"],
                "admission_route": row["route"],
            }
        )
    counts = {
        "reviewed": len(reviews),
        "baseline_no_match": 41,
        "baseline_broad_only": 55,
        "outcomes": dict(Counter(r["outcome"] for r in reviews)),
        "exact_human_studies_in_96": sum(r["exact_human_paper_coverage"] for r in reviews),
        "with_citations_in_96": sum(bool(r["citations"]) for r in reviews),
        "core_papers": len(core),
        "catalog_core_papers": sum(
            sources[k]["kind"] in {"historical_catalog", "program_catalog"} for k in core
        ),
        "individual_core_papers": sum(
            sources[k]["kind"] in {"individual_review", "previous_individual_review"} for k in core
        ),
        "supplementary_sources": len(sources) - len(core),
        "covered_proteins": len(all_covered),
        "covered_ctas": len(all_covered & cta_ids),
        "covered_unlisted_ctas": len(all_covered & unlisted_ids),
        "ctas_without_exact_paper_coverage": len(cta_ids - all_covered),
        "historical_catalog_proteins": len(set().union(*old_sets.values())),
        "historical_catalog_core": compact_cover(old_sets),
        "historical_cta_core": compact_cover({k: v & cta_ids for k, v in old_sets.items()}),
    }
    inputs = [
        Path(__file__),
        GUIDE / "data.json.gz",
        ASSETS / "individual-review-notes.json",
        ASSETS / "review-bibliography.json",
        ASSETS / "literature-source-policy.json",
        ASSETS / "program-source-rows.json.gz",
        ASSETS / "literature-review.json",
        ASSETS / "catalog-review.json",
        OUT / "search-ledger.json",
        OUT / "review-population.json",
        ROOT / "docs/audits/cta-source-walkthrough-20261028/human-memberships.csv.gz",
        ROOT / "oncoref/data/cta-publication-sources.csv",
        *sorted((ASSETS / "literature").glob("*")),
    ]
    return {
        "schema": "oncoref.cta-literature",
        "version": 1,
        "date": "2026-11-07",
        "baseline_release": "1.31.0",
        "classification_effect": "Evidence annotation only; findings requiring adjudication are not yet applied to the frozen selection.",
        "roles": policy["roles"],
        "coverage_rule": policy["coverage_rule"],
        "core_rule": policy["core_rule"],
        "counts": counts,
        "outcome_labels": LABELS,
        "sources": sorted(
            sources.values(),
            key=lambda s: (
                s["collection"],
                {
                    "historical_catalog": 0,
                    "program_catalog": 1,
                    "previous_individual_review": 2,
                    "individual_review": 3,
                    "upstream_program": 4,
                    "measurement_annotation": 5,
                }[s["kind"]],
                s["citation"],
            ),
        ),
        "reviews": reviews,
        "memberships": edges,
        "cta_coverage": coverage,
        "program_rows": mapped_programs,
        "background_screen": {k: v for k, v in programs["screen_only"].items() if k != "rows"}
        | {"gene_rows": len(broad_ids)},
        "measurement_sources": [s for s in guide["sources"] if s["source_id"] not in old_keys],
        "source_files": programs["source_files"],
        "inputs": [
            {"path": str(p.relative_to(ROOT)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in inputs
        ],
    }


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                    for k, v in row.items()
                    if k in fields
                }
            )


def main():
    data = build_registry()
    OUT.mkdir(exist_ok=True)
    write_json(REGISTRY, data)
    write_json(OUT / "registry.json.gz", data)
    # The browser loads a small view; the complete occurrence ledger is a separate download.
    view = {k: v for k, v in data.items() if k not in {"memberships", "program_rows"}}
    write_json(OUT / "review.json", view)
    write_csv(
        OUT / "individual-reviews.csv",
        data["reviews"],
        [
            "protein",
            "gene_ids",
            "baseline_gap",
            "outcome_label",
            "summary",
            "identity",
            "pmids",
            "exact_human_paper_coverage",
            "review_disposition",
        ],
    )
    write_csv(
        OUT / "sources.csv",
        data["sources"],
        [
            "source_id",
            "citation",
            "title",
            "url",
            "collection",
            "roles",
            "inspection",
            "covered_proteins",
            "covered_ctas",
            "covered_unlisted_ctas",
            "unique_proteins",
            "unique_ctas",
            "covered_by",
            "retention_reason",
        ],
    )
    write_csv(
        OUT / "cta-paper-coverage.csv",
        data["cta_coverage"],
        [
            "protein",
            "gene_ids",
            "historical_sources",
            "exact_paper_sources",
            "other_sources",
            "biology_categories",
            "admission_route",
        ],
    )
    for name in ["report.html", "review.js", "style.css"]:
        shutil.copyfile(ASSETS / "literature" / name, OUT / name)
    write_figures(data)
    write_index(data)
    integrate_evidence_guide()
    write_json(
        OUT / "manifest.json",
        {
            "date": data["date"],
            "inputs": data["inputs"],
            "counts": data["counts"],
            "outputs": {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(OUT.iterdir())
                if p.is_file() and p.name != "manifest.json"
            },
        },
    )
    print(json.dumps(data["counts"], indent=2))


def integrate_evidence_guide():
    """Refresh only the presentation overlay, preserving frozen measurements."""
    for name in ["report.html", "explorer.js"]:
        shutil.copyfile(ASSETS / name, GUIDE / name)
    manifest = read(GUIDE / "manifest.json")
    for name in ["report.html", "explorer.js"]:
        manifest["inputs"][str((ASSETS / name).relative_to(ROOT))] = hashlib.sha256(
            (ASSETS / name).read_bytes()
        ).hexdigest()
        manifest["outputs"][name] = hashlib.sha256((GUIDE / name).read_bytes()).hexdigest()
    (GUIDE / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")) + "\n")


def write_figures(data):
    # A static vector figure can be exported without a browser or plotting dependency.
    c = data["counts"]
    colors = ["#216b72", "#a14632", "#96651d", "#5f6681", "#718073", "#889298"]
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 360" role="img"><title>Individual review of 96 CTA proteins</title><rect width="1000" height="360" fill="#f7f5ef"/><g font-family="system-ui,sans-serif" fill="#1d3033"><text x="30" y="40" font-size="25">What the individual review found</text><text x="30" y="68" font-size="15">41 previously unmatched + 55 broad-table-only proteins · 7 November 2026</text>'
    ]
    for i, (key, label) in enumerate(LABELS.items()):
        n = c["outcomes"][key]
        y = 103 + i * 40
        parts.append(
            f'<text x="30" y="{y + 18}" font-size="15">{html.escape(label)}</text><rect x="345" y="{y}" width="{n * 14}" height="25" rx="3" fill="{colors[i]}"/><text x="{355 + n * 14}" y="{y + 18}" font-size="16">{n}</text>'
        )
    parts.append("</g></svg>")
    (OUT / "review-outcomes.svg").write_text("".join(parts))
    flow = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 480" role="img"><title>Distinct roles of measurement and literature sources</title><defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8" fill="#216b72"/></marker></defs><rect width="1200" height="480" fill="#f7f5ef"/><g font-family="system-ui,sans-serif" fill="#203338"><text x="30" y="40" font-size="25">How the sources contribute</text><text x="30" y="68" font-size="15">Classification is based on measurements and adjudicated findings. Paper coverage is recorded separately.</text>
<g fill="#e6eeeb" stroke="#a8c2b7"><rect x="30" y="105" width="250" height="115" rx="5"/><rect x="320" y="105" width="250" height="115" rx="5"/><rect x="610" y="105" width="290" height="115" rx="5"/><rect x="940" y="105" width="230" height="115" rx="5"/></g>
<g font-size="17" font-weight="650"><text x="48" y="135">Reference universe</text><text x="338" y="135">Presence and association</text><text x="628" y="135">Admission and exclusion</text><text x="958" y="135">Released classifications</text></g>
<g font-size="14"><text x="48" y="165">19,863 protein records</text><text x="48" y="190">Locus + selected coding sequence</text><text x="338" y="165">Cancer expression</text><text x="338" y="190">Reproductive / developmental RNA</text><text x="628" y="165">Normal cells, coverage, translation</text><text x="628" y="190">Protein evidence + sequence overlap</text><text x="958" y="165">602 CTA · 1,006 Other</text><text x="958" y="190">18,255 Unsafe · release 1.31.0</text></g>
<g stroke="#216b72" stroke-width="2" marker-end="url(#arrow)"><path d="M282,163 H313"/><path d="M573,163 H603"/><path d="M903,163 H933"/></g>
<g fill="#fff" stroke="#cbd5cf"><rect x="30" y="290" width="350" height="145" rx="5"/><rect x="410" y="290" width="350" height="145" rx="5"/><rect x="790" y="290" width="380" height="145" rx="5"/></g>
<g font-size="18" font-weight="650"><text x="48" y="323">Catalog and program papers</text><text x="428" y="323">Individual studies</text><text x="808" y="323">Core and supplementary sets</text></g>
<g font-size="14"><text x="48" y="355">Who nominated or studied the candidate?</text><text x="48" y="380">Preserve exact source rows and gene identity.</text><text x="48" y="405">Membership does not grant admission.</text><text x="428" y="355">What was measured, where, and in which model?</text><text x="428" y="380">Record normal findings and identity conflicts.</text><text x="428" y="405">Adjudicated findings enter the selection rules.</text><text x="808" y="355">Core preserves mapped paper coverage.</text><text x="808" y="380">Supplementary retains additional evidence.</text><text x="808" y="405">Neither label ranks safety or evidence strength.</text></g></g></svg>"""
    (OUT / "source-roles.svg").write_text(flow)


def write_index(data):
    c = data["counts"]
    (OUT / "index.md").write_text(f"""# CTA literature and source roles

[Open the filterable report](report.html). The review covers all **96 proteins**: 41 previously unmatched and 55 previously found only in a broad screening table.

{c["with_citations_in_96"]} have linked primary or related papers; {c["exact_human_studies_in_96"]} have an exact-human study credited under the declared coverage rule. Evidence types are shown separately: a DNA association is not a protein experiment.

The harmonized collection has {c["core_papers"]} core papers and {c["supplementary_sources"]} supplementary sources. Its mapped coverage reaches {c["covered_ctas"]} of the 602 released CTA records. {c["ctas_without_exact_paper_coverage"]} remain without qualifying paper coverage. The expression and reference datasets remain the basis of selection.

The 1.31.0 classifications are frozen. Normal-expression and identity findings in this review require adjudication before the next classification release.

- [Individual reviews](individual-reviews.csv)
- [Core and supplementary sources, roles and overlap](sources.csv)
- [Paper coverage for all 602 CTAs](cta-paper-coverage.csv)
- [Complete registry, including source-row occurrences](registry.json.gz)
- [Search queries and retrieval receipts](search-ledger.json)
- [Figure: review outcomes](review-outcomes.svg)
- [Reproduction manifest](manifest.json)

Rebuild offline with `.venv/bin/python scripts/build_cta_literature_review.py`. Refresh indexed searches with `scripts/review_unlisted_cta_literature.py` and cited records with `scripts/fetch_cta_individual_review.py`; neither script automatically accepts a finding. The pinned workbook importer is `scripts/import_cta_program_tables.py`.
""")


if __name__ == "__main__":
    main()
