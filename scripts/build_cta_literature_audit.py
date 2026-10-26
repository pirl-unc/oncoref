#!/usr/bin/env python3
"""Join every nomination to search coverage, source flags and explicit adjudications.

This is a review overlay. It never mutates production gene/protein categories and
never promotes an unreviewed search hit into an accepted expression finding.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "docs/audits/cta-literature-20261025"
NO_REVIEW = "No individual clinical-outcome review in this report"
OFF_TARGET = r"off.target|cross.reactiv|unspecific|non.specific"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path):
    return pd.read_csv(path, keep_default_na=False)


def load_adjudications(path, symbols):
    result = {}
    for record in yaml.safe_load(path.read_text())["records"]:
        for gene in record["genes"]:
            if gene not in symbols or gene in result:
                raise ValueError(f"Unknown/duplicate adjudication: {gene}")
            if record["proposed"] not in {"CTA", "Neither", "Unsafe", "Unchanged"}:
                raise ValueError(f"Invalid proposed category: {gene}")
            result[gene] = record
    return result


def primary_signals(atlas, source):
    primary = read_table(source)
    # normal_signal already includes the reviewed anatomy/ignored-reagent logic.
    normal = primary.normal_signal.astype(str).str.lower().eq("true")
    retained = set(atlas.loc[atlas.Category.eq("CTA"), "Gene ID"])
    selected = primary.loc[normal & primary.gene_id.isin(retained)].copy()
    summaries = {}
    for gene_id, group in selected.groupby("gene_id"):
        strong = group.category.eq("vital") & group.level.isin(["Medium", "High"])
        summaries[gene_id] = {
            "Primary normal staining rows": len(group),
            "Primary normal antibodies": group.Antibody.nunique(),
            "Primary vital Medium/High rows": int(strong.sum()),
            "Off-target warning in profile": bool(
                group["Profile notes"].str.contains(OFF_TARGET, case=False, regex=True).any()
            ),
            "Primary staining contexts": " | ".join(
                sorted(
                    {
                        f"{r.tissue}/{r.cell_type}: {r.level} ({r.Antibody})"
                        for r in group.itertuples()
                    }
                )
            ),
        }
    return selected, summaries


def search_evidence(ledger, cache, verify=True):
    papers = {}
    for row in ledger.to_dict("records"):
        if row["Status"] != "Retrieved":
            continue
        path = cache / row["Cache file"]
        if verify and digest(path) != row["SHA256"]:
            raise ValueError(f"Search cache checksum mismatch: {path}")
        with gzip.open(path, "rt") as handle:
            record = json.load(handle)
        for paper in record["response"]["resultList"]["result"]:
            key = row["Gene ID"], paper["id"]
            if key not in papers:
                papers[key] = {
                    "Gene": row["Gene"],
                    "Gene ID": row["Gene ID"],
                    "Category": row["Category"],
                    "PMID": paper["id"],
                    "Title": re.sub(r"<[^>]*>", "", paper["title"]),
                    "Year": paper.get("pubYear", ""),
                    "DOI": paper.get("doi", ""),
                    "URL": f"https://pubmed.ncbi.nlm.nih.gov/{paper['id']}/",
                    "Queries": set(),
                    "Review": "Retrieved lead; gene identity, species, assay and conclusions not adjudicated",
                }
            papers[key]["Queries"].add(row["Search"])
    rows = []
    for record in papers.values():
        record["Queries"] = " | ".join(sorted(record["Queries"]))
        rows.append(record)
    return pd.DataFrame(rows)


def baseline_assessment(row):
    if row["Category"] == "CTA":
        return (
            "Provisional screen pass",
            "Meets the frozen measured-data screen; independent literature adjudication is incomplete.",
        )
    if row["Category"] == "Neither":
        return (
            "Provisional hold consistent with policy",
            "Keep outside targets and default forbidden proteins pending the stated evidence gap; this is not a safety endorsement.",
        )
    return (
        "Exclusion consistent with policy",
        "The recorded normal-expression or safety gate supports exclusion under this policy; RNA alone does not prove clinical injury.",
    )


def build_rows(atlas, ledger, decisions, signals):
    searches = dict(iter(ledger.groupby("Gene ID")))
    rows = []
    for original in atlas.to_dict("records"):
        row = dict(original)
        gene, gene_id = row["Gene"], row["Gene ID"]
        local = searches.get(gene_id, pd.DataFrame(columns=ledger.columns))
        successful = local.loc[local.Status.eq("Retrieved")]
        completed = {name for name in successful.Search if name in {"normal", "therapy"}}
        assessment, implication = baseline_assessment(row)
        signal = signals.get(gene_id, {})
        if row["Category"] == "CTA" and signal:
            assessment = "Reconcile primary antibodies"
            implication = "Curated pass conflicts with individual normal antibody records; inspect specificity before treating raw staining as accepted protein."
        if row["Category"] == "CTA" and row["Literature conflict"]:
            assessment = "Resolve existing literature conflict"
            implication = row["Literature follow-up"]
        row.update(
            {
                "Assessment": assessment,
                "Proposed category": row["Category"],
                "Change proposed": False,
                "Assessment basis": "Source flags and policy only; new papers not individually adjudicated",
                "New finding": "",
                "Implication": implication,
                "Next action": row["Reconsideration evidence"],
                "New review PMIDs": "",
                "New review sources": "",
                "Search coverage": "Both searches retrieved"
                if len(completed) == 2
                else "Incomplete",
                "Successful searches": len(successful),
                "Failed searches": int(local.Status.eq("Failed").sum()),
                "Search result count": int(pd.to_numeric(successful.Hits, errors="coerce").sum()),
                "Truncated searches": int(
                    successful.Truncated.astype(str).str.lower().eq("true").sum()
                ),
                "Zero-hit topic searches": int(
                    pd.to_numeric(
                        successful.loc[successful.Search.isin(["normal", "therapy"]), "Hits"],
                        errors="coerce",
                    )
                    .eq(0)
                    .sum()
                ),
                "Primary normal staining rows": 0,
                "Primary normal antibodies": 0,
                "Primary vital Medium/High rows": 0,
                "Off-target warning in profile": False,
                "Primary staining contexts": "",
            }
        )
        row.update(signal)
        if gene in decisions:
            decision = decisions[gene]
            proposed = decision["proposed"]
            row.update(
                {
                    "Assessment": decision["assessment"],
                    "Proposed category": row["Category"] if proposed == "Unchanged" else proposed,
                    "Assessment basis": decision["depth"],
                    "New finding": decision["finding"],
                    "Implication": decision["implication"],
                    "Next action": decision["next"],
                    "New review PMIDs": "|".join(map(str, decision["pmids"])),
                    "New review sources": " | ".join(
                        f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in decision["pmids"]
                    ),
                }
            )
            row["Change proposed"] = row["Proposed category"] != row["Category"]
        rows.append(row)
    result = pd.DataFrame(rows)
    # These checks detect accidental relabelling or a broken exclusion overlay.
    if not result["Gene ID"].is_unique or len(result) != 2540:
        raise ValueError("Review must contain every mapped nomination exactly once")
    maximum = (
        result[["Normal RNA maximum (TPM)", "Normal RNA maximum (nTPM)"]]
        .apply(pd.to_numeric, errors="coerce")
        .max(axis=1)
    )
    if ((maximum > 30) & result.Category.ne("Unsafe")).any():
        raise ValueError("An existing non-Unsafe gene exceeds the approved hard normal-RNA limit")
    return result


def write_explorer(result, summary, path):
    columns = [
        "Gene",
        "Gene ID",
        "Category",
        "Detail",
        "Proposed category",
        "Assessment",
        "Assessment basis",
        "Reproductive RNA (nTPM)",
        "Normal RNA maximum (nTPM)",
        "Normal RNA maximum (TPM)",
        "Highest normal nTPM context",
        "Highest normal TPM context",
        "Reason",
        "New finding",
        "Implication",
        "Next action",
        "New review sources",
        "Search coverage",
        "Truncated searches",
        "Primary normal staining rows",
        "Primary vital Medium/High rows",
        "Off-target warning in profile",
        "Primary staining contexts",
        "Literature conflict",
        "Literature conflict source",
        "Curated tissue biology",
        "Biological function",
        "Tissue PubMed IDs",
        "Reviewed literature",
        "Reviewed literature sources",
        "Boundary review",
        "Cancer RNA at >=5 pTPM",
        "Cancer IHC Medium/High",
        "Cancer RNA cohorts",
        "Cancer IHC groups",
        "Immune RNA findings",
        "Cell RNA findings",
        "Antibody findings",
    ]
    payload = json.dumps(result[columns].to_dict("records"), ensure_ascii=False).replace(
        "<", "\\u003c"
    )
    page = r"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CTA literature audit · 25 October 2026</title><style>
body{font:16px/1.55 system-ui,sans-serif;color:#192b3c;background:#f5f6f7;margin:0}main{max-width:1150px;margin:auto;padding:32px 24px}
h1{font-size:36px;line-height:1.15}h2{font-size:23px}a{color:#096e77}header p{max-width:950px}.notice{background:#fff2d5;border-left:5px solid #b86d20;padding:18px}
.stats{display:flex;gap:25px;flex-wrap:wrap;margin:25px 0}.stats b{font-size:30px;display:block}.filters{position:sticky;top:0;background:#f5f6f7;padding:15px 0;display:flex;gap:12px;flex-wrap:wrap;border-bottom:1px solid #bbc5ca}
input,select,button{font:inherit;padding:9px;border:1px solid #a9b5bb;border-radius:4px}input{flex:1;min-width:200px}button{cursor:pointer;background:white}
article{background:white;margin:16px 0;padding:20px;border:1px solid #dae0e3;border-radius:6px}article h2{margin:0}.tag{font-size:13px;background:#e7eff1;padding:4px 8px;border-radius:4px;display:inline-block;margin:8px 5px 0 0}
.change{background:#ffe6dd}dl{display:grid;grid-template-columns:220px 1fr;gap:8px 18px}dt{font-weight:650}dd{margin:0;overflow-wrap:anywhere}details{margin-top:14px}summary{cursor:pointer;font-weight:650}.muted{color:#596a74}pre{white-space:pre-wrap}footer{margin-top:25px}
@media(max-width:700px){dl{display:block}dt{margin-top:10px}.stats{gap:14px}main{padding:20px 12px}h1{font-size:28px}}
</style><main><header><p>OncoRef · research review overlay · 25 October 2026</p><h1>Which CTA decisions withstand closer review?</h1>
<p class="notice"><b>All 2,540 genes have a review row. Search coverage is not completed paper adjudication.</b> Published target sets are unchanged. CTA means a research candidate, not established clinical safety; Unsafe includes precautionary expression exclusions, not necessarily observed injury.</p>
<div class="stats"><span><b>399</b>CTA</span><span><b>561</b>Neither</span><span><b>1,580</b>Unsafe</span><span><b>__CHANGES__</b>proposed changes</span><span><b>__SEARCHED__</b>genes searched</span></div>
<p>Five proposed exclusions: <b>ABCC11, ACTRT1, CCDC158, FTMT, YBX2</b>. Three proposed holds: <b>LGALS13, MAGEC2, PLCZ1</b>. Primary antibody discrepancies affect 180 retained gene assignments; some involve explicitly discounted off-target staining or shared family reagents.</p>
<p><a href="all-genes.csv">All-gene review table</a> · <a href="proposed-changes.csv">Proposed changes</a> · <a href="papers.csv">Retrieved paper leads</a> · <a href="searches.csv">Exact search log</a> · <a href="index.md">Methods and caveats</a></p></header>
<div class="filters"><input id="q" aria-label="Search genes or evidence" placeholder="Search gene, tissue, paper ID or evidence"><select id="category" aria-label="Current category"><option value="">All current categories</option><option>CTA</option><option>Neither</option><option>Unsafe</option></select><select id="scope" aria-label="Review subset"><option value="">All evidence levels</option><option value="changes">Proposed changes</option><option value="manual">Specific review notes</option><option value="antibody">Primary antibody conflicts</option><option value="boundary">Threshold boundaries</option></select></div>
<p id="count" class="muted"></p><section id="rows"></section><button id="more">Show 50 more</button>
<footer>Raw abstracts remain in the local cache. Titles and search hits are leads: verify the gene, paralog, species, isoform and experimental context before using a paper as evidence. No-hit searches do not establish absence or safety.</footer>
<script id="data" type="application/json">__DATA__</script><script>
const rows=JSON.parse(document.getElementById('data').textContent);let limit=50;
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function links(s){return String(s).split(' | ').filter(Boolean).map(u=>/^https:\/\//.test(u)?`<a href="${esc(u)}" target="_blank" rel="noopener">${esc(u)}</a>`:esc(u)).join('<br>')}
function pairs(r,keys){return '<dl>'+keys.filter(k=>r[k]!==''&&r[k]!=null).map(k=>'<dt>'+esc(k)+'</dt><dd>'+(k.endsWith('sources')||k.endsWith('source')?links(r[k]):esc(r[k]))+'</dd>').join('')+'</dl>'}
function draw(){const q=document.getElementById('q').value.toLowerCase(),c=document.getElementById('category').value,s=document.getElementById('scope').value;
let selected=rows.filter(r=>(!c||r.Category===c)&&(!q||JSON.stringify(r).toLowerCase().includes(q))&&(!s||s==='changes'&&r.Category!==r['Proposed category']||s==='manual'&&r['New finding']||s==='antibody'&&r['Primary normal staining rows']>0||s==='boundary'&&r['Boundary review']));
document.getElementById('count').textContent=`${selected.length.toLocaleString()} matching genes · showing ${Math.min(limit,selected.length)}`;
document.getElementById('rows').innerHTML=selected.slice(0,limit).map(r=>`<article><h2>${esc(r.Gene||r['Gene ID'])}</h2><span class="tag">${esc(r.Category)} · ${esc(r.Detail)}</span><span class="tag ${r.Category!==r['Proposed category']?'change':''}">${esc(r.Assessment)}</span><p>${esc(r['New finding']||r.Implication)}</p>${pairs(r,['Proposed category','Assessment basis','Reproductive RNA (nTPM)','Normal RNA maximum (nTPM)','Normal RNA maximum (TPM)','Highest normal nTPM context','Highest normal TPM context','New review sources'])}<details><summary>Reason, next step and evidence limitations</summary>${pairs(r,['Reason','Implication','Next action','Search coverage','Truncated searches','Primary normal staining rows','Primary vital Medium/High rows','Off-target warning in profile','Primary staining contexts','Literature conflict','Literature conflict source','Boundary review','Immune RNA findings','Cell RNA findings','Antibody findings'])}</details><details><summary>Biology, cancer occurrence and prior clinical review</summary>${pairs(r,['Curated tissue biology','Biological function','Tissue PubMed IDs','Cancer RNA at >=5 pTPM','Cancer IHC Medium/High','Cancer RNA cohorts','Cancer IHC groups','Reviewed literature','Reviewed literature sources'])}</details></article>`).join('');
document.getElementById('more').hidden=selected.length<=limit;}
['q','category','scope'].forEach(id=>document.getElementById(id).addEventListener('input',()=>{limit=50;draw()}));document.getElementById('more').onclick=()=>{limit+=50;draw()};draw();</script></main></html>"""
    path.write_text(
        page.replace("__CHANGES__", str(summary["proposed_changes"]))
        .replace("__SEARCHED__", f"{summary['genes_with_both_searches']:,}")
        .replace("__DATA__", payload)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT)
    parser.add_argument(
        "--atlas", type=Path, default=ROOT / "docs/audits/cta-atlas-20261024/all-genes.csv"
    )
    parser.add_argument("--cache", type=Path, default=ROOT / "tmp/cta-literature-20261025/raw")
    parser.add_argument(
        "--allow-incomplete", action="store_true", help="Build a visibly incomplete progress report"
    )
    args = parser.parse_args()
    atlas = read_table(args.atlas)
    ledger = read_table(args.out / "searches.csv")
    decisions = load_adjudications(args.out / "adjudications.yml", set(atlas.Gene))
    source = ROOT / "docs/audits/cta-normal-cells-20261023/primary-antibodies.csv.gz"
    selected, signals = primary_signals(atlas, source)
    result = build_rows(atlas, ledger, decisions, signals)
    complete = int(result["Search coverage"].eq("Both searches retrieved").sum())
    if not args.allow_incomplete and (complete != 2540 or result["Failed searches"].sum()):
        raise ValueError(
            "Search retrieval is incomplete; use --allow-incomplete only for an explicitly provisional report"
        )
    papers = search_evidence(ledger, args.cache)
    reviewed_pmids = defaultdict(set)
    for gene, decision in decisions.items():
        reviewed_pmids[gene].update(map(str, decision["pmids"]))
    for index, paper in papers.iterrows():
        if paper.PMID in reviewed_pmids[paper.Gene]:
            papers.at[index, "Review"] = (
                "Discussed in explicit gene adjudication; see its review depth and limitations"
            )
    summary = {
        "date": "2026-10-25",
        "genes": len(result),
        "current_categories": result.Category.value_counts().to_dict(),
        "genes_with_both_searches": complete,
        "searches": len(ledger),
        "failed_searches": int(result["Failed searches"].sum()),
        "gene_paper_leads": len(papers),
        "unique_paper_leads": papers.PMID.nunique(),
        "truncated_searches": int(result["Truncated searches"].sum()),
        "genes_with_no_hits_in_either_topic_query": int(
            result["Zero-hit topic searches"].eq(2).sum()
        ),
        "genes_with_specific_assessment_notes": len(decisions),
        "genes_with_new_primary_review": sum(
            d["depth"].startswith("Primary") for d in decisions.values()
        ),
        "genes_without_new_individual_primary_review": len(result)
        - sum(d["depth"].startswith("Primary") for d in decisions.values()),
        "proposed_changes": int(result["Change proposed"].sum()),
        "proposed_categories_if_accepted": result["Proposed category"].value_counts().to_dict(),
        "retained_primary_normal_signal_genes": len(signals),
        "retained_primary_normal_signal_rows": len(selected),
        "retained_vital_medium_high_signal_genes": int(
            result["Primary vital Medium/High rows"].gt(0).sum()
        ),
        "retained_signal_genes_with_off_target_profile_notes": int(
            result["Off-target warning in profile"].sum()
        ),
        "status": "Search retrieval and source-flag audit complete; individual primary-paper adjudication remains incomplete"
        if complete == 2540
        else "Search retrieval in progress; individual primary-paper adjudication incomplete",
        "production_categories_changed": False,
        "inputs": {
            str(p.relative_to(ROOT)): digest(p)
            for p in [args.atlas, source, args.out / "adjudications.yml", args.out / "searches.csv"]
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    for name, frame in [
        ("all-genes.csv", result),
        ("papers.csv", papers),
        ("proposed-changes.csv", result.loc[result["Change proposed"]]),
        ("retained-primary-normal-signals.csv", selected),
    ]:
        frame.to_csv(args.out / name, index=False)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_explorer(result, summary, args.out / "explore.html")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
