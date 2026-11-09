#!/usr/bin/env python3
"""Create the filterable rebuilt catalog, migrations and count figure."""

from __future__ import annotations

import hashlib
import html
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_literature, target_selection  # noqa: E402

OUT = ROOT / "docs/audits/target-selection-20261109/rebuild"


def build():
    data = target_selection.results()
    meta = target_selection.metadata()
    summary = meta["summary"]
    sources = {}
    for row in cta_literature.attestations().to_dict("records"):
        sources.setdefault(
            row["source_id"], {k: row[k] for k in ["source_url", "citation", "title"]}
        )
    columns = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "representative_gene_id",
        "category",
        "classification",
        "released_category",
        "decision_reason",
        "normal_reference_tpm",
        "normal_reference_contexts",
        "tumor_max_tpm",
        "tumor_top_sample",
        "tumor_top_cancer_code",
        "independent_patients_above_normal",
        "unresolved_patient_samples_above_normal",
        "cancer_evidence",
        "cancer_routes",
        "reproductive_routes",
        "reproductive_ntpm",
        "reproductive_cell_ntpm",
        "developmental_mean_tpm",
        "published_nomination",
        "nomination_source_ids",
        "paper_source_ids",
        "representative_transcripts",
    ]
    rows = data[columns].astype(object).where(data[columns].notna(), None).to_dict("records")
    for row in rows:
        for key in [
            "member_gene_ids",
            "nomination_source_ids",
            "paper_source_ids",
            "normal_reference_contexts",
        ]:
            row[key] = json.loads(row[key]) if row[key] else []
    payload = json.dumps(rows, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
    source_payload = json.dumps(sources, allow_nan=False).replace("<", "\\u003c")
    counts = summary["categories"]
    count_rows = "".join(
        f"<tr><td>{html.escape(c)}</td><td>{counts.get(c, 0):,}</td><td>{summary['published_nomination_by_category'].get(c, 0):,}</td></tr>"
        for c in target_selection.CATEGORIES
    )
    migrations = "".join(
        f"<tr><td>{r['released_category']}</td><td>{r['category']}</td><td>{r['proteins']:,}</td></tr>"
        for r in summary["changes"]
        if r["released_category"] != r["category"]
    )
    options = "".join(f"<option>{html.escape(c)}</option>" for c in target_selection.CATEGORIES)
    downloads = " · ".join(
        f'<a href="{c.lower().replace(" ", "-")}.csv.gz">{html.escape(c)}</a>'
        for c in target_selection.CATEGORIES
    )
    members = summary["memberships"]
    reproductive_only = data[data.category.eq("Reproductive Antigen only")]
    unresolved_reproductive = int(reproductive_only.cancer_evidence.eq("unassessed").sum())
    other_counts = ", ".join(
        f"{summary['classifications'].get('Other:' + name, 0):,} {name.lower()}"
        for name in ["Review", "Sequence", "Expression"]
    )
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Rebuilt antigen catalog — OncoRef</title>
<style>body{font:16px/1.5 system-ui,sans-serif;color:#20333d;margin:0}main{max-width:1550px;margin:auto;padding:28px}
h1{font-size:34px}h2{margin-top:2em}a{color:#126487}.summary{max-width:850px}.note{padding:14px;background:#f4f6f7;border-left:4px solid #b17a31}
table{border-collapse:collapse;width:100%}th,td{padding:9px;text-align:left;border-bottom:1px solid #dbe1e4;vertical-align:top}th{background:#eef3f5;position:sticky;top:0}
.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin:20px 0}input,select,button{font:inherit;padding:7px}input[type=search]{min-width:330px}
.scroll{overflow:auto;max-height:78vh}.catalog{min-width:1250px;font-size:14px}small{color:#576771}details{max-width:480px}details ul{padding-left:18px}
.badge{font-weight:600;white-space:nowrap}.Unsafe{color:#a33532}.CTA{color:#20613c}#count{font-weight:bold}.pager{margin-top:12px;display:flex;gap:16px}
</style><main><h1>Rebuilt antigen catalog</h1>
<p>__TOTAL__ proteins assessed under <code>__POLICY__</code>. Identical selected reference proteins are counted together; every member locus remains listed.</p>
<div class="summary"><table><thead><tr><th>Exclusive category</th><th>Proteins</th><th>With a curated target nomination</th></tr></thead><tbody>__COUNTS__</tbody></table>
<p><strong>Cancer Antigens: __CA__ · Reproductive Antigens: __RA__ · CTA intersection: __CTA__ · Target candidates union: __UNION__.</strong></p></div>
<p>Other: __OTHER_COUNTS__. Of the __RA_ONLY__ Reproductive Antigens only, __RA_UNRESOLVED__ have potentially qualifying repeated tumor measurements with unresolved patient identities. The remaining __RA_INSUFFICIENT__ do not meet either cancer RNA rule in the measured samples.</p>
<p>The shared eligibility screen leaves __ELIGIBLE__ proteins. Cancer RNA needs two confirmed independent patients each ≥1 TPM and above the eligible normal maximum N, or one sample ≥5 × max(N, 1). Reproductive/developmental RNA association uses ≥1 with the established read/cell support. The reproductive concentration score no longer excludes candidates.</p>
<p class="note">The comparison uses native gene TPM with resolved identities and source quality checks. N is the maximum supported tissue/cell summary, not the largest value in every healthy donor. Cross-study batch effects are not corrected. Other assay scales retain their independent normal-exclusion rules. <a href="comparison.json">Exact comparison scope</a> · <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5903355/">Cross-study methods</a>.</p>
<p><strong>NUTM coverage:</strong> the stored UNC three-sample matrix fails the nonnegative-TPM check and contributes no qualifying samples. Its failed source remains in the <a href="tumor-source-coverage.csv">coverage ledger</a>; <a href="https://github.com/pirl-unc/oncoref/issues/604">#604</a> and <a href="https://github.com/pirl-unc/oncoref/issues/665">#665</a> track the repair. No raw reads were re-quantified for this rebuild.</p>
<h2>Changes from the released catalog</h2><div class="summary"><table><thead><tr><th>Released category</th><th>New category</th><th>Proteins</th></tr></thead><tbody>__MIGRATIONS__</tbody></table></div>
<p><a href="category-changes.csv">All changed proteins and their evidence</a> · <a href="summary.json">Counts</a> · <a href="gene-decisions.csv.gz">All member genes</a> · <a href="protein-decisions.csv.gz">Full protein decisions</a></p>
<p>Category downloads: __DOWNLOADS__.</p>
<h2>Find a protein</h2><p>“Only” describes qualified set membership. Missing evidence is shown as —, never zero. The cancer values are recalculated for the __ELIGIBLE__ proteins passing shared eligibility; other rows retain their exclusion or review decision.</p>
<div class="controls"><input id="search" type="search" placeholder="Protein, Ensembl ID, reason, paper or program" aria-label="Search proteins">
<select id="category" aria-label="Category"><option value="">All categories</option>__OPTIONS__</select>
<label><input type="checkbox" id="changed"> Changed category</label><label><input type="checkbox" id="published"> Published target nomination</label>
<span id="count"></span></div>
<div class="scroll"><table class="catalog"><thead><tr><th>Protein and source links</th><th>New category</th><th>Released</th><th>Normal maximum<br>TPM</th><th>Tumor maximum<br>TPM</th><th>Confirmed patients<br>above normal</th><th>Reproductive bulk<br>nTPM</th><th>Reproductive cells<br>nTPM</th><th>Developmental mean<br>TPM</th><th>Evidence and papers</th></tr></thead><tbody id="body"></tbody></table></div>
<div class="pager"><button id="previous">Previous</button><span id="page"></span><button id="next">Next</button></div>
<p><a href="normal-comparison.csv.gz">All comparison normal measurements</a> · <a href="tumor-witnesses.csv.gz">Tumor witnesses</a> · <a href="paper-attestations.json.gz">Every retained protein–paper claim</a> · <a href="inputs.json">Input hashes</a> · <a href="../report.html">Ontology figures and rules</a></p>
<script id="records" type="application/json">__DATA__</script><script id="sources" type="application/json">__SOURCES__</script>
<script>
const records=JSON.parse(document.getElementById('records').textContent),sources=JSON.parse(document.getElementById('sources').textContent);
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=x=>x===null?'—':Number(x).toLocaleString(undefined,{maximumFractionDigits:3});
const link=(url,label)=>`<a href="${esc(url)}" target="_blank" rel="noopener">${esc(label)}</a>`;
let page=0,selected=[];const size=100;
function papers(ids){return ids.map(id=>{let s=sources[id];return s?link(s.source_url,s.citation)+(s.title?' — '+esc(s.title):''):esc(id)}).join('<br>')}
function render(reset=true){if(reset)page=0;let q=document.getElementById('search').value.toLowerCase(),cat=document.getElementById('category').value;
selected=records.filter(r=>(!cat||r.category===cat)&&(!document.getElementById('changed').checked||r.category!==r.released_category)&&(!document.getElementById('published').checked||r.published_nomination)&&(!q||JSON.stringify(r).toLowerCase().includes(q)||r.paper_source_ids.some(id=>JSON.stringify(sources[id]||{}).toLowerCase().includes(q))));
const pages=Math.max(1,Math.ceil(selected.length/size));page=Math.min(page,pages-1);
document.getElementById('count').textContent=selected.length.toLocaleString()+' proteins';
document.getElementById('body').innerHTML=selected.slice(page*size,(page+1)*size).map(r=>`<tr><td><strong>${esc(r.protein_group_label)}</strong><br><small>${r.member_gene_ids.map(g=>link('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+g,g)).join('<br>')}<br>${link('https://www.proteinatlas.org/'+r.representative_gene_id,'HPA')}</small></td>
<td class="badge ${r.category==='Unsafe'?'Unsafe':r.category==='CTA'?'CTA':''}">${esc(r.classification)}</td><td>${esc(r.released_category)}</td><td>${fmt(r.normal_reference_tpm)}</td><td>${fmt(r.tumor_max_tpm)}<br><small>${esc(r.tumor_top_cancer_code)}</small></td><td>${fmt(r.independent_patients_above_normal)}</td><td>${fmt(r.reproductive_ntpm)}</td><td>${fmt(r.reproductive_cell_ntpm)}</td><td>${fmt(r.developmental_mean_tpm)}</td>
<td><details><summary>${r.published_nomination?'Published nomination':'No catalog nomination'} · ${r.paper_source_ids.length} paper sources</summary><p>${esc(r.decision_reason)}</p><p>Cancer assessment: ${esc(r.cancer_evidence)}<br>${esc(r.cancer_routes)}<br>Samples above normal with unresolved patient IDs: ${fmt(r.unresolved_patient_samples_above_normal)}</p><p>Reproductive/developmental routes: ${esc(r.reproductive_routes)}</p><p>Peak tumor sample: ${esc(r.tumor_top_sample)}</p><p>Selected transcript(s): ${esc(r.representative_transcripts)}</p><p>Normal maximum witnesses (up to 10):</p><ul>${r.normal_reference_contexts.map(x=>`<li>${esc(x.tissue)} / ${esc(x.cell_type)}: ${fmt(x.value)} TPM · ${esc(x.source_id)} · ${esc(x.statistic)}</li>`).join('')}</ul><p><strong>Target nominations</strong><br>${papers(r.nomination_source_ids)||'None recorded'}</p><p><strong>All paper roles, including counterevidence</strong><br>${papers(r.paper_source_ids)||'None recorded'}</p></details></td></tr>`).join('');
document.getElementById('page').textContent=`Page ${page+1} of ${pages}`;document.getElementById('previous').disabled=page===0;document.getElementById('next').disabled=page+1>=pages;}
for(const id of ['search','category','changed','published'])document.getElementById(id).addEventListener(id==='search'?'input':'change',()=>render());
document.getElementById('previous').onclick=()=>{page--;render(false)};document.getElementById('next').onclick=()=>{page++;render(false)};render();
</script></main></html>"""
    replace = {
        "TOTAL": f"{len(data):,}",
        "POLICY": target_selection.POLICY_ID,
        "COUNTS": count_rows,
        "CA": f"{members['is_cancer_antigen']:,}",
        "RA": f"{members['is_reproductive_antigen']:,}",
        "CTA": f"{members['is_cta']:,}",
        "UNION": f"{members['is_target_candidate']:,}",
        "ELIGIBLE": str(summary["shared_eligible"]),
        "OTHER_COUNTS": other_counts,
        "RA_ONLY": str(len(reproductive_only)),
        "RA_UNRESOLVED": str(unresolved_reproductive),
        "RA_INSUFFICIENT": str(len(reproductive_only) - unresolved_reproductive),
        "MIGRATIONS": migrations,
        "DOWNLOADS": downloads,
        "OPTIONS": options,
        "DATA": payload,
        "SOURCES": source_payload,
    }
    for key, value in replace.items():
        page = page.replace("__" + key + "__", value)
    (OUT / "report.html").write_text(page)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.2, 1]})
    labels = ["CTA", "Cancer Antigen only", "Reproductive Antigen only"]
    axes[0].barh(
        labels[::-1], [counts[x] for x in labels[::-1]], color=["#c99a46", "#3e8bba", "#438a69"]
    )
    for i, label in enumerate(labels[::-1]):
        axes[0].text(counts[label] + 5, i, f"{counts[label]:,}", va="center")
    axes[0].set_xlim(0, max(counts[x] for x in labels) * 1.17)
    axes[0].set_title(f"{members['is_target_candidate']:,} Target candidates", loc="left")
    axes[1].barh(
        ["Other", "Unsafe"], [counts["Other"], counts["Unsafe"]], color=["#9aabb5", "#ac5552"]
    )
    for i, label in enumerate(["Other", "Unsafe"]):
        axes[1].text(counts[label] + 130, i, f"{counts[label]:,}", va="center")
    axes[1].set_xlim(0, counts["Unsafe"] * 1.22)
    axes[1].set_title("Remaining inventory", loc="left")
    for ax in axes:
        ax.set_xlabel("Proteins")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "counts.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    new_cta = data[data.category.eq("CTA") & data.released_category.ne("CTA")]
    held_cta = data[data.category.eq("Other") & data.released_category.eq("CTA")]
    migration_table = "\n".join(
        f"| {r['released_category']} | {r['category']} | {r['proteins']} |"
        for r in summary["changes"]
    )
    table_rows = "\n".join(
        f"| {r.protein_group_label} | {r.tumor_max_tpm:.3g} | {r.normal_reference_tpm:.3g} | "
        f"{'Patient identity unresolved' if r.cancer_evidence == 'unassessed' else 'Below cancer rule'} |"
        for r in reproductive_only.itertuples()
    )
    (OUT / "index.md").write_text(
        "# Rebuilt antigen catalog\n\n"
        f"**{counts['CTA']} CTA, {counts['Cancer Antigen only']} Cancer Antigen only, "
        f"{counts['Reproductive Antigen only']} Reproductive Antigen only, {counts['Other']} Other "
        f"and {counts['Unsafe']:,} Unsafe.**\n\n"
        "[Filterable protein table](report.html) · [All category changes](category-changes.csv) · "
        "[Exact comparison scope](comparison.json)\n\n![Rebuilt counts](counts.png)\n\n"
        f"Cancer Antigens include {members['is_cancer_antigen']} proteins; Reproductive Antigens include "
        f"{members['is_reproductive_antigen']}. Their intersection is CTA and their union contains "
        f"{members['is_target_candidate']} Target candidates.\n\n"
        f"Other contains {other_counts}. A curated published target nomination is retained for "
        f"{summary['published_nomination_by_category']['CTA']} CTAs and "
        f"{summary['published_nomination_by_category']['Cancer Antigen only']} Cancer Antigens only. "
        "Broader paper coverage also includes biology, normal-expression findings and other evidence.\n\n"
        "## Category changes\n\n| Released category | Rebuilt category | Proteins |\n|---|---|---:|\n"
        + migration_table
        + "\n\nNew CTAs: "
        + ", ".join(new_cta.protein_group_label)
        + ".\n\nFormer CTAs now requiring review: "
        + ", ".join(held_cta.protein_group_label)
        + ".\n\n## Reproductive Antigens only\n\n"
        + f"{unresolved_reproductive} have potentially qualifying tumor measurements but unresolved patient "
        "independence; none reaches the single-sample threshold. The remaining "
        f"{len(reproductive_only) - unresolved_reproductive} do not meet either cancer RNA rule.\n\n"
        "| Protein | Tumor maximum TPM | Normal comparison TPM | Cancer assessment |\n"
        "|---|---:|---:|---|\n" + table_rows + "\n\n## Evidence and reproduction\n\n"
        "The comparison uses source-native TPM and maxima of supported tissue/cell summaries. "
        "It does not correct cross-study batch effects. NUTM's stored three-sample matrix fails "
        "the nonnegative-TPM check and is listed in the [coverage ledger](tumor-source-coverage.csv).\n\n"
        "[Full protein decisions](protein-decisions.csv.gz) · [All member genes](gene-decisions.csv.gz) · "
        "[Normal comparison values](normal-comparison.csv.gz) · [Tumor witnesses](tumor-witnesses.csv.gz) · "
        "[Paper claims](paper-attestations.json.gz) · [Input hashes](inputs.json).\n\n"
        "With the pinned source files in the paths listed in `inputs.json`, run:\n\n"
        "```sh\npython scripts/rebuild_target_catalog.py\npython scripts/render_target_catalog.py\n"
        "python scripts/build_target_selection_report.py\n```\n\n"
        "The builder verifies the pinned tumor matrices, re-resolves source identities and assesses "
        "all proteins. Temporary tumor matrices stay outside the published report. "
        "The new catalog is available through `target_selection.results()`; the `cta_selection` "
        "API and existing sequence downloads retain the 1.36.0 snapshot.\n"
    )
    implementation = [
        "scripts/rebuild_target_catalog.py",
        "scripts/render_target_catalog.py",
        "scripts/build_cta_candidate_evidence.py",
        "scripts/build_cta_alternative_evidence.py",
        "scripts/cta_report_common.py",
        "oncoref/target_selection.py",
        "oncoref/gene_matching.py",
        "oncoref/cta_literature.py",
        "oncoref/tissue_hazards.py",
    ]
    manifest = {
        "policy_id": target_selection.POLICY_ID,
        "implementation_sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in implementation
        },
        "outputs_sha256": {
            str(path.relative_to(OUT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(OUT.rglob("*"))
            if path.is_file() and path.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    build()
