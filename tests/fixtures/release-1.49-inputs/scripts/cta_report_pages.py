"""Readable report pages and a self-contained, searchable gene explorer."""

from __future__ import annotations

import html
import json

import pandas as pd
from cta_recalibration_report import markdown_table

from oncoref import figure_style


def render(
    out,
    current,
    intake,
    clinical,
    genes,
    funnel,
    evidence_labels,
    result_labels,
    policy,
    sources,
    comparison,
):
    figure_style.apply()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.barh(funnel.Step, funnel.Genes, color="#27638b")
    ax.invert_yaxis()
    for i, value in enumerate(funnel.Genes):
        ax.text(value + 30, i, f"{value:,}", va="center")
    ax.set_xlim(0, funnel.Genes.max() * 1.12)
    ax.set_xlabel("Distinct genes / source identities")
    ax.set_title("From source nominations to included candidates")
    fig.subplots_adjust(left=0.43, bottom=0.13, right=0.98, top=0.90)
    for suffix in ["png", "svg"]:
        figure_style.save(fig, out / f"filtering-funnel.{suffix}", keep_titles=True)
    plt.close(fig)
    retained = int(current.retained.sum())
    gains = int((current.retained & ~current.previous_retained).sum())
    losses = int((~current.retained & current.previous_retained).sum())
    clinical_display = clinical[
        ["Symbol", "Clinical evidence", "Result", "Reason", "source_url", "limitations"]
    ].rename(
        columns={
            "Symbol": "Gene",
            "source_url": "Source",
            "limitations": "What this evidence does not establish",
        }
    )
    clinical_display.to_csv(out / "clinical-results.csv", index=False)
    short = clinical_display[["Gene", "Clinical evidence", "Result"]].sort_values("Gene")
    counts = (
        genes.groupby(["category", "other_subgroup"])
        .size()
        .reset_index(name="Genes")
        .rename(columns={"category": "Group", "other_subgroup": "Other subgroup"})
    )
    rules = pd.DataFrame(
        [
            {
                "Tissue group": group.title(),
                **{
                    label: policy["negative_ihc_rna_limits"][level][group]
                    for label, level in [
                        ("No usable IHC", "Missing"),
                        ("Approved negative", "Approved"),
                        ("Supported negative", "Supported"),
                        ("Enhanced negative", "Enhanced"),
                    ]
                },
            }
            for group in ["vital", "protected", "manageable"]
        ]
    )
    clinical_html = "".join(
        "<tr>" + "".join("<td>" + html.escape(str(value)) + "</td>" for value in row) + "</tr>"
        for row in short.itertuples(index=False, name=None)
    )
    candidates = current.set_index("Ensembl_Gene_ID")
    trial = clinical.set_index("Ensembl_Gene_ID")
    records = []
    for row in intake.itertuples():
        gid = row.Ensembl_Gene_ID
        if gid in candidates.index:
            item = candidates.loc[gid]
            source = trial.loc[gid] if gid in trial.index else None
            records.append(
                {
                    "Gene": item.Symbol,
                    "Gene ID": gid,
                    "Group": item.category,
                    "Result": item.Result,
                    "Reason": item.Reason or "All current screening and review gates pass",
                    "RNA findings": item.rna_exclusions,
                    "Protein findings": item.ihc_exclusions,
                    "Selectivity findings": item.selectivity_exclusions,
                    "Highest reproductive RNA (nTPM)": item.reproductive_rna_max_ntpm,
                    "Publication DOIs": item.source_dois,
                    "Clinical evidence": item["Clinical evidence"],
                    "Also known as": source.source_target_label if source is not None else "",
                    "Clinical source": source.source_url if source is not None else "",
                    "Evidence limits": source.limitations if source is not None else "",
                    "Missing RNA contexts": int(item.rna_routine_missing),
                    "Missing IHC contexts": int(item.ihc_routine_missing),
                    "Change": item.Change,
                }
            )
        else:
            mapped = bool(gid)
            records.append(
                {
                    "Gene": row.Symbol,
                    "Gene ID": gid or row.intake_identity,
                    "Group": "Unassessed",
                    "Result": "Not protein-coding" if mapped else "Missing gene mapping",
                    "Reason": "Outside the mapped protein-coding candidate universe"
                    if mapped
                    else "A source nomination exists, but no unique canonical gene mapping is available",
                    "Clinical evidence": "",
                }
            )
    overview = f"**{len(intake):,} source identities → {len(current):,} protein-coding candidates → {retained:,} included.**"
    source_summary = sources[
        ["Publication / list", "Source rows", "Mapped coding candidates", "Included"]
    ]
    text = f"""# CTA curation and filtering

{overview} Updated 21 October 2026, with a minimum outside-tissue RNA ceiling of **2 nTPM**.

## Start here

- **[Find a gene and read its evidence](explore.html#genes)** — search names, aliases or IDs; filter by result and expression group.
- **[Clinical targets](#clinical-targets)** — readable results for antigens with clinical evidence.
- **[Publication sources](#publication-sources)** — what each source contributes to the current set.
- **[Filtering rules](#filtering-rules)** — the RNA, protein-staining and selectivity rules.
- **[Gene and protein sets](../normal-sequence-sets-20261021/index.md)** — disjoint lists and FASTA files for k-mer comparisons.
- **[Downloads](#downloads)** — tables, figures and reproducibility records.

## What changed

**{gains:,} newly included; {losses:,} no longer included**, compared with the previous [648-gene report](../cta-curation-20261021/index.md).

**XAGE1A and XAGE1B are both included.** XAGE1B's liver 1.5 and lung 1.2 nTPM measurements now pass the 2 nTPM baseline through the general rule. Shared protein sequences are resolved separately in the gene/protein sets.

All absolute RNA ceilings were doubled. Reproductive admission remains RNA >0 or positive IHC. Selectivity percentages and the exclusion for any positive normal-tissue IHC remain unchanged. The [complete change list](changes.csv) records every affected gene.

## Filtering funnel

{markdown_table(funnel)}

![Filtering funnel](filtering-funnel.png)

Each row applies an additional rule to the preceding row. Individual genes can fail several rules; the downloadable findings preserve all of them.

## Compare the screens

{markdown_table(comparison)}

**Fixed** uses the baseline ceilings (Vital 2, Protected 6, Manageable 10 nTPM) and baseline selectivity percentages for every tissue. **IHC-adaptive** is the current result: complete local negative IHC can relax both thresholds. **Zero measured outside RNA** adds a stricter requirement of zero measured RNA in all screened normal tissues. All three keep the same positive-IHC veto, reproductive admission, evidence-coverage requirements and curation reviews. Zero measured RNA does not establish protein absence or coverage of unmeasured tissues.

[Compare individual genes](threshold-comparison.csv) · [Every tissue's RNA, reference and IHC support](tissue-evidence.csv.gz).

## Publication sources

{markdown_table(source_summary)}

Sources overlap, so these rows must not be added together. Counts refer to distinct mapped coding candidates within each source; published source rows can include historical noncoding or unmapped entries. The clinical evidence audit adds targeted primary evidence and nominated comparator genes independently of these publication lists.

[Download exact citations, source tables and current counts](publication-results.csv).

## Clinical targets

{markdown_table(short)}

**Clinical evidence describes the source, not a safety grade.** A trial registration names an intended treatment target; it does not establish that treatment occurred or worked. A published clinical study reports patient experience. An institutional report is a primary treatment account without necessarily being a peer-reviewed trial paper. An approved therapy has a product- and indication-specific approval. Shared targets and antigen-family labels do not resolve every gene or isoform.

**Result describes this screen.** Included means the current rules pass. Normal protein detected means IHC staining was observed outside Reproductive/Disposable. RNA above limit and Selectivity below limit identify numerical failures requiring review. Safety hold records external evidence that prevents automatic inclusion. Needs review covers evidence or curation gaps. None of these labels alone predicts clinical outcomes.

### Important exclusions

- **MAGEA12:** an explicit hold based on published normal-brain expression and the anti-MAGE trial's neurological injuries, despite an HPA screen pass. Unique causal attribution is not established. [Primary report](https://pubmed.ncbi.nlm.nih.gov/23377668/)
- **SPAG9 / SP17:** measured normal-tissue protein staining remains an exclusion. Trial use does not erase that evidence.
- **LY6K:** remaining RNA failures include midbrain 2.4, hypothalamus/retina 2.3 and basal ganglia 2.2 nTPM against the 2 nTPM baseline.
- **CLDN6:** the adult reproductive reference is weak and fetal tissue coverage is missing. Its selectivity failure is a reference-coverage concern, not demonstrated clinical toxicity.
- **MAGEA3:** included, with a separate receptor-specific titin cross-reactivity review. That clinical event does not establish cardiac MAGEA3 expression. [Primary report](https://pmc.ncbi.nlm.nih.gov/articles/PMC3743463/)

[Download the clinical results with sources and evidence limits](clinical-results.csv). The 30 reviewed loci include related oncofetal targets and safety comparators; this is not a claim of 30 independently approved therapies or exhaustive worldwide trial coverage.

## Filtering rules

**Maximum RNA, in HPA v23 nTPM:**

{markdown_table(rules)}

IHC means antibody staining in normal tissue. A negative result means “Not detected”; Approved, Supported and Enhanced are HPA reliability labels. A higher ceiling requires negative results across every routine cell type measured in that same tissue, using the weakest reliability. Missing or uncertain IHC uses the baseline column. **Any positive IHC, including Low, excludes in Vital, Protected and Manageable tissues.**

Selectivity = highest reproductive RNA / (highest reproductive RNA + that tissue's RNA). Minimum percentages, in the same column order: **Vital 75/70/65/60%; Protected 70/65/60/55%; Manageable 65/60/55/50%.** Both RNA and selectivity must pass for every screened tissue. Disposable tissue contributes neither positive reproductive admission nor an exclusion.

These are project research-screening parameters. [Full tissue definitions and citations](../../normal-tissue-scopes.md) · [Exact policy](screening-policy.json).

## BAGE and GAGE2B: what is missing?

These are real, published genes. The problem is reference mapping and gene-specific expression coverage.

| Gene | What the current Ensembl check found | Consequence |
|---|---|---|
| BAGE | No current symbol lookup or symbol cross-reference in Ensembl 116. HGNC has no Ensembl ID in its retrieved record. | Keep its NCBI/HGNC identity in intake; HPA v23 does not supply a usable BAGE-specific RNA/IHC record. |
| GAGE2B | Historical ENSG00000231850 ended at GRCh37/release 75. A current symbol cross-reference points to **GAGE2A, ENSG00000189064**, whose record links several GAGE2 genes. | A shared cross-reference is not a unique GAGE2B locus or gene-specific expression measurement. Keep GAGE2B unresolved rather than substitute GAGE2A measurements. |

[Official Ensembl checks and responses](ensembl-identity-checks.json) · [GAGE2B ID history](https://rest.ensembl.org/archive/id/ENSG00000231850?content-type=application/json) · [Current cross-reference](https://rest.ensembl.org/xrefs/symbol/homo_sapiens/GAGE2B?content-type=application/json) · [GAGE2A record](https://rest.ensembl.org/lookup/id/ENSG00000189064?content-type=application/json).

## Reproductive, Normal and Other

{markdown_table(counts)}

These groups cover all 20,097 canonical coding genes, including genes without a CTA publication. Other contains Disposable and Silent under the approved below-threshold definition. Missing routine evidence remains Unassessed. [Expression-group definitions](../normal-expression-categories-20261021-r2/index.md).

For sequence work, [the disjoint gene and protein sets](../normal-sequence-sets-20261021/index.md) account for exact shared proteins and retain unresolved translations. Distinct proteins can still share k-mers; sequence-level subtraction is still required.

## Downloads

| Purpose | Download |
|---|---|
| Read results for every mapped candidate | [Readable gene results](gene-results.csv.gz) |
| Inspect clinical evidence and limitations | [Clinical results](clinical-results.csv) |
| See gains and losses | [Changes](changes.csv) |
| Compare fixed, adaptive and zero-outside-RNA screens | [Per-gene comparison](threshold-comparison.csv) |
| Check the original tissue measurements and IHC confidence | [Tissue evidence](tissue-evidence.csv.gz) |
| See each publication's contribution | [Publication results and citations](publication-results.csv) |
| Reuse gene / protein / sequence sets | [Disjoint sets and FASTA files](../normal-sequence-sets-20261021/index.md) |
| Export the funnel figure | [PNG](filtering-funnel.png) / [SVG](filtering-funnel.svg) / [Counts](funnel.csv) |
| Inspect source nominations, including unresolved IDs | [Complete intake](intake-decisions.csv.gz) |
| Reproduce detailed decisions | [Candidate data](candidate-decisions.csv.gz) / [Clinical audit](clinical-target-audit.csv) / [Checksums](run-manifest.json) |

## Sources and earlier reports

The intake combines the {len(sources)} publication lists or targeted source records above with reviewed primary clinical evidence. [Earlier source-import audit](../cta-recalibration-20261021/index.md) · [Clinical registry snapshot and source evidence](../cta-curation-20261021/index.md#downloads-and-reproducibility) · [Earlier threshold comparison](../cta-recalibration-20261021/index.md) · [Human/dog/mouse tissue mappings](../normal-tissue-registry-20261021-r2/index.md). Earlier numerical reports remain dated archives.
"""
    (out / "index.md").write_text(text)
    page = (
        PAGE.replace("__INCLUDED__", f"{retained:,}")
        .replace("__CANDIDATES__", f"{len(current):,}")
        .replace("__GAINS__", f"{gains:,}")
    )
    page = page.replace("__CLINICAL__", clinical_html).replace(
        "__RULES__", rules.to_html(index=False, border=0)
    )
    page = page.replace("__SOURCES__", source_summary.to_html(index=False, border=0)).replace(
        "__COMPARISON__", comparison.to_html(index=False, border=0)
    )
    page = page.replace(
        "__DATA__", json.dumps(records, ensure_ascii=True, allow_nan=False).replace("<", "\\u003c")
    )
    (out / "explore.html").write_text(page)


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CTA curation: find a gene</title><style>
:root{--ink:#173b5e;--text:#25384a;--blue:#27638b;--pale:#edf3f7;--line:#d1dce5;--green:#21665d;--rust:#9d432b}*{box-sizing:border-box}body{margin:0;background:white;color:var(--text);font:15px/1.65 Verdana,sans-serif}a{color:var(--blue)}a:focus-visible,button:focus-visible,input:focus-visible,select:focus-visible,summary:focus-visible{outline:3px solid #e49a3b;outline-offset:3px}aside{position:fixed;inset:0 auto 0 0;width:225px;background:var(--pale);padding:32px 24px;border-right:1px solid var(--line)}aside strong{font:24px Georgia,serif;color:var(--ink)}nav a{display:block;padding:9px 0;text-decoration:none}nav{margin-top:28px}main{margin-left:225px;max-width:1300px;padding:42px 48px 80px}h1,h2{font-family:Georgia,serif;color:var(--ink);line-height:1.2}h1{font-size:40px;margin:0 0 16px}h2{font-size:28px;margin:0 0 16px}h3{font-size:17px;margin:22px 0 8px}p{max-width:78ch}.intro{font-size:17px;max-width:75ch}.overview{border-left:5px solid var(--green);padding:5px 20px;margin:25px 0;background:#f3f8f6}.overview p{margin:9px 0}.note{font-size:13px;color:#526576}section{padding-top:30px;margin-top:24px;scroll-margin-top:20px;border-top:1px solid var(--line)}.tools{display:grid;grid-template-columns:minmax(190px,2fr) 1fr 1.3fr;gap:14px;background:var(--pale);padding:20px}.tools label{font-size:13px;font-weight:bold}input[type=search],select{display:block;width:100%;padding:10px;margin-top:6px;border:1px solid #889bac;border-radius:3px;background:white;color:var(--text);font:inherit}button{font:inherit;cursor:pointer;border:1px solid #889bac;background:white;border-radius:3px;padding:7px 13px;color:var(--ink)}button:hover{background:var(--pale)}.toolbar{display:flex;flex-wrap:wrap;align-items:center;gap:16px;margin:14px 0}.toolbar #status{flex:1}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th{text-align:left;border-bottom:2px solid #8297a9;color:var(--ink);padding:10px;white-space:nowrap}td{border-bottom:1px solid var(--line);padding:12px 10px;vertical-align:top}td:first-child{font-weight:bold;min-width:120px}tbody tr:nth-child(even){background:#f8fafc}.included{color:var(--green);font-weight:bold}.hold{color:var(--rust);font-weight:bold}details summary{cursor:pointer;color:var(--blue);font-size:12px;font-weight:normal;margin-top:7px}dl{font-weight:normal;min-width:240px;max-width:650px;font-size:12px}dt{font-weight:bold;margin-top:10px}dd{margin:2px 0;overflow-wrap:anywhere}.pager{display:flex;justify-content:flex-end;gap:12px;align-items:center;margin:18px 0}.legend{border-left:3px solid var(--blue);padding-left:18px}.links{display:flex;flex-wrap:wrap;gap:12px}.links a{border:1px solid var(--line);padding:10px 15px;text-decoration:none}footer{margin-top:42px;font-size:12px;color:#526576}button:disabled{opacity:.45;cursor:default}@media(max-width:850px){aside{position:static;width:auto;padding:18px 22px;border-bottom:1px solid var(--line)}nav{margin-top:12px;display:flex;gap:18px;flex-wrap:wrap}nav a{padding:0}main{margin:0;padding:28px 22px}h1{font-size:32px}.tools{grid-template-columns:1fr}.table-wrap table{min-width:670px}}@media print{aside,.tools,.toolbar,.pager{display:none}main{margin:0;padding:0}section{break-inside:avoid}}
</style></head><body><aside><strong>oncoref</strong><nav aria-label="Report sections"><a href="#overview">Overview</a><a href="#genes">Find a gene</a><a href="#clinical">Clinical targets</a><a href="#sources">Publication sources</a><a href="#rules">Filtering rules</a><a href="#identities">BAGE / GAGE2B</a><a href="#sets">Gene and protein sets</a><a href="#downloads">Downloads</a><a href="./">Full written report</a></nav></aside><main>
<header id="overview"><h1>CTA curation</h1><p class="intro">Which candidates pass the normal-tissue screen, what excludes the others, and where the evidence is incomplete.</p><p class="note">21 October 2026 · Human HPA v23 · RNA ceilings start at 2 nTPM</p><div class="overview"><p><strong>__INCLUDED__ included</strong> from __CANDIDATES__ mapped protein-coding candidates. <strong>__GAINS__ newly included</strong> under the revised RNA ceilings.</p><p><strong>XAGE1A and XAGE1B are both included.</strong> XAGE1B's liver 1.5 and lung 1.2 nTPM measurements now pass the general baseline.</p></div><img src="filtering-funnel.svg" alt="Filtering funnel from 3,908 source identities to the included CTA set" style="width:100%;max-width:1000px;height:auto"><p>Use the gene search to inspect results. “Included” means this research screen passes; clinical use and safety evidence are shown separately.</p></header>
<section id="genes"><h2>Find a gene</h2><div class="tools"><label>Gene, alias or ID<input id="search" type="search" placeholder="Try XAGE1B, PRAME or KK-LC-1" autocomplete="off"></label><label>Expression group<select id="group"><option value="">All groups</option><option>Reproductive</option><option>Normal</option><option>Other</option><option>Unassessed</option></select></label><label>Result<select id="result"><option value="">All results</option></select></label></div><div class="toolbar"><label><input type="checkbox" id="clinical-only"> Clinical evidence only</label><button id="clear" type="button">Clear filters</button><span id="status" role="status" aria-live="polite"></span></div><div class="table-wrap"><table><thead><tr><th scope="col">Gene</th><th scope="col">Expression group</th><th scope="col">Result</th><th scope="col">Reason</th><th scope="col">Clinical evidence</th></tr></thead><tbody id="results"></tbody></table></div><div class="pager"><button id="prev" type="button">Previous</button><span id="page"></span><button id="next" type="button">Next</button></div><p class="note">Open “Evidence” beneath a gene for tissue measurements, coverage gaps and primary clinical sources. Unmapped source nominations are included in this search.</p></section>
<section id="sources"><h2>Where the candidates came from</h2><p>Sources overlap. These counts show the current outcome for each source and must not be added together.</p><div class="table-wrap">__SOURCES__</div><p><a href="publication-results.csv">Download citations, exact source tables and counts</a></p><h3>Compare the screens</h3><div class="table-wrap">__COMPARISON__</div><p>Fixed uses the baseline RNA and selectivity thresholds. IHC-adaptive allows the approved relaxation for complete local negative IHC. Zero measured outside RNA additionally requires zero measured RNA in every screened normal tissue. All three retain the other screening and review gates; missing evidence never becomes zero.</p><p><a href="threshold-comparison.csv">Compare individual genes</a> · <a href="tissue-evidence.csv.gz">All tissue measurements and IHC support</a></p></section>
<section id="clinical"><h2>Clinical targets</h2><p>Evidence describes the source. Result describes the current screen. A trial registration records an intended treatment target; it does not prove that patients were treated or that treatment worked.</p><div class="table-wrap"><table><thead><tr><th>Gene</th><th>Clinical evidence</th><th>Result</th></tr></thead><tbody>__CLINICAL__</tbody></table></div><p class="legend"><strong>Safety hold:</strong> external normal-tissue evidence prevents automatic inclusion. <strong>RNA above limit</strong> and <strong>Selectivity below limit:</strong> numerical failures that require review. <strong>Normal protein detected:</strong> positive IHC outside permitted tissues.</p><p><a href="clinical-results.csv">Download clinical results, sources and limitations</a></p></section>
<section id="rules"><h2>Filtering rules</h2><p>Reproductive RNA must be above zero, or reproductive protein staining must be positive. Every screened normal tissue must pass both its absolute RNA ceiling and reproductive selectivity threshold.</p><h3>Maximum RNA (nTPM)</h3><div class="table-wrap">__RULES__</div><p>IHC is antibody staining in normal tissue. Negative means “Not detected.” The confidence labels are HPA reliability categories. Relaxation requires negative results in all routine cell types measured in that tissue. <strong>Any positive IHC, including Low, excludes in Vital, Protected or Manageable tissues.</strong></p><h3>Minimum reproductive selectivity</h3><p>Highest reproductive RNA ÷ (highest reproductive RNA + that tissue's RNA). In the same column order: Vital <strong>75 / 70 / 65 / 60%</strong>; Protected <strong>70 / 65 / 60 / 55%</strong>; Manageable <strong>65 / 60 / 55 / 50%</strong>.</p><p>These are project screening parameters. <a href="../../normal-tissue-scopes/">Read the full tissue definitions and citations.</a></p></section>
<section id="identities"><h2>BAGE and GAGE2B</h2><p>The genes are present in source nominations. The missing piece is a unique usable reference mapping and gene-specific expression data.</p><p><strong>BAGE:</strong> no current Ensembl 116 symbol lookup or cross-reference was returned. Its HGNC/NCBI identity remains recorded.</p><p><strong>GAGE2B:</strong> the historical ENSG00000231850 record ended at GRCh37/release 75. Current cross-references point to GAGE2A's record and link several GAGE2 genes. This is not a unique GAGE2B locus, so GAGE2A expression is not silently substituted.</p><p><a href="ensembl-identity-checks.json">Inspect the official Ensembl responses</a></p></section>
<section id="sets"><h2>Gene and protein sets</h2><p>The three main expression groups are Reproductive, Normal and Other. Other contains Disposable and Silent below-threshold expression. Unassessed remains separate.</p><p>Exact protein sequences are assigned once, with Normal taking precedence over Unassessed, then Reproductive, then Other. All complete isoforms are considered. Distinct proteins can still share k-mers, so the exclusion step must compare their actual substrings.</p><div class="links"><a href="../normal-sequence-sets-20261021/">Disjoint gene lists and FASTA files</a><a href="../normal-expression-categories-20261021-r2/">Expression groups and coverage</a></div></section>
<section id="downloads"><h2>Downloads</h2><div class="links"><a href="gene-results.csv.gz">Readable gene results</a><a href="clinical-results.csv">Clinical results</a><a href="changes.csv">Changes</a><a href="filtering-funnel.svg">Funnel figure</a><a href="run-manifest.json">Checksums</a></div><p><a href="./#sources-and-earlier-reports">Publication sources, earlier reports and complete audit files</a></p></section><footer>Source versions, coverage gaps and decisions are preserved. Protein groups inherit gene-level expression evidence; isoform-specific expression and clinical safety are not established by this screen.</footer>
</main><script type="application/json" id="gene-data">__DATA__</script><script>
'use strict';
const data=JSON.parse(document.getElementById('gene-data').textContent).sort((a,b)=>a.Gene.localeCompare(b.Gene)||a['Gene ID'].localeCompare(b['Gene ID']));
const ids=['search','group','result','clinical-only'];const controls=Object.fromEntries(ids.map(id=>[id,document.getElementById(id)]));let page=0;const size=50;
for(const result of [...new Set(data.map(r=>r.Result))].sort()){const option=document.createElement('option');option.textContent=result;controls.result.append(option)}
function element(tag,text){const node=document.createElement(tag);node.textContent=String(text);return node}
function update(){const q=controls.search.value.trim().toLowerCase();const filtered=data.filter(r=>(!q||[r.Gene,r['Gene ID'],r['Also known as']||''].join(' ').toLowerCase().includes(q))&&(!controls.group.value||r.Group===controls.group.value)&&(!controls.result.value||r.Result===controls.result.value)&&(!controls['clinical-only'].checked||r['Clinical evidence']));const pages=Math.max(1,Math.ceil(filtered.length/size));page=Math.min(page,pages-1);const body=document.getElementById('results');body.replaceChildren();for(const r of filtered.slice(page*size,(page+1)*size)){const tr=document.createElement('tr');const gene=element('td',r.Gene||'Unnamed source');const details=document.createElement('details');details.append(element('summary','Evidence'));const dl=document.createElement('dl');for(const [key,value] of Object.entries(r)){if(['Gene','Group','Result'].includes(key)||value==='')continue;dl.append(element('dt',key));const dd=document.createElement('dd');if(key==='Clinical source'&&String(value).startsWith('https://')){const link=element('a','Open primary source');link.href=value;dd.append(link)}else{dd.textContent=String(value)}dl.append(dd)}details.append(dl);gene.append(details);tr.append(gene,element('td',r.Group));const result=element('td',r.Result);if(r.Result==='Included')result.className='included';if(r.Result==='Safety hold')result.className='hold';tr.append(result,element('td',r.Reason),element('td',r['Clinical evidence']||'—'));body.append(tr)}if(!filtered.length){const tr=document.createElement('tr');const td=element('td','No matching genes. Try a shorter name or clear the filters.');td.colSpan=5;tr.append(td);body.append(tr)}document.getElementById('status').textContent=filtered.length.toLocaleString()+' matching genes / source identities';document.getElementById('page').textContent='Page '+(page+1)+' of '+pages;document.getElementById('prev').disabled=page===0;document.getElementById('next').disabled=page===pages-1}
for(const control of Object.values(controls)){control.addEventListener(control.type==='search'?'input':'change',()=>{page=0;update()})}document.getElementById('clear').addEventListener('click',()=>{for(const control of Object.values(controls)){if(control.type==='checkbox')control.checked=false;else control.value=''}page=0;update();controls.search.focus()});document.getElementById('prev').addEventListener('click',()=>{page--;update()});document.getElementById('next').addEventListener('click',()=>{page++;update()});update();
</script></body></html>"""
