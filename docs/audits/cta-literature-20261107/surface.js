"use strict";
function frequency(n, d) { return d ? `${n}/${d} (${(100*n/d).toFixed(1)}%)` : "—"; }
function surfaceRNA(r) {
  const cohort=$("surface-cohort").value||"TCGA", cancer=$("surface-cancer").value;
  return cancer ? data.surface_cancer_rna.filter(x=>x.entity_id===r.entity_id&&x.cohort===cohort&&x.cancer===cancer) : r.top_rna[cohort];
}
function renderSurface() {
  const q=$("surface-query").value.toLowerCase().trim(),sort=$("surface-sort").value;
  const rows=data.surface_domains.filter(r=>!q||JSON.stringify(r).toLowerCase().includes(q));
  const prevalence=r=>Math.max(-1,...surfaceRNA(r).map(x=>x.samples?x.expressed_samples_ptpm_ge_1/x.samples:-1));
  rows.sort((a,b)=>(sort==="size"?(b.extracellular_longest_aa??-1)-(a.extracellular_longest_aa??-1):sort==="prevalence"?prevalence(b)-prevalence(a):0)||a.protein.localeCompare(b.protein));
  $("surface-count").textContent=`${rows.length} of 37 membrane candidates · RNA and IHC are separate cohorts`;
  $("surface-table").innerHTML=slice("surface",rows,renderSurface).map(r=>{
    const staining=data.surface_cancer_ihc.filter(x=>x.entity_id===r.entity_id&&x.measurement_status==="measured").sort((a,b)=>b.prevalence_medium_high-a.prevalence_medium_high||a.cancer.localeCompare(b.cancer));
    const match=r.sequence_match==="exact_sequence"?"Exact selected-sequence match":r.sequence_match==="mature_chain_identical_anchor_tail_differs"?"Mature chain matches; anchor tail differs":"Different from selected sequence";
    const features=r.topology_features.map(f=>`<p>${esc(f.type)}: ${esc(f.description||"")} ${esc(f.location.start.value??"?")}–${esc(f.location.end.value??"?")}<br><span class="meta">${(f.evidences||[]).map(e=>esc(e.evidenceCode)).join("; ")||"No evidence code"} ${nativeEvidence(f)}</span></p>`).join("");
    return `<tr><td>${gene(r)}<p><a href="${esc(r.uniprot_url)}">${esc(r.accession)}</a> · ${r.uniprot_length} aa</p><p class="meta">Selected: ${r.selected_lengths.join(", ")} aa<br>${esc(match)}</p></td>
    <td><b>${r.extracellular_total_aa===null?"Unknown":r.extracellular_total_aa+" aa total; "+r.extracellular_longest_aa+" aa longest"}</b><p>${r.extracellular_intervals.map(x=>x.join("–")).join(", ")||"Outside-facing boundaries unannotated"}</p><p class="meta">${esc(r.basis)}<br>${r.extracellular_evidence_codes.map(esc).join("; ")||"No extracellular evidence code"}</p>${r.notes.map(n=>`<p>${esc(n)}</p>`).join("")}${r.sequence_match==="different_sequence"?'<p><b>Selected-sequence size: unknown.</b></p>':""}<details><summary>Native feature coordinates and evidence</summary>${features}</details></td>
    <td>${surfaceRNA(r).map(x=>`<p><b>${esc(x.cancer_type)}</b><br>≥1: ${frequency(x.expressed_samples_ptpm_ge_1,x.samples)}<br>≥5: ${frequency(x.expressed_samples_ptpm_ge_5,x.samples)}</p>`).join("")||"— No measurement for this cohort/cancer."}<p class="meta">Top-three mode ranks by the fraction at ≥1 pTPM; ties use ≥5 pTPM, then cancer name.</p></td>
    <td>${staining.length?`<details><summary>IHC: ${staining.length} cancer categories</summary><p>Medium/High staining among scored patients:</p>${staining.map(x=>`<p>${esc(x.cancer)}: ${frequency(x.medium_high,x.total)}<br><span class="meta">High ${x.high}; Medium ${x.medium}; Low ${x.low}; Not detected ${x.not_detected}</span></p>`).join("")}</details>`:"— No HPA cancer IHC measurements."}${data.concerns.filter(c=>c.entity_id===r.entity_id).map(c=>`<p><b>${esc(c.category)}</b><br>${esc(c.observation)} <a href="?protein=${encodeURIComponent(r.protein)}#concerns">Full finding and papers</a></p>`).join("")}</td></tr>`;
  }).join("")||'<tr><td colspan="4">No matching membrane candidates.</td></tr>';
}
function initSurface() {
  page.surface=0;
  const params=new URLSearchParams(location.search);
  if(params.get("protein"))$("surface-query").value=params.get("protein");
  const populate=()=>{
    $("surface-cancer").replaceChildren();$("surface-cancer").add(new Option("Top three for each protein",""));
    for(const name of [...new Set(data.surface_cancer_rna.filter(r=>r.cohort===($("surface-cohort").value||"TCGA")).map(r=>r.cancer))].sort())$("surface-cancer").add(new Option(name,name));
    $("surface-cancer").value="";
  };
  populate();bind(["surface-query","surface-cancer","surface-sort"],"surface",renderSurface);
  $("surface-cohort").addEventListener("input",()=>{page.surface=0;populate();renderSurface();});
  const s=data.surface_sources;
  $("surface-methods").innerHTML=`<p>Topology: <a href="${esc(s.uniprot_receipt.url)}">UniProt ${esc(s.uniprot_receipt.headers["X-UniProt-Release"])}</a>, retrieved ${esc(s.uniprot_receipt.retrieved)}. The separate 602-protein localization snapshot remains version 2026_03. Coordinates are 1-based and inclusive.</p><p>ECO:0000255: sequence prediction; ECO:0000250: similarity; ECO:0000305: curator inference; ECO:0000269: experimental evidence. Each feature retains its own evidence code and citation; an experimentally detected protein does not make every feature experimental.</p><p>Cancer measurements: <a href="https://www.proteinatlas.org/humanproteome/cancer/data">HPA ${esc(s.hpa_source.hpa_release)}</a>, native Ensembl ${esc(s.hpa_source.ensembl_release)} gene IDs matched exactly to the released candidate IDs. RNA denominators count measured samples. IHC denominators count scored patients, not cells. TCGA, validation and IHC are never pooled.</p><p><a href="${esc(s.hpa_source.artifacts.hpa_cancer_rna.url)}">Pinned genome-wide RNA source</a> · <a href="${esc(s.hpa_source.artifacts.hpa_cancer_ihc.url)}">Pinned IHC source</a> · <a href="https://www.uniprot.org/help/sequence_annotation">UniProt feature definitions</a></p>`;
  renderSurface();
}
