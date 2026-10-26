#!/usr/bin/env node
/** Reproduce the reviewed 2026-10-24 CTA logic deck; policy edits require review. */
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {parseArgs} from 'node:util';
import {createReportDeck} from './report_deck.mjs';

const {values:args}=parseArgs({options:{
  'out-dir':{type:'string'}, 'artifact-tool':{type:'string'}, help:{type:'boolean'},
}});
if(args.help){
  console.log('Usage: node scripts/slides/cta_logic.mjs [--out-dir DIR] [--artifact-tool MODULE]');
  console.log('Requires @oai/artifact-tool 2.8.99 (or a compatible runtime). See scripts/slides/README.md.');
  process.exit(0);
}
const ROOT=fileURLToPath(new URL('../../',import.meta.url));
const OUT=path.resolve(args['out-dir']??path.join(ROOT,'tmp/pdfs/cta-logic-slides'));
const modulePath=args['artifact-tool']??process.env.ONCOREF_ARTIFACT_TOOL_MODULE;
const {Presentation,PresentationFile}=await import(modulePath?pathToFileURL(path.resolve(modulePath)).href:'@oai/artifact-tool');
const release=JSON.parse(await fs.readFile(path.join(ROOT,'docs/audits/cta-release-1.19.1/summary.json'),'utf8'));
const review=JSON.parse(await fs.readFile(path.join(ROOT,'docs/audits/cta-neither-review-20261024/summary.json'),'utf8'));
assert.deepEqual(release.counts,{Unsafe:1580,Neither:561,CTA:399});
assert.equal(Object.values(release.counts).reduce((a,b)=>a+b),2540);
assert.equal(release.details.filter(x=>x.Category==='Neither').reduce((a,b)=>a+b.Genes,0),561);
assert.equal(release.details.filter(x=>x.Category==='Unsafe').reduce((a,b)=>a+b.Genes,0),1580);
assert.equal(Object.values(review.genome).reduce((a,b)=>a+b),20097);
assert.equal(Object.values(review.proteins).reduce((a,b)=>a+b),129985);
assert.equal(37+88,125); assert.equal(436-37,399); assert.equal(649-88,561);

const metadata={policy:release.policy,commit:'2e40a0f2e6e2b18c919d1a656e7594f0afbe784c'};
const {p,pages,C,text,slide,note,table,box,connect,bullet,bars}=createReportDeck(Presentation,metadata);
// 1
{
const s=slide('CTA selection and peptide availability','OncoRef / 24 October 2026','Policy oncoref.target-sets.2026-10-24.1 | Code snapshot 2e40a0f | Package 1.19.2');
text(s,'One broad modeling universe.\nA separate, evidence-based target screen.',64,185,1140,145,45,C.ink,true);
text(s,'399',64,365,350,115,100,C.teal,true);
text(s,'CTA genes retained from 2,540 mapped nominations',400,389,795,95,34,C.ink);
note(s,'Current rules, two funnel figures, review holds and the limits of the evidence.',565, C.muted,27);
}
// 2
{
const s=slide('Include broadly for modeling; select targets separately','Purpose','[1,2] Target policy | [6,7] Processing, abundance and presentation');
table(s,[['Current label','Peptide-availability model','CTA target set'],['CTA','Include with sequence + assay context','Eligible under this research screen'],['Neither','Include; preserve unresolved evidence','Held for further evidence'],['Unsafe','Include; useful normal-tissue examples','Excluded; forbidden background']],{h:300,widths:[210,480,462],size:25});
text(s,'The six reviewed proteins belong in the modeling universe.\nTheir target labels do not establish peptide availability.',64,500,1152,90,29,C.ink,true);
note(s,'An Unsafe protein can yield presented peptides. An unobserved peptide is not a proven zero.',602,C.muted,22);
}
// 3
{
const s=slide('Define the quantity before comparing model classes','Modeling recommendation','[6] O’Donnell et al., 2020 | [7] Bassani-Sternberg et al., 2015');
const a=box(s,'Precursor\n+ cell context',64,180,245,105,C.ink,29);
const b=box(s,'Processing\nand turnover',359,180,245,105,C.teal,29);
const c=box(s,'Available\npeptide pool',654,180,245,105,C.teal,29);
const d=box(s,'Surface\npeptide-MHC',949,180,267,105,C.ink,29);
connect(s,a,b,'right','left');connect(s,b,c,'right','left');connect(s,c,d,'right','left');
bullet(s,'Per-precursor yield','Factoring out starting abundance gives a conditional yield, not absolute peptide concentration.',335);
bullet(s,'Absolute concentration','Needs abundance, uptake, turnover, compartment, time and assay calibration; surface presentation is another endpoint.',445);
note(s,'Compare RL and differentiable models on the same endpoint, data and leakage-controlled splits.',582,C.teal,26);
}
// 4
{
const s=slide('The policy combines distinct, pinned evidence sources','Evidence','[2,3] Release review and normal-cell evidence guide | HPA v23');
table(s,[['Evidence layer','What it contributes','Keep distinct'],['HPA bulk RNA + curated IHC','Tissue ceilings, selectivity, protein gate','Consensus RNA is not paired with IHC donors'],['Primary antibody records','Staining, validation and conflicts','Primary staining vs curated consensus'],['HPA single-cell + immune RNA','Cell identity and rare-cell review','nTPM; reads; source-native summaries'],['DICE + adult corneal atlas','Immune activation and eye coverage','TPM vs CP10K; no invented donor pairing'],['Reviewed publications + sequences','Admission support, normal protein, identity','Cancer evidence does not waive normal gates']],{h:380,widths:[325,430,397],size:23});
note(s,'44 added source artifacts use a shared cache, checksums and uniform observations API.',584,C.muted,24);
}
// 5
{
const s=slide('Five anatomy categories drive the normal-tissue screen','Tissue policy','[1,2,3] Tissue registry and current policy; examples shown, not an exhaustive anatomy list');
table(s,[['Category','Current examples','Role'],['Reproductive','Testis, ovary, tracts, mammary, placenta, extraembryonic and explicitly fetal/embryonic tissue','Provides reproductive evidence'],['Disposable','Hair follicles, cutaneous melanocytes, tonsil, appendix, gallbladder','Permitted; no normal-expression penalty'],['Manageable','Other skin, spleen, thymus, thyroid, adipose','RNA + selectivity + protein gate'],['Protected','Salivary glands, lymph nodes, cartilage, stomach/colon/rectum','RNA + selectivity + protein gate'],['Vital','Heart, lung, CNS, eye, liver, kidney, pancreas, muscle, blood/marrow, small intestine','RNA + selectivity + protein gate']],{h:410,widths:[190,620,342],size:23,rowHeights:[55,79,69,62,69,76]});
note(s,'Cell identity matters: testicular immune cells are not permitted merely because of their location.',605,C.muted,22);
}
// 6
{
const s=slide('The broad reproductive union has explicit boundaries','Anatomy and species','[1,3] Anatomy registry and normal-cell mapping safeguards');
bullet(s,'Developmental context','Explicit fetal liver is Reproductive. Adult maternal liver remains Vital; pregnancy does not relabel maternal organs.',180);
bullet(s,'Sex and fertility','The broad union is not a pregnancy-safe or fertility-preserving policy. Those require narrower permitted compartments.',300);
bullet(s,'Cross-species use','Anatomy can map across human, mouse and dog. These current measurements and numerical rules are human only.',420);
note(s,'Unresolved anatomy is not silently permitted. Missing measurements never become zero.',574,C.teal,27);
}
// 7
{
const s=slide('Decision precedence: exclusion evidence wins','Classification logic','[1,3] target_sets, normal_cell_screen and normal_rna_review; logical precedence, not runtime order');
const rows=[
 ['Decisive normal protein, safety hold, or any non-permitted RNA >30','Unsafe',C.rust],
 ['Reviewed donor-variable vital RNA with negative local IHC','Neither / Variable',C.amber],
 ['Other failures of absolute tissue or donor RNA ceilings','Unsafe',C.rust],
 ['Unresolved admission, coverage, selectivity, cell or antibody evidence','Neither',C.amber],
 ['Nominated candidate satisfies all remaining requirements','CTA',C.teal]
];
rows.forEach(([label,out,color],i)=>{const y=170+i*80;const a=box(s,`${i+1}. ${label}`,64,y,848,65,C.pale,24);a.text.style={fontSize:24,color:C.ink,bold:false};const b=box(s,out,969,y,247,65,color,25);connect(s,a,b,'right','left');});
note(s,'Read top to bottom. Variable cannot override confirmed normal protein, safety or the >30 rule.',594,C.muted,23);
}
// 8
{
const s=slide('Admission requires affirmative reproductive evidence','Candidate admission','[1,3,4] Current target policy; CT55 and XAGE5 reviews');
text(s,'Start with a mapped curated CTA nomination, usable coverage\nand reproductive bulk RNA R > 0.',64,170,1152,90,30,C.ink,true);
table(s,[['At least one additional admission route','Example or interpretation'],['Reproductive bulk RNA R >=5 nTPM','The current RNA admission floor'],['Reliable positive reproductive protein','The 5 nTPM floor no longer applies'],['Reviewed reproductive + tumor/recognition evidence','CT55 qualifies through cited support']],{y:285,h:260,widths:[600,552],size:25});
note(s,'Admission does not bypass normal-tissue gates. A nomination citation alone is not a rescue.',584,C.teal,25);
}
// 9
{
const s=slide('Absolute RNA ceilings depend on local IHC evidence','Threshold grid','[1,2] Current bulk screen; values are project research parameters, not clinical safety limits');
table(s,[['nTPM ceiling','Baseline','Approved\nnegative','Supported\nnegative','Enhanced\nnegative'],['Vital','2','4','6','12'],['Protected','6','8','12','24'],['Manageable','10','20','30','30 effective*']],{h:285,widths:[280,205,220,220,227],size:28,rowHeights:[87,66,66,66]});
text(s,'Baseline = missing, uncertain or incomplete negative IHC,\nor absence of a reliable positive reproductive control.',64,495,1152,80,27,C.ink);
note(s,'* The underlying Enhanced ceiling is 60; the independent >30 rule overrides it.',599,C.muted,23);
}
// 10
{
const s=slide('Selectivity uses the reproductive maximum','Selectivity grid','[1,2] cta_screen and bulk-tissue rule export');
table(s,[['Minimum share','Baseline','Approved\nnegative','Supported\nnegative','Enhanced\nnegative'],['Vital','80%','75%','70%','65%'],['Protected','75%','70%','65%','60%'],['Manageable','70%','65%','60%','55%']],{h:285,widths:[280,205,220,220,227],size:28,rowHeights:[87,66,66,66]});
text(s,'R = highest reproductive bulk nTPM\nN = the individual normal tissue’s bulk nTPM',64,490,1152,86,29,C.ink);
note(s,'This is a tissue-by-tissue screen. It is not a whole-profile breadth or entropy score.',597,C.teal,25);
}
// 11
{
const s=slide('All three bulk RNA inequalities must hold','Worked example','[1,2] Implemented formula; fixed 0.1 nTPM allowance is not an uncertainty interval');
text(s,'N <= C       N <= R × (1 - s) / s + 0.1       N < R',64,178,1152,75,37,C.teal,true);
text(s,'C is the category/IHC ceiling; s is the minimum selectivity fraction.',64,260,1152,55,26,C.muted);
table(s,[['Example: R = 5; Vital baseline','Normal RNA','Tissue result'],['Absolute ceiling C = 2; s = 0.80','1.0 nTPM','Passes these RNA checks'],['Selectivity limit = 1.35 nTPM','1.5 nTPM','Fails selectivity only'],['All three inequalities must hold','2.5 nTPM','Fails absolute ceiling']],{y:340,h:245,widths:[575,230,347],size:25});
note(s,'Selectivity-only failure normally means Neither; an absolute RNA failure normally means Unsafe.',606,C.muted,22);
}
// 12
{
const s=slide('Negative IHC relaxes RNA under explicit conditions','Protein evidence','[1,3] Curated IHC vs primary antibody rules | Six-gene source audit, 24 October 2026');
bullet(s,'For relaxation','Every routine measured cell compartment in that tissue is negative, with qualifying confidence and a reliable reproductive positive control.',176);
bullet(s,'For exclusion','Accepted positive curated IHC, including Low, in Vital/Protected/Manageable tissue excludes. Cited normal protein can also exclude.',309,C.rust);
bullet(s,'Primary-record caveat','Uncertain-profile normal staining causes an Antibody hold. Reliable-profile primary positives can remain discordant with negative consensus.',442,C.amber);
note(s,'The control is not a paired-donor RNA-to-protein calibration. YBX2 and ACTRT1 need conflict review.',603,C.muted,22);
}
// 13
{
const s=slide('Cell-resolved RNA adds source-specific review triggers','Normal-cell review','[3] normal_cell_screen; native units remain separate');
table(s,[['Source','Review trigger (before stronger exclusions)'],['HPA tissue/cell clusters','>=5 nTPM, >=20 reads, normal >= same-source reproductive-cell maximum'],['Original DICE','>=5 TPM in at least 3 reported values within a cell/state distribution'],['HPA / Monaco immune summaries','>=5 in that source’s native TPM or nTPM'],['Adult corneal atlas','Mean CP10K >=1; >=5 positive cells; >=20 counts; in >=2 sample groups'],['Primary antibody records','Positive non-permitted staining with Uncertain profile reliability']],{h:380,widths:[360,792],size:25,rowHeights:[55,80,65,55,70,55]});
note(s,'These normally produce Neither; activated immune cells have no automatic exemption.\nSignificant unresolved anatomy or reproductive-cell coverage also requires review.',581,C.muted,22);
}
// 14
{
const s=slide('The >30 normal-RNA rule is an independent hard stop','High RNA','[1,2,3] normal_rna_review; applied to prior CTA/Neither and reference loci in the audit');
text(s,'Any measured non-permitted value >30 TPM or >30 nTPM\ncauses Unsafe / High RNA.',64,175,1152,110,36,C.rust,true);
bullet(s,'Scope','Outside explicitly mapped Reproductive or Disposable contexts, including unsited observations. Exactly 30 does not trigger it.',324);
bullet(s,'No rescue','One high observation is sufficient. Negative IHC, high reproductive expression and a Variable hold do not override it.',442);
note(s,'TPM and nTPM are checked separately. Corneal CP10K keeps its own rule; it is not TPM.',593,C.muted,24);
}
// 15
{
const s=slide('Variable records a path to reconsideration','Neither / Variable','[1,2] Reviewed vital-tissue donor RNA; currently CTAG2 only');
table(s,[['Requirement','Meaning'],['Same gene, vital tissue, cohort and version','Reported zero in one donor and a threshold failure in another'],['Negative local IHC','Required protein context; donor RNA and IHC are not asserted to be paired'],['Threshold follows current category / IHC policy','There is no universal >10 nTPM definition of Variable'],['No stronger exclusion','Normal protein, a safety hold or >30 RNA can still make it Unsafe']],{h:315,widths:[460,692],size:26});
text(s,'CTAG2 remains Variable: excluded from CTA targets and from\nthe default forbidden background, pending better evidence.',64,530,1152,90,30,C.amber,true);
}
// 16
{
const s=slide('Current funnel: 399 of 2,540 nominations retained','Final classification funnel','[2] Release 1.19.1 frozen gene table; same policy in package 1.19.2');
const a=box(s,'2,540 mapped nominations',64,171,728,88,C.ink,33);
const b=box(s,'960 without an Unsafe assignment',144,322,568,88,C.ink,28);
const c=box(s,'399 CTA',270,474,316,88,C.teal,39);
connect(s,a,b);connect(s,b,c);
const u=box(s,'1,580 Unsafe\n62.2% of nominations',880,246,336,100,C.rust,28);
const n=box(s,'561 Neither\n22.1% of nominations',880,441,336,100,C.amber,28);
connect(s,a,u,'right','left');connect(s,b,n,'right','left');
note(s,'CTA = 15.7%. BAGE and GAGE2B remain unmapped intake, outside this denominator.',601,C.muted,24);
}
// 17
{
const s=slide('Unsafe: most principal reasons are protein or bulk RNA','Exclusion funnel detail','[2] Mutually exclusive assigned reasons; a gene can have several underlying findings');
bars(s,['Normal IHC','RNA ceiling','High RNA >30','Literature','Safety hold'],[845,606,125,3,1],C.rust,{h:375});
note(s,'1,580 total. Unsafe is a screening exclusion, not proof that a specific therapy causes injury.',585,C.muted,25);
}
// 18
{
const s=slide('Neither preserves uncertainty without assuming absence','Review funnel detail','[2] Mutually exclusive assigned subgroups; detailed questions remain in follow-up tables');
bars(s,['Antibody','Low evidence','Cellular','Screening','Unassessed','Selectivity','Variable'],[199,191,70,59,31,10,1],C.amber,{h:395});
note(s,'561 total. Neither is neither target-eligible nor forbidden by default.',597,C.muted,26);
}
// 19
{
const s=slide('Last update: 37 retained nominations moved to Unsafe','Policy-change funnel','[1,2] Before/after mapped nominations only; do not mix with genome-wide counts');
text(s,'Previous',64,165,270,38,25,C.muted,true);text(s,'Current',890,165,310,38,25,C.muted,true);
const a=box(s,'436 CTA',64,223,296,88,C.teal,34);
const b=box(s,'649 Neither',64,389,296,88,C.amber,34);
const c=box(s,'399 CTA',890,223,326,88,C.teal,34);
const d=box(s,'561 Neither',890,389,326,88,C.amber,34);
const e=box(s,'125 new Unsafe',490,312,300,88,C.rust,31);
connect(s,a,c,'right','left');connect(s,b,d,'right','left');connect(s,a,e,'right','left');connect(s,b,e,'right','left');
text(s,'37 to Unsafe',64,320,296,45,27,C.rust,true);text(s,'88 to Unsafe',64,487,296,45,27,C.rust,true);
note(s,'Unsafe increased from 1,455 to 1,580. Across the full genome, 290 genes moved to Unsafe.',575,C.muted,25);
}
// 20
{
const s=slide('Gene, protein and k-mer sets have different scopes','Sequence exclusion','[1] Current Ensembl 115 protein inventory and exclusion recipe');
table(s,[['Scope','CTA','Neither','Unsafe','Total'],['Mapped nominations','399','561','1,580','2,540'],['Canonical human genes','399','1,999','17,699','20,097'],['Unique complete proteins','723','2,878','126,384','129,985']],{h:280,widths:[390,160,180,200,222],size:27});
text(s,'Shared-protein priority: Unsafe > Neither > CTA.\nDefault forbidden = Unsafe complete + Unsafe incomplete translations.',64,495,1152,95,28,C.ink,true);
note(s,'Full protein sets are disjoint; k-mers can still overlap. Additional reference loci extend gene scope.',605,C.muted,22);
}
// 21
{
const s=slide('Benchmarks show what the screen preserves and loses','Clinical / requested benchmarks','[1,2] Current gene decisions; examples do not establish clinical efficacy or safety');
table(s,[['Current outcome','Examples / reason'],['CTA retained','PRAME, MAGEA4, CTAG1A/B, XAGE1A/B, CT83, CT55, SAGE1, CTCFL, GAGE1, SSX4/5'],['Neither','CTAG2: Variable; FOXR2 and XAGE5: Antibody; HORMAD1, TDRD1, BRDT, DAZL: Cellular'],['New benchmark exclusions','MAGEA3: colonic cluster 53.1 nTPM; ROPN1: salivary cluster 32.9 nTPM'],['Other high-RNA exclusions','ACRBP and TEX101 now Unsafe / High RNA']],{h:350,widths:[340,812],size:25,rowHeights:[55,85,90,70,50]});
note(s,'19 of the previous 21 retained benchmark genes remain. Published targeting is not an automatic pass.',583,C.muted,24);
}
// 22
{
const s=slide('Six retained proteins need evidence-specific caution','Retained does not mean settled','[2] Retained literature review | Six-gene all-source audit, 24 October 2026');
table(s,[['Gene','Reproductive\nbulk nTPM','Max audited normal\nnTPM / TPM*','Main unresolved point'],['ABCC11','45.3','7.7 / 2.10','Reported neural protein'],['MAGEC2','19.0','1.9 / 0.28','Published Purkinje staining'],['YBX2','535.8','27.0 / 2.07','17/146 primary normal records positive'],['CCDC158','33.3','18.8 / 0.80','Normal isoforms; brain cell RNA'],['ACTRT1','4.0','0.4 / 0.36','5/48 primary normal records positive'],['LGALS13','10.9','6.5 / 0.43','Extraplacental protein claims']],{h:378,widths:[180,220,280,472],size:23,rowHeights:[78,50,50,50,50,50,50]});
note(s,'All six remain CTA. *Separate source maxima, not equivalent units or donor-positive fractions.',595,C.muted,23);
}
// 23
{
const s=slide('What the present screen still cannot establish','Limits and next decisions','[1,2,3] Frozen results; proposed next work below is not an implemented policy change');
bullet(s,'Whole-profile specificity','A high reproductive maximum can conceal diffuse lower expression. A continuous breadth metric is not yet implemented.',174);
bullet(s,'Protein concordance','Negative consensus can coexist with positive primary records. ABCC11, MAGEC2, YBX2 and CCDC158 remain priority reviews.',289);
bullet(s,'Coverage and identity','Rare cells, cell state, donor variation and isoform-specific protein production are incompletely resolved.',404);
note(s,'Use the broad dataset for modeling; keep review holds, assay uncertainty and safety labels separate.',565,C.teal,28);
}
// 24
{
const s=slide('Use the current target-set API','Reproducibility','[1,3,4] Public API and pinned source contract');
text(s,'from oncoref import target_sets, normal_cell_data as cells\n\nresults = target_sets.results()\ncta = target_sets.gene_ids("CTA")\nunsafe = target_sets.gene_ids("Unsafe")\nneither = target_sets.gene_ids("Neither")\n\nevidence = cells.observations("hpa_rna_consensus",\n                              genes=["CTAG2"])',64,174,1152,335,28,C.ink);
note(s,'Current: policy 2026-10-24.1. Earlier cta / cta_screen accessors preserve dated definitions.',545,C.teal,25);
text(s,'Raw cache: ~/.cache/oncoref/sources/<source>/<version>/<filename>',64,605,1152,36,22,C.muted);
}
// 25
{
const s=slide('Sources and navigation','References / 1 of 2','All report links refer to published OncoRef documentation; policy snapshot 24 October 2026');
const refs=[
 ['[1] Current target sets and protein exclusion recipe','https://pirl-unc.github.io/oncoref/target-sets/'],
 ['[2] Full release review, thresholds, reasons and caveats','https://pirl-unc.github.io/oncoref/audits/cta-release-1.19.1/'],
 ['[3] Normal-cell sources, units, API and rule details','https://pirl-unc.github.io/oncoref/normal-cell-evidence/'],
 ['[4] XAGE5: historical rescue vs current antibody hold','https://pirl-unc.github.io/oncoref/cta-xage5-review/'],
 ['[5] Full CTA atlas: individual genes, cancer occurrence and literature','https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/']
];
refs.forEach(([label,url],i)=>{const y=173+87*i;text(s,label,64,y,1152,35,25,C.ink,true);text(s,url,64,y+38,1152,32,21,C.teal);pages.at(-1).links.push({url,x:64,y:y+38,w:1152,h:32});});
}
// 26
{
const s=slide('Primary literature for the modeling distinction','References / 2 of 2','Literature supports separate processing, abundance and presentation factors; no clinical threshold validation claimed');
text(s,'[6] O’Donnell, Rubinsteyn & Laserson (2020)',64,177,1152,45,30,C.ink,true);
text(s,'MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented\nPeptides by Incorporating Antigen Processing. Cell Systems.',64,230,1152,80,27,C.ink);
text(s,'https://doi.org/10.1016/j.cels.2020.06.010',64,320,1152,35,24,C.teal);
pages.at(-1).links.push({url:'https://doi.org/10.1016/j.cels.2020.06.010',x:64,y:320,w:1152,h:35});
text(s,'[7] Bassani-Sternberg et al. (2015)',64,395,1152,45,30,C.ink,true);
text(s,'Mass Spectrometry of Human Leukocyte Antigen Class I Peptidomes\nReveals Strong Effects of Protein Abundance and Turnover on\nAntigen Presentation. Molecular & Cellular Proteomics.',64,448,1152,115,27,C.ink);
text(s,'https://doi.org/10.1074/mcp.M114.042812',64,583,1152,35,24,C.teal);
pages.at(-1).links.push({url:'https://doi.org/10.1074/mcp.M114.042812',x:64,y:583,w:1152,h:35});
}

await fs.mkdir(path.join(OUT,'rendered'),{recursive:true});
await fs.writeFile(path.join(OUT,'manifest.json'),JSON.stringify({policy:release.policy,commit:'2e40a0f2e6e2b18c919d1a656e7594f0afbe784c',slideSize:{width:1280,height:720},pages},null,2));
for(let i=0;i<p.slides.items.length;i++){
  const s=p.slides.items[i],stem=`slide-${String(i+1).padStart(2,'0')}`;
  const png=await p.export({slide:s,format:'png',scale:2});
  await fs.writeFile(path.join(OUT,'rendered',`${stem}.png`),new Uint8Array(await png.arrayBuffer()));
  const layout=await s.export({format:'layout'});
  await fs.writeFile(path.join(OUT,'rendered',`${stem}.layout.json`),typeof layout==='string'?layout:JSON.stringify(layout));
  console.log(`${i+1}/${p.slides.items.length} ${pages[i].title}`);
}
// Editable native source accompanies the previews and manifest.
const pptx=await PresentationFile.exportPptx(p);await pptx.save(path.join(OUT,'source-deck.pptx'));
console.log('All slides rendered and native source saved.');
