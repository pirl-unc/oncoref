// Editable release deck from presentation-data.json; scientific scatter from the audit.
// Required: RUNTIME_NODE_MODULES, ONCOREF_PRESENTATIONS_SKILL, RUNTIME_PYTHON.
// Run render_cta_source_walkthrough.py first. Output goes to the versioned audit.
import fs from 'node:fs/promises';
import {readFileSync} from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
import {fileURLToPath, pathToFileURL} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const audit=path.join(root,'docs/audits/cta-source-walkthrough-20261028');
const build=path.resolve(root,process.env.CTA_SLIDE_BUILD_DIR||'tmp/cta-source-slides-expanded');
const data=JSON.parse(await fs.readFile(path.join(audit,'presentation-data.json'),'utf8'));
const trials=JSON.parse(await fs.readFile(path.join(root,'oncoref/data/tcr-trial-evidence.json'),'utf8'));
const definition=JSON.parse(await fs.readFile(path.join(root,'docs/audits/cta-definition-20261104/report-data.json'),'utf8'));
const inventory=JSON.parse(await fs.readFile(path.join(root,'docs/audits/cta-inventory-20261105/summary.json'),'utf8'));
const sequences=JSON.parse(await fs.readFile(path.join(root,'docs/audits/sequence-bundle-20261104/manifest.json'),'utf8'));
const skill=process.env.ONCOREF_PRESENTATIONS_SKILL;
const runtime=process.env.RUNTIME_NODE_MODULES;
if(!skill||!runtime) throw Error('Set presentation skill and runtime paths');
const req=createRequire(path.join(runtime,'__runtime__.cjs'));
const {Presentation,PresentationFile}=await import(pathToFileURL(req.resolve('@oai/artifact-tool')).href);
const {finalizePresentation,resolvePresentationFont,applyPresentationChartFont}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const font=resolvePresentationFont();
const p=Presentation.create({slideSize:{width:1280,height:720}});
const C={ink:'#173B55',muted:'#536675',teal:'#087E8B',amber:'#A76708',rust:'#AC493B',pale:'#F3F6F8',line:'#D9E2E8',white:'#FFFFFF'};
const URL='https://pirl-unc.github.io/oncoref/audits/cta-source-walkthrough-20261028/';
const pages=[],tableOwners=[],chartOwners=[];let current;
const linksByName=new Map([
 ...data.papers.map(r=>[r.citation,r.source_url]),
 ...(data.paper_titles||[]).flatMap(r=>[[r.citation,'https://doi.org/'+r.doi],[r.title,'https://doi.org/'+r.doi]]),
 ['HPA v23','https://v23.proteinatlas.org/about/download'],
 ['GTEx v11','https://gtexportal.org/home/downloads/adult-gtex'],
 ['HPA v25.1','https://www.proteinatlas.org/about/releases'],
 ['Tyser et al. (2021)','https://doi.org/10.1038/s41586-021-04158-y'],
 ['Petropoulos et al. (2016)','https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results'],
 ['DICE','https://dice-database.org/'],
 ['GSE155683','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683'],
 ['Shaw et al. (2024)','https://doi.org/10.1038/s41467-024-47649-y'],
 ['TRIM71 liver study (2024)','https://pmc.ncbi.nlm.nih.gov/articles/PMC11388079/'],
 ['Cardoso-Moreira et al. (2019)','https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/'],
 ['Clinical trial and SAE review','https://pirl-unc.github.io/oncoref/audits/cta-tcr-trials-20261029/'],
]);
await fs.mkdir(path.join(build,'rendered'),{recursive:true});
await fs.mkdir(path.join(audit,'slides'),{recursive:true});
function linkedText(owner,value,box){
 const content=String(value),matches=[];
 for(const [label,url] of linksByName){const at=content.indexOf(label);if(at>=0)matches.push({at,label,url});}
 matches.sort((a,b)=>a.at-b.at||b.label.length-a.label.length);
 let end=0;const runs=[],urls=new Set();
 for(const m of matches){if(m.at<end)continue;if(m.at>end)runs.push({run:content.slice(end,m.at)});runs.push({run:m.label,style:{color:C.teal,underline:'sng'},link:{uri:m.url,isExternal:true}});end=m.at+m.label.length;urls.add(m.url);}
 if(!runs.length){owner.text.set(content);return;}
 if(end<content.length)runs.push({run:content.slice(end)});
 owner.text.set([runs]);
 current.links.push({url:urls.size===1?[...urls][0]:URL+'references/',...box});
}
function text(s,t,x,y,w,h,size=27,color=C.ink,bold=false){
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 linkedText(sh,t,{x,y,w,h});sh.text.style={typeface:font,fontSize:size,color,bold};sh.text.insets=0;current.text.push(String(t));return sh;
}
function slide(title,section,notes=''){
 const s=p.slides.add();s.background.fill=C.white;current={title,text:[],links:[{url:URL,x:64,y:670,w:1100,h:24}]};pages.push(current);
 text(s,section.toUpperCase(),64,28,1152,25,16,C.teal,true);
 text(s,title,64,72,1152,85,42,C.ink,true);
 text(s,'OncoRef CTA evidence walkthrough   •   5 November 2026',64,670,1100,26,16,C.muted);
 text(s,String(pages.length).padStart(2,'0'),1174,668,42,30,19,C.muted,true);
 s.speakerNotes.text=`${title}\n${notes}\nReport: ${URL}\nPolicy: ${data.summary.policy}\nDecision snapshot SHA256: ${data.report_sha256}`;
 return s;
}
function note(s,t,y=598,size=23){text(s,t,64,y,1152,62,size,C.muted);}
function table(s,values,{x=64,y=170,w=1152,h=390,widths,size=24}={}){
 tableOwners.push(pages.length);
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnWidths:widths,values});
 t.borders.assign({style:'solid',fill:C.line,width:.5});
 t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({margins:{top:2,bottom:2,left:10,right:10}});
 values.forEach((row,r)=>{t.rows[r].height=h/values.length;row.forEach((_,c)=>{const cell=t.getCell(r,c);if(r>0)linkedText(cell,row[c],{x:x+(widths||Array(row.length).fill(w/row.length)).slice(0,c).reduce((a,b)=>a+b,0),y:y+r*h/values.length,w:(widths||Array(row.length).fill(w/row.length))[c],h:h/values.length});cell.fill=r===0?C.ink:(r%2===0?C.pale:C.white);cell.text.style={typeface:font,fontSize:size,color:r===0?C.white:C.ink,bold:r===0||c===0};});current.text.push(row.join(' | '));});return t;
}
function chart(s,categories,series,{x=64,y=170,w=1152,h=390,stacked=false,size=22,labels=false}={}){
 chartOwners.push(pages.length);
 if(series.some(r=>r.values.length!==categories.length||r.values.some(v=>!Number.isFinite(v))))throw Error('Chart labels and finite data must align');
 const ch=s.charts.add('bar',{position:{left:x,top:y,width:w,height:h},categories:[...categories].reverse(),series:series.map(a=>({...a,values:[...a.values].reverse()})),barOptions:{direction:'bar',grouping:stacked?'stacked':'clustered',gapWidth:40},hasLegend:false,legend:{position:'bottom',textStyle:{fontSize:20,fill:C.ink}},xAxis:{textStyle:{fontSize:size,fill:C.ink},majorGridlines:null},yAxis:{textStyle:{fontSize:18,fill:C.muted}},dataLabels:{showValue:labels,position:'outEnd',textStyle:{fontSize:22,fill:C.ink}},chartFill:C.white,plotAreaFill:C.white});
 applyPresentationChartFont(ch,{fontFamily:font});return ch;
}
function figure(s,name,x,y,w,h){s.images.add({blob:new Uint8Array(readFileSync(path.join(audit,'figures',name+'.png'))),contentType:'image/png',fit:'contain',position:{left:x,top:y,width:w,height:h},alt:name.replaceAll('-',' ')});}
const n=v=>Number(v??0).toLocaleString('en-US');


const counts=data.summary.categories;
const release=readFileSync(path.join(root,'oncoref/version.py'),'utf8').match(/__version__\s*=\s*["']([^"']+)/)[1];
{
const s=slide('Human and dog CTA evidence','Current selection','Every human protein has a category. Other includes unmet identity, sequence and review-hold gates. Assessment status records remaining evidence limitations.');
text(s,`${n(counts.CTA)} human CTA proteins`,64,202,1152,100,66,C.teal,true);
text(s,`${n(data.summary.selected_member_loci)} loci   /   ${n(counts.Unsafe)} Unsafe   /   ${n(counts.Other)} Other`,64,335,1152,65,32);
text(s,`Every human protein classified. Dog selection remains unassessed.`,64,446,1152,95,32);
note(s,'Sources, expression profiles, single-cutoff filtering and measured tissue evidence.',580,26);
}
{
const s=slide('Safety, cancer and Reproductive evidence','Independent questions','Only safety supplies Unsafe evidence. Reproductive is the umbrella for adult reproductive, placental, extraembryonic, embryonic and fetal contexts.');
table(s,[['Assessment','Question','Current rule'],['Safety','Concerning normal expression?','Absolute RNA / IHC and clinical rules'],['Cancer expression','Expressed in a measured cancer?','One sample >0 TPM/pTPM or reviewed tumor protein'],['Reproductive evidence','Expression in a qualifying branch?','Bulk/cell RNA ≥5, reliable reproductive protein, or embryonic RNA'],['Coverage','Are required panels measured?','Combined HPA bulk / brain and GTEx for every locus']],{widths:[300,452,400],size:25,h:380,y:170});
note(s,'Only safety establishes Unsafe. Coverage and normal review holds remain explicit in assessment_status.',580,25);
}
{
const s=slide('Diffuse expression is an identity question','HPA profile identity','The old normal-maximum veto is removed. Diffuseness averages N/(N+R) across fixed weighted normal contexts, including Disposable. The active 10% adult/placental identity boundary requires complete HPA measurements and an RNA- or IHC-supported positive reference. Thymus is excluded from concentration; its hazard tier remains Protected. Supported embryonic evidence retains its separate route with unassessed concentration.');
table(s,[['Term','Meaning'],['N','HPA v23 normal-tissue nTPM'],['R','Highest Reproductive nTPM in HPA v23 (adult/placental)'],['N / (N + R)','Average across 29 organs / 38 contexts; thymus excluded'],['Active boundary','Average >10%, with a supported positive R and complete comparison']],{widths:[335,817],size:27,h:320,y:170});
const example=data.concentration_example;
note(s,`${example.label}: ${example.classification} (${example.unsafe_reasons}); diffuseness ${example.diffuseness_pct.toFixed(2)}%. Embryonic concentration remains unassessed.`,565,26);
}
{
const s=slide('Required and optional RNA coverage','Coverage','Combined bulk HPA consensus, HPA brain and GTEx coverage is required for every member locus. Cell, immune and corneal panels contribute measured evidence without imposing coverage requirements.');
table(s,[['Requirement','Interpretation'],['Combined normal bulk','HPA consensus / brain and GTEx; measured zero counts'],['All member loci measured','Missing identities are not replaced by paralogs'],['Missing required coverage','Blocks CTA; assessment status records the gap'],['Cell / immune / cornea panels','Missing coverage allowed; measured evidence counts']],{widths:[350,802],size:26,h:330});
note(s,'The fixed GTEx core is an inspectable profile, not a separate admission gate. Rare-cell coverage remains incomplete.',575,25);
}
{
const s=slide('ZNF730: panel composition matters','Calculation audit','Full tissue tables: '+URL.replace('cta-source-walkthrough-20261028/','cta-unassigned-20261029/znf730/')+' . Original panels differ; scores are not a controlled assay comparison.');
table(s,[['Calculation','HPA v23','GTEx v11'],['Reproductive reference (testis)','4.9 nTPM','4.91414 TPM'],['Full source panel','9.53%: 29 anatomy units','1.65%: 23 anatomy units'],['Same 20 anatomy units','6.36%','1.71%']],{widths:[460,346,346],size:28,h:280,y:180});
text(s,'Thymus is excluded from both sides of the concentration score',64,505,1152,70,28,C.teal,true);
note(s,'Thymus exclusion changes HPA 10.34% → 9.53%. Independent safety evidence keeps ZNF730 Unsafe.',592,25);
}
{
const s=slide('Published candidate lists','Human annotations',data.papers.map(r=>`${r.citation}: ${r.source_url}`).join('\n'));
table(s,[['Paper','Proteins','CTA','Unsafe','Other','Unassigned'],...data.papers.map(r=>[r.citation,n(r.mapped_entities),r.CTA,r.Unsafe,r.Other,r.Unassigned])],{widths:[390,150,140,160,150,162],size:19,h:420,y:158});
note(s,`Union: ${n(data.paper_union_groups)} proteins / ${n(data.paper_union_loci)} loci. ${data.annotated} selected proteins have an annotation; ${data.unannotated} do not.`,600,23);
}
{
const s=slide('Broader published nominations','Publication scope','Sasaroli et al. (2011) https://doi.org/10.4161/cbt.12.3.16260 . Hirbo et al. (2015) https://doi.org/10.1371/journal.pone.0144155 . Exact extracted rows and matching receipts accompany the report.');
table(s,[['Paper','Published scope','Identity accounting'],['Sasaroli et al. (2011)','13 detailed tumor-vascular candidates','13 matched genes'],['Hirbo et al. (2015)','86 placenta-enriched genes, including 30 fast-evolving','82 matched, 4 unresolved'],['Compact source view','Bai, Boël and Bradley covered by larger sets','All memberships retained'],['Garg et al. (2007)','Still uniquely contributes SPAG9','Retained in main view']],{widths:[290,570,292],size:25,h:350});
note(s,'Vascular and placental nominations retain their scope. Paper membership never bypasses the quantitative funnel.',580,26);
}
{
const s=slide('Overlap between published lists','Human annotations','Tables and author-recognized prior antigens within the same paper are unioned. Evidence roles, source IDs and pairwise counts accompany the report.');
figure(s,'published-list-overlaps',36,152,1208,500);
}
{
const s=slide('Reproductive evidence: measured branches','Positive expression contexts','Tyser et al. (2021): https://doi.org/10.1038/s41586-021-04158-y . HPA v23 https://v23.proteinatlas.org/about/download .');
table(s,[['Qualifying route','Entrants','Final CTA'],...data.routes.map(r=>[r.route,n(r.entrants),n(r.CTA)])],{widths:[620,266,266],h:280,size:29});
text(s,'RNA and reliable reproductive IHC establish association',64,495,1152,50,32,C.teal,true);
note(s,'Entrants count positive expression before the other admission gates. Reproductive presence and profile concentration are separate questions.',565,26);
}
{
const s=slide('Normal and cancer measurement sources','Human evidence','HPA https://v23.proteinatlas.org/about/download ; DICE https://dice-database.org/ ; cornea https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683 ; Tyser https://doi.org/10.1038/s41586-021-04158-y . Full URLs and sample lineage accompany the report.');
table(s,[['Evidence','Source','Quantity'],['Normal / reproductive bulk','HPA v23 consensus / brain + GTEx v11','nTPM / TPM'],['Normal tissue cells','HPA v23 clusters, at least 20 reads','nTPM'],['Immune cells','HPA v23 / Monaco / original DICE','Native nTPM / TPM'],['Adult cornea','Collin et al., GSE155683','CP10K'],['Normal protein','HPA v23 + curated external IHC','Staining / reagent confidence'],['Cancer','132 compatible OncoRef cohort views + HPA','TPM / protein-coding TPM'],['Developmental','Tyser gastrulation + Petropoulos preimplantation','Stage / lineage mean TPM']],{widths:[310,570,272],size:22,h:405});
note(s,'DICE and corneal coverage are optional. Measured counterevidence still contributes.',600,24);
}
{
const s=slide('RNA cutoffs','Current policy','The same nominal values apply across native TPM/nTPM assays. CP10K uses nominal /100, not a TPM conversion. Values are screening-policy choices, not measured toxicity limits.');
table(s,[['Ceiling','Missing / unusable','Approved −','Supported −','Enhanced −'],['Vital',4,8,16,32],['Protected',8,16,32,64],['Manageable',16,32,64,128]],{widths:[310,245,199,199,199],size:26,h:285});
text(s,'One boundary: above-cutoff Vital / Protected is Unsafe',64,500,1152,48,29,C.teal,true);
note(s,'Differentiated bowel cells and fibroblasts use Protected. Crypt progenitors use Vital. Equality passes.',568,25);
}
{
const s=slide('Normal-tissue hazard ontology','Bulk and cellular consequence','Central registry and full evidence ledger: https://pirl-unc.github.io/oncoref/audits/tissue-hazards-1.29.0/ . Adult human qualitative policy. AT2 evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC4821031/ . Club-cell depletion: https://pmc.ncbi.nlm.nih.gov/articles/PMC3056226/ . Ciliated-cell depletion: https://pmc.ncbi.nlm.nih.gov/articles/PMC4521991/ .');
table(s,[['Tier','Examples','Key interpretation'],['Vital','AT1/AT2, club/basal, crypt progenitors','Function and regenerative reserve'],['Protected','Ciliated/goblet, fibroblasts, skin','Persistent injury and uncertain control'],['Manageable','Thyroid, mature B cells','Defined compensation precedent'],['Disposable','Cutaneous pigment, hair follicles','Narrow anatomical allowance']],{widths:[230,455,467],size:25,h:350});
note(s,'Bulk retains organ concern. Recognition and reproductive membership are separate. Protected rare-cell findings create a review hold.',570,25);
}
{
const s=slide('Clinical evidence adjusts RNA ceilings','TCR-T evidence','One step for reviewed treated-cohort evidence without identified grade 4-5 specificity-related toxicity; cap at Enhanced. Unknown or suspected specificity injury blocks product credit. Exact peptide/HLA/receptor limits remain. Clinical trial and SAE review: https://pirl-unc.github.io/oncoref/audits/cta-tcr-trials-20261029/');
table(s,[['With clinical credit','Missing / unusable','Approved −','Supported −','Enhanced −'],['Vital',8,16,32,32],['Protected',16,32,64,64],['Manageable',32,64,128,128]],{widths:[310,245,199,199,199],size:25,h:285});
text(s,'IHC grades stay unchanged. Trials do not stack.',64,488,1152,48,29,C.teal,true);
note(s,'Protein and clinical-injury exclusions remain. This is protein-profile credit with tested epitope limitations, not clinical clearance.',550,25);
}
{
const s=slide('Serious events and antigen specificity','Clinical safety',Object.values(trials.sources).map(r=>`${r.title}: ${r.url}`).join('\n'));
table(s,[['Targeted CTA epitope','Important event','Specificity attribution'],['MAGE-A3 / a3a','Two fatal cardiac injuries','Titin cross-recognition implicated'],['MAGE-A3/A9 / NCI receptor','Three neurological cases; two deaths','MAGE-A12 / EPS8L2 attribution unresolved'],['NY-ESO-1 / LAGE-1A, c259','Late cardiac death; two GBS cases','Peptide-specific mechanism not established'],['NY-ESO-1 / TBI-1301','Grade 3 lung injury','Unresolved; no credit for this product'],['MAGE-A4 / afami-cel','Aplasia and stroke deaths in phase 1','No demonstrated antigen-specific mechanism'],['MAGE-A10 / ADP-A2M10','Grade 4 CRS; fatal aplasia after repeat treatment','Conditioning favored for marrow injury']],{widths:[355,410,387],size:21,h:405});
note(s,'Clinical trial and SAE review. Seriousness differs from severity; this is a selected summary, not a pooled event rate.',602,23);
}
{
const s=slide('Other records an unmet admission gate','Categorical decisions','Safety and coverage are assessed independently. Other reasons use reproductive presence, cancer presence, then concentration. No gene is admitted or excluded by an abundance label.');
const meanings={Nonreproductive:'No qualifying reproductive evidence',Reproductive:'No qualifying cancer evidence',Diffuse:'Reproductive concentration fails',Translation:'Selected primary translation requires review',Held:'Target gates pass; normal review hold remains',Nested:'Nested short reference-CDS product',Contained:'At least 95% of 9-mers occur in Unsafe proteins'};
table(s,[['Reason','Meaning','Proteins'],...Object.entries(data.other_subtypes).filter(([,v])=>v>0).map(([k,v])=>[k,meanings[k]||k,v])],{widths:[290,722,140],size:23,h:380});
note(s,'Normal RNA level, tumor support and Disposable expression remain separate attributes.',580,25);
}
{
const s=slide('Evidence attributes stay independent','Inspect each protein','The profiles() API preserves native RNA units, supporting samples, adult Disposable pattern and developmental evidence independently of the final decision.');
table(s,[['Question','Fields / interpretation'],['How high in cancer?','Native TPM / HPA pTPM maxima; two identified tumors >=10 TPM is a support flag'],['How high in adult normal?','Separate nTPM / TPM / CP10K summaries; raw unsupported observations remain visible'],['Adult Disposable pattern?','Permitted adult expression with low measured adult expression elsewhere'],['Developmental expression?','Separate lineage RNA and measured-stage coverage; never changes adult tissue permission']],{widths:[350,802],size:23,h:365});
note(s,'cta_selection.profiles() exposes these fields. Low / Weak / Normal are no longer competing classes.',570,25);
}
{
const s=slide('Reproductive cell evidence','Positive association','HPA v23 https://v23.proteinatlas.org/about/download ; passenger transcripts: https://pmc.ncbi.nlm.nih.gov/articles/PMC7459098/ ; normal interferon induction: https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/');
table(s,[['Route','Positive evidence','Additional requirements'],['Identified lineage','≥5 nTPM with ≥20 reads','Mapped reproductive lineage; all protein loci measured'],['Indirect site','Generic/immune cell in reproductive tissue','Complete outside cluster panel; every value <5 nTPM; no normal-induction conflict'],['Examples','CAGE1 208.9; PRDM9 61.8; FOXR2 5.2 nTPM','Spermatids / spermatocytes / spermatogonia'],['Counterexample','IFNA2 and IFNA14 in breast macrophages','Reviewed normal immune protein → Unsafe']],{widths:[275,415,462],size:23,h:375});
note(s,'Indirect support identifies the sample site, not the producing cell. Normal and bulk-concentration gates still apply.',575,24);
}
{
const s=slide('Reviewed normal-cell knowledge','Literature evidence','OncoRef 1.26.0: https://pirl-unc.github.io/oncoref/audits/cta-other-literature-20261030/ ; pDC protein: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full ; human olfactory atlas: https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/');
table(s,[['Evidence','Consequence'],['Established normal production','Vital/Protected protein or localized RNA plus orthogonal cell function can establish Unsafe'],['Unresolved normal findings','Sparse RNA, uncertain species/stage or isoforms remain review holds'],['Reproductive protein','Supplies presence below the RNA floor; safety and concentration still apply'],['Search and curation','All 611 baseline Other genes searched; accepted findings remain separate from search leads']],{widths:[340,812],size:25,h:370});
note(s,'Rare normal cells and activation states count. Expression exclusions do not assert clinical injury.',590,25);
}
{
const s=slide('Historical findings retain explicit dispositions','Evidence migration','All historical findings and prior adjudications: https://pirl-unc.github.io/oncoref/audits/cta-evidence-migration-20261104/ . MAGE-A12/EPS8L2 counterevidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC7946352/ .');
table(s,[['Case','Current result','Reason retained'],['MAGEC2','Other:Held','Purkinje-cell reagent specificity remains unresolved'],['LGALS13','Other:Held','Adult/fetal and galectin identity remain unresolved'],['PLCZ1','Other:Held','Pancreatic RNA is unresolved; CFTR staining is not PLCZ1 protein'],['MAGEA12','Unsafe','Independent rare-cell exclusion; clinical peptide attribution remains unresolved']],{widths:[220,245,687],size:24,h:360});
note(s,'A completed review can retain uncertainty. Unresolved evidence requires review without establishing Unsafe.',590,24);
}
{
const s=slide('Every candidate receives a classification','Complete inventory','Full record-level audit and historical omission ledger: https://pirl-unc.github.io/oncoref/audits/cta-inventory-20261105/ . Qualifying coding annotation is not experimental proof of protein production.');
table(s,[['Candidate records','Count','Current treatment'],['Previously classified',n(inventory.input_records-inventory.restored_records),'All previous classifications unchanged'],['Restored: Other:Translation',n(inventory.restored_classifications['Other:Translation']),'Translation admission fails; sequence review required'],['Restored: Unsafe',n(inventory.restored_classifications.Unsafe),'Normal exclusion takes precedence'],['Complete inventory',n(inventory.input_records),'All records enter the same assessment']],{widths:[355,145,652],size:25,h:360});
note(s,'Translation is a CTA admission gate. Insufficient support cannot hide normal-tissue evidence.',590,25);
}
{
const q=inventory.preliminary_expression_profile,z=inventory.stored_zero_rna_screen;
const s=slide('Expression-screen matches and counterevidence','Cancer restriction','Record-level matches, review statuses, counterevidence and stored-zero RNA audit: https://pirl-unc.github.io/oncoref/audits/cta-inventory-20261105/ . Normal findings for IFNW1, IFNA4, MRGPRX4, CNGA2, SCN10A and GFRAL are retained with original sources.');
table(s,[['Numerical screen','Count','Interpretation'],['Low normal expression profile',n(q.matches),'Preliminary matches; every independent exclusion still applies'],['Normal-cell counterevidence',n(q.disposition_counts.excluded_by_normal_evidence),'Unsafe despite a routine-panel profile match'],['Stored adult/reproductive RNA = 0',n(z.adult_and_reproductive_bulk_zero),'Cancer present and required normal bulk coverage complete'],['Also developmental RNA = 0',n(z.also_developmental_zero),`${n(z.developmental_measurement_missing)} adult-zero records lack developmental measurements`]],{widths:[425,115,612],size:24,h:370});
note(s,'Stored RNA zeros do not establish protein absence. CTA counts are not cancer-only protein counts.',590,24);
}
{
const s=slide('Selected translations require qualification','Reference sequence','Exact reference records and qualification rules: https://pirl-unc.github.io/oncoref/audits/translation-qualification-20261104/ . Conservative AA20 support does not establish absent biological production.');
const names=['SERPINA2','MROH5','SELENOV'];
table(s,[['Record','Current classification','Selected translation'],...names.map(name=>{const r=definition.worked.find(x=>x.label===name);return [name,r.classification,name==='SERPINA2'?'Incomplete primary LoF record; alternate locus cannot rescue it':name==='MROH5'?'Interrupted primary record; normal exclusion still takes precedence':'Selenocysteine requires review under the AA20 contract']})],{widths:[220,330,602],size:25,h:310});
note(s,'Qualification reviews the selected reference sequence. Original interrupted and alternate records remain available.',560,26);
}
{
const s=slide('Sequence roles preserve every occurrence','Verified sequence bundle','Verified loader, coverage gaps and historical transition: https://pirl-unc.github.io/oncoref/sequence-bundle/ . Bundle identity: '+sequences.bundle_id);
table(s,[['Role','Source records','Meaning'],['Exclude',n(sequences.record_roles.exclude),'Applicable normal exclusion takes precedence'],['Review',n(sequences.record_roles.review),'Unresolved decision, isoform or translation'],['Allow',n(sequences.record_roles.allow),'Qualified assessed representative under this policy']],{widths:[230,250,672],size:26,h:300});
note(s,`${n(sequences.record_count)} original records and all source occurrences survive. AA20 segments retain coordinates and never bridge ambiguity.`,555,26);
}
{
const s=slide('High tumor expression does not establish CTA identity','Four requested proteins','Pinned normal atlas peaks are source-supported summaries, not exhaustive rare-cell or donor maxima. OR10A4 also has normal human nasal RNA in two primary studies.');
table(s,[['Protein','Cancer TPM','Adult atlas nTPM','Current reasoning'],['DBX1','297.5','0.8; thymic T cells','No qualifying reproductive evidence'],['OR2B2','246.3','0.5; pons','No qualifying reproductive evidence'],['KRTAP20-2','147.6','1.4; skin','Skin expression; reproductive identity unestablished'],['OR10A4','138.0','1.9; basal ganglia','Normal nasal RNA requires review']],{widths:[220,175,335,422],size:24,h:360});
note(s,'Disposable describes allowed adult tissues. Developmental expression is a separate measurement.',575,25);
}
{
const s=slide('Other accumulates as identity gates fail','Coarse category funnel','Normal assessment is displayed first. Identity requirements move proteins into Other without creating Unsafe evidence. All final categories match the released classifier.');
figure(s,'category-funnel',64,158,1152,475);
}
{
const s=slide('Every filtering layer can contribute','Human funnel','The gray segment contains proteins pending final assignment. All layers evaluate every protein. Final human categories contain no unassigned proteins.');
const shortLabels={'Reproductive / embryonic association':'Reproductive evidence','Cancer expression':'Cancer expression','HPA normal protein':'HPA protein','Required RNA coverage':'Coverage','Reproductive concentration':'Concentration','Normal bulk':'Bulk RNA','HPA tissue clusters':'Tissue clusters','Monaco immune':'Monaco','DICE immune':'DICE','Adult cornea':'Cornea','Accepted normal-protein evidence':'Curated protein','Reviewed normal-cell knowledge':'Normal-cell knowledge','Clinical safety exclusion':'Clinical safety','External validated normal protein':'External IHC','Rare vital-tissue cell expression':'Rare Vital cells','Rare Protected-cell expression review':'Rare Protected cells','CTA and Other assignments':'Final categories'};
const labels=data.funnel.map(r=>shortLabels[r.layer]||r.layer);
chart(s,labels,[{name:'Eligible so far',values:data.funnel.map(r=>r.CTA_eligible_so_far),fill:C.teal},{name:'Unsafe',values:data.funnel.map(r=>r.Unsafe),fill:C.rust},{name:'Other',values:data.funnel.map(r=>r.Other),fill:'#A28C68'},{name:'Pending',values:data.funnel.map(r=>r.Unassigned||0),fill:'#B9C3CD'}],{h:400,y:158,stacked:true,size:18});
[['Eligible so far',C.teal],['Unsafe',C.rust],['Other','#A28C68'],['Pending','#8795A5']].forEach(([label,color],i)=>text(s,label,330+i*210,586,205,28,20,color,true));
note(s,`Final: ${n(counts.CTA)} CTA / ${n(counts.Unsafe)} Unsafe / ${n(counts.Other)} Other.`,625,21);
}
for (let page=0;page<2;page++){
const s=slide(`Removals at each filter: ${page+1} of 2`,'Human funnel','Sequential first removals depend on the displayed order; the independent-layer table reports overlaps.');
table(s,[['Layer','First removals','Eligible','Unsafe'],...data.funnel.slice(page*10,(page+1)*10).map(r=>[r.layer.replace('Reproductive / embryonic association','Association'),n(r.new_candidate_removals),n(r.CTA_eligible_so_far),n(r.Unsafe)])],{widths:[575,190,190,197],size:22,h:380,y:170});
note(s,'Sequential removals overlap. Final categories use every applicable evidence layer.',612,23);
}
{
const s=slide('IHC matching to measured RNA compartments','Anatomy and evidence','HPA colon https://www.proteinatlas.org/learn/dictionary/normal/colon ; bronchus https://www.proteinatlas.org/learn/dictionary/normal/bronchus . Parent-compartment credit is a declared anatomical interpretation, not proof that all rare cells were sampled.');
table(s,[['Context','Evidence and current ceiling'],['PRAME colon undifferentiated cells','20.3 nTPM / 96 reads; Enhanced-negative glandular IHC; Vital 32'],['PRAME adrenal bulk','3.6 nTPM; incomplete IHC coverage; Vital 4 → 8 with IMA203 credit'],['CTAG2 bronchial basal cells','24.9 nTPM / 58 reads; Enhanced-negative respiratory IHC; Vital 32'],['Unrelated cells in the same tissue','No negative-IHC credit without exact or declared compartment coverage']],{widths:[365,787],size:25,h:350});
note(s,'PRAME and CTAG2 pass this policy. Parent-compartment sensitivity and PRAME’s Low adrenal protein remain explicit limitations.',590,23);
}
for(const [name,title] of [['selected-reproductive-embryonic','Selected CTA: adult/placental versus embryonic RNA'],['selected-reproductive-embryonic-immune','Selected CTA: healthy RNA in TPM'],['selected-reproductive-embryonic-cornea','Selected CTA: normal corneal RNA']]){
const s=slide(title,'Expression distributions','Color uses the highest supported non-disposable healthy expression within the stated native unit. There is no mixed-unit maximum. Labels include widely studied antigens with mapped values. Companion panels label all selected proteins ≥5 on both axes. Missing embryonic measurements appear in a separate categorical band; every selected protein is represented without zero filling.');
figure(s,name+'-slide',48,160,1184,475);
text(s,'Gray: healthy measurement unavailable. The lower band has no embryonic-expression value.',64,639,1152,25,16,C.muted);
}
for(const number of [1,2,3]){
const s=slide(`High adult/placental and embryonic RNA: labels ${number} of 3`,'Selected CTA detail','The three panels use identical axes and together label all selected proteins with reproductive nTPM ≥5 and measured embryonic lineage-mean TPM ≥5. This descriptive threshold does not replace developmental cell-support requirements. Color is the highest non-disposable healthy nTPM.');
figure(s,`selected-high-both-${number}-slide`,48,160,1184,475);
}
{
const s=slide('Complementary embryonic measurements','Developmental coverage','Tyser E-MTAB-9388 https://doi.org/10.1038/s41586-021-04158-y . Petropoulos E-MTAB-3929 https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results . Sources retain their own stage, cell support and donor counts.');
table(s,[['Source','Measured cells / embryos','Developmental scope'],['Tyser et al. (2021)','1,195 cells / 1 embryo','Gastrulation, Carnegie stage 7'],['Petropoulos et al. (2016)','1,519 cells / 88 embryos','Preimplantation, days 3–7'],['Recovered coverage','683 of 1,006 earlier missing proteins','Measured TPM, including familiar CTAs']],{widths:[360,420,372],size:26,h:320});
text(s,'PRAME, MAGE-A4, CTAG1A/B, XAGE1A/B and CT83 now have TPM rows',64,525,1152,70,28,C.teal,true);
note(s,'Each source keeps its own stages. Neither atlas covers all fetal or pediatric development.',610,23);
}
for(const [name,title] of [['reproductive','Widely studied antigens: adult/placental RNA'],['embryonic','Widely studied antigens: embryonic RNA']]){
const s=slide(title,'Current assignments','MAGE-A1 and MAGE-A3 are comparators, not admission controls. Full-width figures highlight PRAME, NY-ESO-1 (CTAG1A/B), XAGE1A/B, MAGE-A4, CT83, SSX2/B and TBXT. Missing embryonic measurements are in a separate band.');
figure(s,'widely-studied-antigens-'+name+'-slide',48,160,1184,475);
}
{
const s=slide('Assessment status is separate from category','Reading the results','A category can be assigned despite gaps in optional or irrelevant sources. Screened does not mean clinically safe or fully measured. MAGE-A3 has measured colon counterevidence. An independent developmental measurement cannot erase it.');
table(s,[['API status','Meaning','Human proteins'],['screened','Classification with a resolved applicable screen',n(data.assessment_statuses.screened)],['incomplete_evidence','Classified with disclosed evidence gaps',n(data.assessment_statuses.incomplete_evidence)],['normal_evidence_hold','Required normal coverage passes; a normal-screen hold remains',n(data.assessment_statuses.normal_evidence_hold)],['translation_review','Selected primary translation needs qualification',n(data.assessment_statuses.translation_review)],['unassessed','No species-specific selection policy applied','0']],{widths:[360,622,170],size:22,h:355});
note(s,'Other includes review holds and sequence rules. No human protein has a null category.',585,26);
}
{
const s=slide('Which sources have incomplete evidence','Coverage audit','Source versions, complete missing tissue/cell contexts, member ENSG IDs, historical IDs, synonyms and tumor cohort coverage: '+URL+'coverage/ . Counts overlap and do not sum to a total.');
table(s,[['Source or unresolved evidence','Proteins'],...data.normal_gap_counts.map(r=>[r.source,r.groups]),['Combined Tyser + Petropoulos TPM',data.evidence_gap_counts.embryonic_measurement_missing],['No tumor measurement: OncoRef + HPA',data.evidence_gap_counts.tumor_measurement_missing],['Measured tumor, incomplete member coverage',data.evidence_gap_counts.tumor_group_profile_incomplete],['HPA v23 Uncertain normal-tissue IHC',data.evidence_gap_counts.unresolved_evidence]],{widths:[970,182],size:22,h:395,y:170});
note(s,`Overlapping gaps within ${n(data.assessment_statuses.incomplete_evidence)} proteins with incomplete evidence. DICE and corneal coverage are optional.`,600,24);
}
{
const s=slide('Dog evidence and integration gaps','Species-specific evidence','Nemec et al. (2019) https://doi.org/10.1111/vco.12475 ; Hindriks et al. (2025) https://doi.org/10.1016/j.vetimm.2025.110946 . Issues https://github.com/pirl-unc/oncoref/issues/574 and https://github.com/pirl-unc/oncoref/issues/585 .');
table(s,[['Layer','Available evidence / remaining gap'],['Nominations','5 + 10 − 1 shared = 14 native nominations'],['Identity','13 unique mappings; one ambiguous MAGE-family assay'],['Normal / developmental RNA','Bgee TPM and BarkBase TMM-CPM need reference reconciliation'],['Tumor RNA','GSE304066 count/FPKM and reference conflicts remain unresolved'],['Normal protein','Limited study IHC; Pfizer pooled iBAQ is not cell-resolved IHC'],['Selection','All 14 unassigned / unassessed; no canine quantitative policy']],{widths:[290,862],size:24,h:390});
note(s,'Human thresholds, antibody grades and normal-tissue absence never transfer by orthology.',605,24);
}
{
const s=slide('Curated oncofetal sources','Developmental nominations','Shaw et al. (2024) https://doi.org/10.1038/s41467-024-47649-y ; TRIM71 liver study https://pmc.ncbi.nlm.nih.gov/articles/PMC11388079/ ; Cardoso-Moreira et al. (2019) https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/ . Newly identified sources are not in current selection counts.');
table(s,[['Source','Scope and integration status'],['Tyser + Petropoulos','Current gastrulation and preimplantation TPM sources'],['Shaw et al. (2024)','157 cancer-specific targets; six oncofetal annotations. Exon-specific scope retained'],['TRIM71 liver study (2024)','Liver oncofetal signature; requires identity and measurement ingestion'],['Cardoso-Moreira et al. (2019)','Developmental organ trajectories with child/adult comparisons; not integrated']],{widths:[380,772],size:25,h:350});
note(s,'Published nominations need independent normal/tumor filtering and postnatal coverage.',590,26);
}
{
const s=slide('Common API with native evidence underneath','API','Human, dog and mouse share schema and explicit unavailable states. Data are offline and checksum-verified. Frozen release reports remain unchanged.');
text(s,'results = cta_selection.results(species="human")\naudit = cta_selection.audit(species="human")\naudit.sources    audit.overlaps    audit.funnel\naudit.explain("PRAME")',64,180,1152,200,31,C.teal,true);
table(s,[['Field','Meaning'],['category','CTA / Unsafe / Other; classification gives the named category'],['assessment_status','Screened, evidence/translation review, or unassessed'],['classification / subtype','One categorical label; Other subtype always populated']],{y:400,h:190,widths:[335,817],size:25});
note(s,'Source-specific native units, tissue/cell context, confidence and coverage remain inspectable.',610,24);
}
{
const s=slide('Every source resolves to a concrete version','Data API','Data version API: https://pirl-unc.github.io/oncoref/data-versions/ . GTEx versions and pipelines: https://github.com/broadinstitute/gtex-pipeline/blob/master/rnaseq/README.md . HPA current RNA snapshot: https://github.com/pirl-unc/oncoref/releases/tag/reference-hpa-v25.1-20261028 .');
text(s,'catalog.ensure("gtex_tissue_rna")               # v11\ncatalog.ensure("gtex_tissue_rna", version="v10")\ncatalog.provenance("gtex_tissue_rna", "v11")',64,173,1152,145,28,C.teal,true);
table(s,[['Shared choice','Source-specific meaning'],['Default / latest','Newest registered release; saved runs keep explicit pins'],['Concrete version','GTEx v11; HPA v25.1; a study snapshot or OncoRef build'],['Provenance','Version, content checksum, native unit, annotation and pipeline'],['Species','The same selector applies to native dog and mouse assets']],{y:340,h:235,widths:[310,842],size:24});
note(s,'This analysis pins GTEx v11 and HPA v23 normal panels. A newer public default does not rewrite it.',605,23);
}
for(const [i,s] of p.slides.items.entries()){
 await fs.writeFile(path.join(build,'rendered',`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await(await p.export({slide:s,format:'png',scale:1.5})).arrayBuffer()));
}
await fs.writeFile(path.join(build,'manifest.json'),JSON.stringify({policy:data.summary.policy,commit:'report-sha256:'+data.report_sha256,slideSize:{width:1280,height:720},pages},null,2));
const candidate=path.join(build,'candidate.pptx');
await(await PresentationFile.exportPptx(p)).save(candidate);
const finalPath=path.join(audit,'slides',process.env.CTA_SLIDE_FILENAME||`cta-evidence-walkthrough-${release}-validated-final.pptx`);
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath,pythonExecutable:process.env.RUNTIME_PYTHON,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tableOwners.flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:tableOwners,requiredNativeChartOwnerSlides:chartOwners,materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,receiptPath:path.join(build,'validation.json')});
console.log(JSON.stringify({slides:pages.length,finalPath,font,tableOwners,chartOwners}));
