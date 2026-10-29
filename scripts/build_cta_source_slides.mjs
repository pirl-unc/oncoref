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
const build=process.env.CTA_SLIDE_BUILD_DIR||path.join(root,'tmp/cta-source-slides-expanded');
const data=JSON.parse(await fs.readFile(path.join(audit,'presentation-data.json'),'utf8'));
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
await fs.mkdir(path.join(build,'rendered'),{recursive:true});
await fs.mkdir(path.join(audit,'slides'),{recursive:true});
function text(s,t,x,y,w,h,size=27,color=C.ink,bold=false){
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=String(t);sh.text.style={typeface:font,fontSize:size,color,bold};sh.text.insets=0;current.text.push(String(t));return sh;
}
function slide(title,section,notes=''){
 const s=p.slides.add();s.background.fill=C.white;current={title,text:[],links:[{url:URL,x:64,y:670,w:1100,h:24}]};pages.push(current);
 text(s,section.toUpperCase(),64,28,1152,25,16,C.teal,true);
 text(s,title,64,72,1152,85,42,C.ink,true);
 text(s,'OncoRef CTA evidence walkthrough   •   28 October 2026',64,670,1100,26,16,C.muted);
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
 values.forEach((row,r)=>{t.rows[r].height=h/values.length;row.forEach((_,c)=>{const cell=t.getCell(r,c);cell.fill=r===0?C.ink:(r%2===0?C.pale:C.white);cell.text.style={typeface:font,fontSize:size,color:r===0?C.white:C.ink,bold:r===0||c===0};});current.text.push(row.join(' | '));});return t;
}
function chart(s,categories,series,{x=64,y=170,w=1152,h=390,stacked=false,size=22,labels=false}={}){
 chartOwners.push(pages.length);
 const ch=s.charts.add('bar',{position:{left:x,top:y,width:w,height:h},categories:[...categories].reverse(),series:series.map(a=>({...a,values:[...a.values].reverse()})),barOptions:{direction:'bar',grouping:stacked?'stacked':'clustered',gapWidth:40},hasLegend:stacked,legend:{position:'bottom',textStyle:{fontSize:20,fill:C.ink}},xAxis:{textStyle:{fontSize:size,fill:C.ink},majorGridlines:null},yAxis:{textStyle:{fontSize:18,fill:C.muted}},dataLabels:{showValue:labels,position:'outEnd',textStyle:{fontSize:22,fill:C.ink}},chartFill:C.white,plotAreaFill:C.white});
 applyPresentationChartFont(ch,{fontFamily:font});return ch;
}
function figure(s,name,x,y,w,h){s.images.add({blob:new Uint8Array(readFileSync(path.join(audit,'figures',name+'.png'))),contentType:'image/png',fit:'contain',position:{left:x,top:y,width:w,height:h},alt:name.replaceAll('-',' ')});}
const n=v=>Number(v).toLocaleString('en-US');


const counts=data.summary.categories;
{
const s=slide('Human and dog CTA evidence','Current selection','One RNA ceiling. Neither requires an affirmative expression profile. Unassigned is an assessment status, not a biological category.');
text(s,`${n(counts.CTA)} human CTA groups`,64,202,1152,100,66,C.teal,true);
text(s,`${n(data.summary.selected_member_loci)} loci   /   ${n(counts.Unsafe)} Unsafe   /   ${n(counts.Neither)} Neither`,64,335,1152,65,32);
text(s,`${n(counts.Unassigned)} unassigned; dog selection remains unassessed`,64,446,1152,95,32);
note(s,'Sources, expression profiles, single-cutoff filtering and measured tissue evidence.',580,26);
}
{
const s=slide('Published candidate lists','Human annotations',data.papers.map(r=>`${r.citation}: ${r.source_url}`).join('\n'));
table(s,[['Paper','Groups','CTA','Unsafe','Neither','Unassigned'],...data.papers.map(r=>[r.citation,n(r.mapped_entities),r.CTA,r.Unsafe,r.Neither,r.Unassigned])],{widths:[390,150,140,160,150,162],size:19,h:420,y:158});
note(s,`Union: ${n(data.paper_union_groups)} groups / ${n(data.paper_union_loci)} loci. ${data.annotated} selected groups have an annotation; ${data.unannotated} do not.`,600,23);
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
const s=slide('Reproductive and developmental association','Positive expression contexts','Tyser et al. (2021): https://doi.org/10.1038/s41586-021-04158-y . HPA v23 https://v23.proteinatlas.org/about/download .');
table(s,[['Qualifying route','Entrants','Final CTA'],...data.routes.map(r=>[r.route,n(r.entrants),n(r.CTA)])],{widths:[620,266,266],h:280,size:29});
text(s,'Whole-reference RNA measurements create this set',64,495,1152,50,32,C.teal,true);
note(s,'This is not a 12,000-gene literature list. Embryonic presence does not establish developmental restriction.',565,26);
}
{
const s=slide('Normal and cancer measurement sources','Human evidence','HPA https://v23.proteinatlas.org/about/download ; DICE https://dice-database.org/ ; cornea https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683 ; Tyser https://doi.org/10.1038/s41586-021-04158-y . Full URLs and sample lineage accompany the report.');
table(s,[['Evidence','Source','Quantity'],['Normal / reproductive bulk','HPA consensus + HPA brain + direct GTEx','nTPM / TPM'],['Normal tissue cells','HPA clusters, at least 20 reads','nTPM'],['Immune cells','HPA / Monaco / original DICE','Native nTPM / TPM'],['Adult cornea','Collin et al., GSE155683','CP10K'],['Normal protein','HPA + curated external IHC','Staining / reagent confidence'],['Cancer','132 compatible OncoRef cohort views + HPA','TPM / protein-coding TPM'],['Developmental','Tyser gastrulation + Petropoulos preimplantation','Stage / lineage mean TPM']],{widths:[310,570,272],size:22,h:405});
note(s,'DICE and corneal coverage are optional. Measured counterevidence still contributes.',600,24);
}
{
const s=slide('RNA cutoffs','Current policy','The same nominal values apply across native TPM/nTPM assays. CP10K uses nominal /100, not a TPM conversion. Values are screening-policy choices, not measured toxicity limits.');
table(s,[['Ceiling','Missing / unusable','Approved −','Supported −','Enhanced −'],['Vital',4,8,16,32],['Protected',8,16,32,64],['Manageable',16,32,64,128]],{widths:[310,245,199,199,199],size:26,h:285});
text(s,'One boundary: above-cutoff Vital / Protected is Unsafe',64,500,1152,48,29,C.teal,true);
note(s,'Manageable excess blocks CTA. Bowel epithelium and identified fibroblasts use Manageable. Equality passes.',568,25);
}
{
const s=slide('Neither requires a measured expression profile','Category definitions','Low means below the common 2-unit expression floor in every measured source, with required normal panels and tumor/developmental measurements. Missing optional DICE/cornea is disclosed. No detected protein is not proof of absence in every unmeasured tissue.');
table(s,[['Profile','Required evidence','Groups'],['Low expression','RNA below 2 across measured normal / reproductive / embryonic / tumor panels; no detected protein',data.neither_profiles.low_expression||0],['Diffuse Manageable','No association route; RNA ≥2 in at least two Manageable tissue concepts; no Vital / Protected exclusion',data.neither_profiles.diffuse_manageable||0]],{widths:[265,747,140],size:24,h:300});
note(s,'Both require complete normal RNA coverage and resolved evidence. Other failures remain unassigned.',550,28);
}
{
const s=slide('Every filtering layer can contribute','Human funnel','The gray segment is held/unassigned evidence, not Neither. The final profile check assigns Neither. All layers evaluate every group; final categories are independent of display order.');
const labels=['Association','Tumor','HPA protein','Coverage','Bulk RNA','Tissue clusters','HPA immune','Monaco','DICE','Cornea','Curated protein','Clinical safety','External IHC','Neither profiles'];
chart(s,labels,[{name:'Eligible so far',values:data.funnel.map(r=>r.CTA_eligible_so_far),fill:C.teal},{name:'Unsafe',values:data.funnel.map(r=>r.Unsafe),fill:C.rust},{name:'Neither',values:data.funnel.map(r=>r.Neither),fill:'#A28C68'},{name:'Unassigned',values:data.funnel.map(r=>r.Unassigned),fill:'#B9C3CD'}],{h:435,y:158,stacked:true,size:20});
note(s,`Final: ${n(counts.CTA)} CTA / ${n(counts.Unsafe)} Unsafe / ${n(counts.Neither)} Neither / ${n(counts.Unassigned)} unassigned.`,610,24);
}
{
const s=slide('Removals at each filter','Human funnel','Sequential first removals depend on the displayed order; the independent-layer table reports overlaps.');
table(s,[['Layer','First removals','Eligible','Unsafe'],...data.funnel.map(r=>[r.layer.replace('Reproductive / embryonic association','Association'),n(r.new_candidate_removals),n(r.CTA_eligible_so_far),n(r.Unsafe)])],{widths:[575,190,190,197],size:18,h:430,y:158});
note(s,'Every group remains available to every later evidence layer.',612,23);
}
{
const s=slide('IHC matching to measured RNA compartments','Anatomy and evidence','HPA colon https://www.proteinatlas.org/learn/dictionary/normal/colon ; bronchus https://www.proteinatlas.org/learn/dictionary/normal/bronchus . Parent-compartment credit is a declared anatomical interpretation, not proof that all rare cells were sampled.');
table(s,[['Context','Evidence and current ceiling'],['PRAME colon epithelial cells','20.3 nTPM / 96 reads; Enhanced-negative glandular IHC; Manageable 128'],['PRAME adrenal bulk','3.6 nTPM; incomplete compartment IHC coverage; Vital 4'],['CTAG2 bronchial basal cells','24.9 nTPM / 58 reads; Enhanced-negative respiratory IHC; Vital 32'],['Unrelated cells in the same tissue','No negative-IHC credit without exact or declared compartment coverage']],{widths:[365,787],size:25,h:350});
note(s,'PRAME and CTAG2 pass this policy. Parent-compartment sensitivity and PRAME’s Low adrenal protein remain explicit limitations.',590,23);
}
for(const [name,title] of [['selected-reproductive-embryonic','Selected CTA: reproductive versus embryonic RNA'],['selected-reproductive-embryonic-immune','Selected CTA: healthy RNA in TPM'],['selected-reproductive-embryonic-cornea','Selected CTA: normal corneal RNA']]){
const s=slide(title,'Expression distributions','Color uses the highest supported non-disposable healthy expression within the stated native unit. There is no mixed-unit maximum. Labels include widely studied antigens with mapped values. Companion panels label all selected groups ≥5 on both axes. Missing embryonic measurements appear in a separate categorical band; every selected group is represented without zero filling.');
figure(s,name+'-slide',48,160,1184,475);
text(s,'Gray: healthy measurement unavailable. The lower band has no embryonic-expression value.',64,639,1152,25,16,C.muted);
}
for(const number of [1,2,3]){
const s=slide(`High reproductive and embryonic RNA: labels ${number} of 3`,'Selected CTA detail','The three panels use identical axes and together label all selected groups with reproductive nTPM ≥5 and measured embryonic lineage-mean TPM ≥5. This descriptive threshold does not replace developmental cell-support requirements. Color is the highest non-disposable healthy nTPM.');
figure(s,`selected-high-both-${number}-slide`,48,160,1184,475);
}
{
const s=slide('Complementary embryonic measurements','Developmental coverage','Tyser E-MTAB-9388 https://doi.org/10.1038/s41586-021-04158-y . Petropoulos E-MTAB-3929 https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results . Sources retain their own stage, cell support and donor counts.');
table(s,[['Source','Measured cells / embryos','Developmental scope'],['Tyser et al. (2021)','1,195 cells / 1 embryo','Gastrulation, Carnegie stage 7'],['Petropoulos et al. (2016)','1,519 cells / 88 embryos','Preimplantation, days 3–7'],['Recovered coverage','683 of 1,006 earlier missing groups','Measured TPM, including familiar CTAs']],{widths:[360,420,372],size:26,h:320});
text(s,'PRAME, MAGE-A4, CTAG1A/B, XAGE1A/B and CT83 now have TPM rows',64,525,1152,70,28,C.teal,true);
note(s,'Each source keeps its own stages. Neither atlas covers all fetal or pediatric development.',610,23);
}
for(const [name,title] of [['reproductive','Widely studied antigens: reproductive RNA'],['embryonic','Widely studied antigens: embryonic RNA']]){
const s=slide(title,'Current assignments','MAGE-A1 and MAGE-A3 are comparators, not admission controls. Full-width figures highlight PRAME, NY-ESO-1 (CTAG1A/B), XAGE1A/B, MAGE-A4, CT83, SSX2/B and TBXT. Missing embryonic measurements are in a separate band.');
figure(s,'widely-studied-antigens-'+name+'-slide',48,160,1184,475);
}
{
const s=slide('Assessment status is separate from category','Reading the results','A category can be assigned despite gaps in optional or irrelevant sources. Screened does not mean clinically safe or fully measured. MAGE-A3 has measured colon counterevidence. An independent developmental measurement cannot erase it.');
table(s,[['API status','Meaning','Human groups'],['screened','CTA, Unsafe or Neither assigned',n(counts.CTA+counts.Unsafe+counts.Neither)],['incomplete_evidence','Unassigned with missing coverage or unresolved evidence',n(data.assessment_statuses.incomplete_evidence)],['outside_category_definitions','Tracked coverage checks pass; no definition fits',n(data.assessment_statuses.outside_category_definitions)],['unassessed','No species-specific selection policy applied','0']],{widths:[360,622,170],size:24,h:335});
note(s,'Neither requires a measured biological profile. Unassigned does not mean Unsafe or low expression.',585,26);
}
{
const s=slide('Which sources have incomplete evidence','Coverage audit','Source versions, complete missing tissue/cell contexts, member ENSG IDs, historical IDs, synonyms and tumor cohort coverage: '+URL+'coverage/ . Counts overlap and do not sum to a total.');
table(s,[['Source or unresolved evidence','Groups'],...data.normal_gap_counts.map(r=>[r.source,r.groups]),['Combined Tyser + Petropoulos TPM',data.evidence_gap_counts.embryonic_measurement_missing],['No tumor measurement: OncoRef + HPA',data.evidence_gap_counts.tumor_measurement_missing],['Measured tumor, incomplete group profile',data.evidence_gap_counts.tumor_group_profile_incomplete],['HPA v23 Uncertain normal-tissue IHC',data.evidence_gap_counts.unresolved_evidence]],{widths:[970,182],size:22,h:395,y:170});
note(s,`Overlapping gaps within ${n(data.assessment_statuses.incomplete_evidence)} unassigned groups. DICE and corneal coverage are optional.`,600,24);
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
table(s,[['Field','Meaning'],['category','CTA / Unsafe / Neither, or null when unassigned'],['assessment_status','Screened, incomplete evidence, outside definitions, or unassessed'],['gene_ids / gene_identities','All group members, stable IDs, history and curated synonyms']],{y:400,h:190,widths:[335,817],size:25});
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
const finalPath=path.join(audit,'slides',process.env.CTA_SLIDE_FILENAME||'cta-evidence-walkthrough-expanded.pptx');
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath,pythonExecutable:process.env.RUNTIME_PYTHON,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tableOwners.flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:tableOwners,requiredNativeChartOwnerSlides:chartOwners,materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,receiptPath:path.join(build,'validation.json')});
console.log(JSON.stringify({slides:pages.length,finalPath,font,tableOwners,chartOwners}));
