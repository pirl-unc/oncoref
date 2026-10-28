// Editable release deck from presentation-data.json; scientific scatter from the audit.
// Required: RUNTIME_NODE_MODULES, ONCOREF_PRESENTATIONS_SKILL, RUNTIME_PYTHON.
// Run prepare_cta_release_presentation.py first. Output goes to the versioned audit.
import fs from 'node:fs/promises';
import {readFileSync} from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
import {fileURLToPath, pathToFileURL} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const audit=path.join(root,'docs/audits/cta-release-1.21.0');
const build=path.join(root,'tmp/cta-release-slides');
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
const URL='https://pirl-unc.github.io/oncoref/audits/cta-release-1.21.0/';
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
 text(s,'OncoRef 1.21.0 · Current research policy · Full evidence and downloads in the release report',64,670,1100,26,16,C.muted);
 text(s,String(pages.length).padStart(2,'0'),1174,668,42,30,19,C.muted,true);
 s.speakerNotes.text=`${title}\n${notes}\nReport: ${URL}\nPolicy: ${data.summary.policy}\nReport manifest SHA256: ${data.report_sha256}`;
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
const s0=slide('CTA and oncofetal selection','Release 1.21.0','All counts are recomputed from the complete reference; no clinical safety claim.');
text(s0,n(data.summary.categories.CTA),64,205,520,160,126,C.teal,true);
text(s0,'selected protein-equivalence groups',64,390,1100,60,36,C.ink,true);
text(s0,`${n(data.summary.selected_member_loci)} gene loci · ${n(data.summary.categories.Unsafe)} Unsafe · ${n(data.summary.categories.Neither)} Neither`,64,475,1152,70,31);
note(s0,'Versioned evidence, uniform rules and explicit coverage for human, dog and mouse.',580,28);
await fs.writeFile(path.join(build,'rendered/slide-01.png'),new Uint8Array(await(await p.export({slide:s0,format:'png',scale:1.5})).arrayBuffer()));

{
const s=slide('The selection rules','Human policy','Reproductive HPA nTPM; Tyser E-MTAB-9388 gene TPM. Embryonic support: lineage >=20 cells, >=5 positive cells and >=10% at >=1 TPM. Candidate equality passes; stronger RNA evidence >=2x. Every layer contributes even after a hold.');
table(s,[['Gate','Rule'],['Association','Reproductive ≥5 nTPM and > adult bulk; OR embryonic mean ≥5 TPM'],['Tumor','≥2 TPM / pTPM in at least one compatible tumor sample'],['Coverage','Required bulk / cell / immune panels; DICE and cornea optional'],['Normal RNA','Tissue / IHC ceiling; > ceiling is a hold, ≥2× is Unsafe evidence'],['Normal protein','Low stays flagged; stronger / unquantified accepted positives exclude'],['Final category','CTA passes all gates; Unsafe dominates; other failures are Neither']],{widths:[245,907],size:24,h:390});
note(s,'Thresholds express research preferences. No named gene or clinical product bypasses a gate.');
}
{
const s=slide('One linked threshold grid','Normal RNA','Native TPM/nTPM retain separate assay identity. CP10K uses nominal ceiling/100, not an empirical conversion. Negative IHC requires exact tissue/cell match, complete group coverage, positive profile controls for every member and no unresolved conflict. https://v23.proteinatlas.org/about/assays%2Bannotation');
table(s,[['Candidate ceiling','Missing / unusable','Approved −','Supported −','Enhanced −'],['Vital',4,8,16,32],['Protected',8,16,32,64],['Manageable',16,32,64,128],['Specified bowel cells',32,64,128,256]],{widths:[340,235,192,192,193],size:26,h:310});
text(s,'Stronger RNA boundary = 2 × candidate ceiling',64,510,1152,45,31,C.teal,true);
note(s,'Bulk keeps its organ tier. Mapped adult HPA fibroblasts use Manageable; specified colon/rectum epithelial cells use 32. Positive Low gives no negative-IHC credit.',570,22);
}
{
const s=slide('Published lists and current assignments','Sources',data.papers.map(r=>`${r.paper}: https://doi.org/${r.doi}`).join('\n')+'\nPaper subsets are unioned before counting; these are groups, not original paper rows.');
table(s,[['Paper','Mapped groups','CTA','Unsafe','Neither'],...data.papers.map(r=>[r.paper,n(r.mapped_groups),r.CTA,r.Unsafe,r.Neither])],{widths:[400,230,170,176,176],size:20,h:420,y:158});
note(s,`${data.papers.length} papers · ${n(data.paper_union)} unique mapped groups. List membership annotates selection; it does not grant admission.`,600,23);
}
{
const s=slide('Overlaps and biological association','Sources','Top six papers by mapped group count. Matrix counts shared selected groups; diagonal is selected count in each paper. Association counts refer to all reference groups and are independent of paper membership.');
const short=data.overlap_order.map(a=>a.replace('da Silva ','da Silva\n').replace(/([a-z])([0-9])/,'$1\n$2'));
table(s,[['Shared selected',...short],...data.overlap.map((r,i)=>[data.overlap_order[i].replace(/([a-z])([0-9])/,"$1 $2"),...r])],{widths:[240,152,152,152,152,152,152],size:21,h:305});
const a=data.association;
text(s,`Association: ${n(a.reproductive_only)} reproductive only  ·  ${n(a.developmental_only)} embryonic only  ·  ${n(a.both)} both`,64,508,1152,63,25,C.teal,true);
note(s,'Embryonic admission uses one gastrulation embryo. Later fetal atlases remain outside this release.',600,23);
}
{
const s=slide('Data used to filter normal expression','Evidence','Pinned sources: HPA v23, Monaco, DICE and Collin adult cornea. DICE public slots are not linked donor identities. HPA DICE summaries are not an independent dataset. Full URLs and integration states are in species-source-inventory.csv.');
table(s,[['Source','Measured quantity','Use'],['HPA consensus','Bulk nTPM','Organ tier, association and coverage'],['HPA tissue clusters','Cell-cluster nTPM; ≥20 reads','Cell compartments and rare-cell signals'],['HPA / Monaco immune','Native nTPM / TPM','Protected unless explicitly classified otherwise'],['Original DICE','TPM; third-highest reported value','Recurring immune RNA; coverage optional'],['Collin adult cornea','CP10K; support in ≥2 sample groups','Ocular signals; coverage optional'],['HPA + external IHC','Staining, confidence, tissue / cell','Positive findings and matched negative credit'],['Curated protein / clinical','Categorical findings + source evidence','Independent exclusions and retained flags']],{widths:[285,405,462],size:22,h:410});
note(s,'Missing measurements remain unknown. These adult panels do not establish pediatric safety.',600,24);
}
{
const s=slide('Every layer can contribute','Filtering funnel','Bars are disjoint categories after each cumulative layer. Earlier CTA counts mean eligible so far. Layer order changes first-removal attribution, not final classification. Full table on next slide.');
const labels=['Association','Tumor','HPA protein','Coverage','Bulk RNA','Tissue clusters','HPA immune','Monaco','DICE','Cornea','Curated protein','Clinical safety','External IHC'];
chart(s,labels,[{name:'CTA / eligible so far',values:data.funnel.map(r=>r.CTA_eligible_so_far),fill:C.teal},{name:'Unsafe',values:data.funnel.map(r=>r.Unsafe),fill:C.rust},{name:'Neither',values:data.funnel.map(r=>r.Neither),fill:'#BEA477'}],{h:440,y:157,stacked:true,size:21});
note(s,`${n(data.universe)} groups at every step. Final: ${n(data.summary.categories.CTA)} CTA / ${n(data.summary.categories.Unsafe)} Unsafe / ${n(data.summary.categories.Neither)} Neither.`,613,22);
}
{
const s=slide('Counts at each filtering choice','Filtering funnel','First-removal counts are sequential. Independent overlapping effects are retained separately in independent-layer-effects.csv. Unsafe means stronger evidence under policy, not demonstrated toxicity.');
const labels=['Association','Tumor ≥2','HPA protein','Required coverage','Normal bulk','HPA tissue clusters','HPA immune','Monaco immune','DICE immune','Adult cornea','Curated protein','Clinical safety','External IHC'];
table(s,[['Layer','New removals','CTA*','Unsafe','Neither'],...data.funnel.map((r,i)=>[labels[i],n(r.new_candidate_removals),n(r.CTA_eligible_so_far),n(r.Unsafe),n(r.Neither)])],{widths:[365,220,167,200,200],size:21,h:430,y:154});
note(s,'*Eligible so far until the final row. A later layer can change Neither to Unsafe.',609,22);
}
{
const s=slide('Effect of the policy choices','Sensitivity','Cumulative effects relative to the previous 343-group policy; admissions are not independent marginal effects. Current reprocessing also corrects grouped IHC levels and all-context assay evaluation.');
const names=['Previous policy','Double the linked RNA grid','Classify fibroblasts as Manageable','Allow specified large-bowel cells','Flag Low HPA protein','Flag explicitly Low imported protein'];
table(s,[['Cumulative rule','CTA','Added'],...data.components.map((r,i)=>[names[i],r.CTA,r.additional_CTA])],{widths:[760,196,196],size:26,h:350});
note(s,`${n(data.summary.review.corrected_group_ihc_levels)} grouped IHC labels corrected; ${data.summary.review.changed_layer_decisions} intermediate decisions changed; ${data.summary.review.changed_final_categories} final category changes versus the prior 463-group proposal.`,575,25);
}
{
const s=slide('Named targets in the current set','Assignments','NY-ESO-1 is CTAG1A/B; LAGE-1 is CTAG2. Flags do not establish toxicity and are not erased by admission. MAGE genes other than MAGEA4 were not calibration controls.');
const flag=r=>r.protein_group_label==='PRAME'?'Low protein; progenitor expression':r.protein_group_label==='TBXT'?'External review; adult disc':r.protein_group_label==='CT83'?'External-IHC review':r.category==='Unsafe'?'Bronchial basal RNA; other findings':'All gates pass';
table(s,[['Target / group','Category','Retained context'],...data.named.map(r=>[r.protein_group_label,r.category,flag(r)])],{widths:[330,180,642],size:23,h:420,y:158});
note(s,'No gene-name exemptions. Other MAGE genes enter only through the same rules.',607,23);
}
{
const s=slide('Current assignments of known CTAs','Current versus historical','Denominator: 213 historically annotated gene loci. Each is mapped to its current protein-group decision. This slide does not preserve old classifications. Every row and reason is in known-cta-decisions.csv.gz.');
chart(s,['CTA','Unsafe','Neither'],[{name:'Historical CTA loci, current decision',values:['CTA','Unsafe','Neither'].map(k=>data.historical_loci[k]||0),fill:C.teal}],{h:320,labels:true,y:190,size:30});
text(s,`${n(data.historical_loci.CTA)} retained  ·  ${n(data.historical_loci.Unsafe)} Unsafe  ·  ${n(data.historical_loci.Neither)} Neither`,64,551,1152,56,31,C.ink,true);
}
{
const s=slide('Expression distributions and named targets','Normal and reproductive RNA','Standalone scientific figure with editable source CSVs. Axes are HPA native nTPM; +0.05 displays zero on logarithmic axes. Color reflects the entire cascade, not a two-axis cutoff. TBXT can qualify through embryonic expression.');
figure(s,'named-antigen-distribution',64,142,1152,510);
}
{
const s=slide('PRAME, CTAG2 and TBXT','Evidence interpretation','PRAME normal protein: https://pubmed.ncbi.nlm.nih.gov/33890816/ ; IMA203: https://www.nature.com/articles/s41591-025-03650-6 ; fibroblasts: https://elifesciences.org/articles/69854 ; CEA colitis: https://pubmed.ncbi.nlm.nih.gov/21157437/ ; TBXT: https://www.oncotarget.com/article/3086/text/ ; adult disc: https://www.nature.com/articles/s41598-017-01567-w');
table(s,[['Target','Measured normal finding','Current treatment'],['PRAME','Adrenal 3.6 / ceiling 4; Low protein','CTA; protein observation retained'],['PRAME','Heart fibroblasts 12 / 16','CTA; stromal tolerance assumption'],['PRAME','Colon undifferentiated cells 20.3 / 32','CTA; progenitor caution retained'],['CTAG2','Bronchial basal respiratory cells 24.9 / 4','Unsafe; stronger boundary is 8'],['TBXT','External antibody disagreement; adult disc','CTA; tissue-specific review retained']],{widths:[155,533,464],size:23,h:340});
note(s,'RNA values are nTPM. Clinical targeting supports product-specific evidence; it does not validate these thresholds or prove gene-wide safety.',555,25);
}
{
const s=slide('The same evidence categories across species','Species coverage','Native source scopes and integration states are in the package inventory. BarkBase: https://pubmed.ncbi.nlm.nih.gov/31181663/ ; DoGA: https://www.nature.com/articles/s41467-024-52798-1 ; Pfizer: https://doi.org/10.1038/s41597-026-07512-9 ; Tabula Muris/Senis: https://www.nature.com/articles/s41586-020-2496-1 ; MOCA: https://www.nature.com/articles/s41586-019-0969-x ; GXD: https://www.informatics.jax.org/mgihome/GXD/aboutGXD.shtml');
table(s,[['Evidence role','Human: used / annotated','Dog: available / registered','Mouse: identified only'],['Nomination','Published lists','14 native nominations','No curated snapshot'],['Reproductive / developmental','HPA / Tyser','Bgee / BarkBase / DoGA*','MOCA / gastrulation / GXD'],['Tumor','132 cohorts + HPA','Study assays; cohort blocked','Native cohort needed'],['Normal bulk','HPA consensus','Bgee / BarkBase / Pfizer*','GXD assay resource'],['Normal cell / immune','HPA / Monaco / DICE / cornea','No equivalent integrated panel','Tabula Muris/Senis'],['Normal protein','HPA / external evidence','Pfizer iBAQ / native IHC','GXD assay resource']],{widths:[255,300,315,282],size:21,h:390});
note(s,'*Metadata only for DoGA and Pfizer RNA. No mouse atlas is integrated; human thresholds and safety do not transfer by orthology.',585,23);
}
{
const s=slide('A simple API with explicit coverage','Release API','See docs/cta-selection.md. Gene identities retain namespaces, reference identity, aliases, historic IDs and unresolved alternatives. Representative protein groups do not assert exhaustive isoform identity. Old exact-protein background products remain independently versioned.');
text(s,'cta_selection.results(species="human", category="CTA")\ncta_selection.coverage(species="dog")\ncta_selection.sources(species="mouse")',64,184,1152,150,31,C.teal,true);
table(s,[['Species','Result entities','Assessment'],['Human','19,863 groups','Current research-policy decisions'],['Dog','14 native nominations','Neither / unassessed; policy unavailable'],['Mouse','0 rows; stable schema','Curation unavailable, not biological absence']],{y:385,h:190,widths:[190,370,592],size:25});
note(s,'IDs, reasons and flags stay attached. Missing history is unknown; ambiguous names never silently collapse.',603,23);
}
{
const s=slide('What the release establishes','Scope','Literature and implementation review: '+URL+'literature-review/ . All selected rows and excluded rows have current assignments. The versioned exact-protein near-self background remains separate.');
text(s,'463 groups under one reproducible human policy',64,182,1152,57,36,C.teal,true);
text(s,'Complete decisions, named-target evidence and source overlaps\nConsistent species schema with honest missing-data states\nProtein observations, donor uncertainty and native units preserved',64,278,1152,160,29);
text(s,'Still unvalidated: numerical toxicity cutoffs, mucosal tolerance, pediatric use and cross-assay equivalence. One early embryo does not provide comprehensive fetal coverage.',64,491,1152,128,29,C.rust);
}
for(const [i,s] of p.slides.items.entries()){
 await fs.writeFile(path.join(build,'rendered',`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await(await p.export({slide:s,format:'png',scale:1.5})).arrayBuffer()));
}
await fs.writeFile(path.join(build,'manifest.json'),JSON.stringify({policy:data.summary.policy,commit:'report-sha256:'+data.report_sha256,slideSize:{width:1280,height:720},pages},null,2));
const candidate=path.join(build,'candidate.pptx');
await(await PresentationFile.exportPptx(p)).save(candidate);
const finalPath=path.join(audit,'slides',process.env.CTA_SLIDE_FILENAME||'cta-selection-1.21.0.pptx');
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath,pythonExecutable:process.env.RUNTIME_PYTHON,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tableOwners.flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:tableOwners,requiredNativeChartOwnerSlides:chartOwners,materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,receiptPath:path.join(build,'validation.json')});
console.log(JSON.stringify({slides:pages.length,finalPath,font,tableOwners,chartOwners}));
