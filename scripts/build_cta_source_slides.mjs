// Editable release deck from presentation-data.json; scientific scatter from the audit.
// Required: RUNTIME_NODE_MODULES, ONCOREF_PRESENTATIONS_SKILL, RUNTIME_PYTHON.
// Run prepare_cta_release_presentation.py first. Output goes to the versioned audit.
import fs from 'node:fs/promises';
import {readFileSync} from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
import {fileURLToPath, pathToFileURL} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const audit=path.join(root,'docs/audits/cta-source-walkthrough-20261028');
const build=path.join(root,'tmp/cta-source-slides');
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
 text(s,'OncoRef 1.21.0 evidence walkthrough   •   28 October 2026',64,670,1100,26,16,C.muted);
 text(s,String(pages.length).padStart(2,'0'),1174,668,42,30,19,C.muted,true);
 s.speakerNotes.text=`${title}\n${notes}\nReport: ${URL}\nPolicy: ${data.summary.policy}\nFrozen decision report SHA256: ${data.report_sha256}`;
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

{
const s=slide('Human and dog CTA evidence','Sources and decisions','Current human selection is unchanged. Dog nomination counts are not a genome-wide selection.');
text(s,'463 human CTA groups',64,202,1152,100,66,C.teal,true);
text(s,'481 loci   /   17,388 Unsafe groups   /   2,012 Neither groups',64,335,1152,65,32);
text(s,'Dog: 14 native nominations, quantitative selection unassessed',64,446,1152,95,32);
note(s,'Source contributions, overlaps, measured contexts and explicitly labeled policy comparisons.',580,26);
}
{
const s=slide('Published candidate lists','Human nominations',data.papers.map(r=>`${r.paper}: https://doi.org/${r.doi}`).join('\n'));
table(s,[['Paper','Groups','CTA','Unsafe','Neither'],...data.papers.map(r=>[r.paper,n(r.mapped_groups),r.CTA,r.Unsafe,r.Neither])],{widths:[440,180,170,180,182],size:20,h:420,y:158});
note(s,'Union: 2,482 groups / 2,538 loci. 418 selected groups have an annotation; 45 do not.',600,23);
}
{
const s=slide('How published lists overlap','Human nominations','Paper tables are unioned. Counts represent mapped protein-equivalence groups. Full source IDs and pairwise matrix: human-overlaps.csv.gz in the walkthrough.');
figure(s,'published-list-overlaps',36,152,1208,500);
}
{
const s=slide('Why 12,139 groups enter the funnel','Positive expression contexts','Tyser et al. (2021): https://doi.org/10.1038/s41586-021-04158-y . The complete reference has 19,863 groups. HPA supplies reproductive nTPM. Developmental presence alone is not restriction.');
table(s,[['Qualifying route','Entrants','Final CTA'],['Reproductive only','1,231','425'],['Developmental only','9,481','29'],['Both','1,427','9'],['Total','12,139','463']],{widths:[620,266,266],h:300,size:29});
text(s,'Whole-reference RNA measurements create this set',64,508,1152,50,32,C.teal,true);
note(s,'One early embryo supplies the developmental route. Later fetal atlases remain outside the current policy.',578,25);
}
{
const s=slide('Normal and cancer expression sources','Human evidence','URLs and versions: oncoref/data/cta-selection-sources.csv. HPA v23 https://v23.proteinatlas.org/about/download ; DICE https://dice-database.org/ ; cornea https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683 ; tumor registry and sample lineage in the frozen report.');
table(s,[['Evidence','Source','Quantity'],['Normal / reproductive bulk','HPA consensus with GTEx','nTPM'],['Normal tissue-cell clusters','HPA tissue clusters','nTPM with read support'],['Normal immune cells','HPA / Monaco / original DICE','Native nTPM / TPM'],['Adult cornea','Collin et al., GSE155683','CP10K'],['Normal protein','HPA + curated external IHC','Staining and reagent confidence'],['Cancer','132 compatible OncoRef cohort views + HPA','TPM / protein-coding TPM'],['Developmental','Tyser et al. (2021), one embryo','Lineage mean TPM']],{widths:[310,570,272],size:22,h:405});
note(s,'14,452 namespace-deduplicated tumor sample IDs. Patient independence and source overlap remain explicit.',600,23);
}
{
const s=slide('Active RNA cutoffs','Released policy','The second RNA boundary is a policy multiplier, not an antibody grade or a measured toxicity limit. Cornea CP10K uses nominal limits /100, not a TPM conversion.');
table(s,[['Candidate ceiling','Unmatched IHC','Approved −','Supported −','Enhanced −'],['Vital',4,8,16,32],['Protected',8,16,32,64],['Manageable',16,32,64,128],['Specified bowel cells',32,64,128,256]],{widths:[340,235,192,192,193],size:26,h:310});
text(s,'Above ceiling: hold. At least 2× ceiling: Unsafe evidence.',64,510,1152,52,30,C.teal,true);
note(s,'Negative-IHC credit requires verified tissue/cell coverage. Same-tissue staining alone cannot establish that match.',578,25);
}
{
const s=slide('Rules applied before final assignment','Released policy','All sources retain native units and missingness. RNA candidate equality passes. Every layer can contribute after a prior hold.');
table(s,[['Gate','Current rule'],['Positive context','Reproductive ≥5 nTPM and > adult maximum, OR supported embryonic mean ≥5 TPM'],['Tumor','≥2 TPM / pTPM in at least one compatible sample'],['Coverage','Required bulk / cell / immune panels. DICE and cornea optional'],['Normal RNA','Tissue- and IHC-dependent ceilings'],['Normal protein','Low retains a flag. Accepted stronger positive evidence excludes'],['Final decision','Unsafe evidence takes precedence. Remaining holds give Neither']],{widths:[260,892],size:24,h:390});
note(s,'These are active research-policy cutoffs. Adult normal panels do not establish pediatric safety.',600,25);
}
{
const s=slide('Every filtering layer can contribute','Human funnel','Eligible so far becomes a final CTA category only at the last stage. Sequential removals depend on display order; final decisions do not.');
const labels=['Association','Tumor','HPA protein','Coverage','Bulk RNA','Tissue clusters','HPA immune','Monaco','DICE','Cornea','Curated protein','Clinical safety','External IHC'];
chart(s,labels,[{name:'Eligible so far',values:data.funnel.map(r=>r.CTA_eligible_so_far),fill:C.teal},{name:'Unsafe',values:data.funnel.map(r=>r.Unsafe),fill:C.rust},{name:'Neither',values:data.funnel.map(r=>r.Neither),fill:'#BEA477'}],{h:435,y:158,stacked:true,size:21});
note(s,'19,863 groups at every stage. Final: 463 CTA / 17,388 Unsafe / 2,012 Neither.',610,24);
}
{
const s=slide('Counts at each filtering choice','Human funnel','These are sequential first-removal counts, not independent marginal effects. Full source decisions remain available after an earlier failure.');
const labels=['Association','Tumor ≥2','HPA protein','Required coverage','Normal bulk','HPA tissue clusters','HPA immune','Monaco immune','DICE immune','Adult cornea','Curated protein','Clinical safety','External IHC'];
table(s,[['Stage','Removed','Eligible*','Unsafe','Neither'],...data.funnel.map((r,i)=>[labels[i],n(r.new_candidate_removals),n(r.CTA_eligible_so_far),n(r.Unsafe),n(r.Neither)])],{widths:[365,220,167,200,200],size:21,h:430,y:154});
note(s,'*Eligible so far until the final row. Later evidence can change Neither to Unsafe.',610,23);
}
{
const s=slide('Widely studied antigens','Current assignments','Comparison genes are not forced positive controls. NY-ESO-1 is CTAG1A/B; LAGE-1 is CTAG2. Evidence includes retained flags.');
const reasons={PRAME:'Adrenal 3.6, heart fibroblasts 12, colon cells 20.3 nTPM',MAGEA1:'Reproductive 4.4 nTPM is below association floor 5',MAGEA3:'Colon undifferentiated cells 53.1 nTPM, ceiling 32',MAGEA12:'Reviewed clinical injury with attribution uncertainty',CTAG2:'Bronchial basal cells 24.9 nTPM, Unsafe boundary 8',TBXT:'Embryonic 245.1 TPM. Adult-disc / IHC review retained',CT83:'Quantitative gates pass. External-IHC review retained'};
table(s,[['Antigen / group','Decision','Reason / retained context'],...data.named.map(r=>[r.protein_group_label.replace('MAGEA','MAGE-A'),r.category,reasons[r.protein_group_label]||'Passes all current gates'])],{widths:[240,160,752],size:20,h:430,y:154});
note(s,'Known CTA loci receive current decisions: 162 CTA / 16 Unsafe / 35 Neither.',610,23);
}
{
const s=slide('Reproductive and developmental expression','Widely studied antigens','Axes include zero using a symmetric-log scale. Missing measurements are omitted. Classification uses the full funnel. Embryonic and reproductive measurements retain separate units.');
figure(s,'widely-studied-antigens',36,148,1208,504);
}
{
const s=slide('PRAME has colon IHC','Cell-type coverage','HPA v23 PRAME: https://v23.proteinatlas.org/ENSG00000185686-PRAME/tissue/colon . The frozen table reports Enhanced Not detected glandular, endothelial and peripheral nerve/ganglion cells. Exact undifferentiated-cluster coverage remains unresolved.');
table(s,[['Observation','Measured evidence'],['Colon RNA','Undifferentiated cells: 20.3 nTPM, 96 reads'],['Colon IHC','Enhanced Not detected: glandular, endothelial, nerve/ganglion'],['Current ceiling','32 without matched negative-IHC credit'],['Bowel moved to Manageable 16','PRAME would become Neither'],['Matched negative IHC, if verified','Approved would give 32; Enhanced would give 128']],{widths:[410,742],size:25,h:350});
note(s,'Same-tissue negatives are retained. Assigning cell-type coverage requires a supported anatomical match.',575,27);
}
{
const s=slide('Dog nominations and funnel coverage','Dog evidence','Nemec et al. (2019): https://doi.org/10.1111/vco.12475 ; Hindriks et al. (2025): https://doi.org/10.1016/j.vetimm.2025.110946 . Nomination rows can represent ambiguous assay families.');
table(s,[['Stage','Evidence / gap'],['Published nominations','Nemec 5 + Hindriks 10 − shared NT5C1B 1 = 14'],['Identity','13 unique nominations, one ambiguous MAGE family assay'],['Normal / developmental RNA','Native matrices available. Reference / scale reconciliation needed'],['Tumor RNA','GSE304066 count/FPKM and annotation conflicts unresolved'],['Normal protein','Limited study IHC and pooled Pfizer iBAQ'],['Quantitative selection','Unassessed. All 14 nominations remain Neither']],{widths:[360,792],size:24,h:390});
note(s,'No genome-wide canine CTA or Unsafe count has been established.',600,26);
}
{
const s=slide('Dog data sources and tracked integration','Dog evidence','Source issues https://github.com/pirl-unc/oncoref/issues/574 and https://github.com/pirl-unc/oncoref/issues/585 . Pfizer https://doi.org/10.1038/s41597-026-07512-9 ; DoGA https://doi.org/10.1038/s41467-024-52798-1 ; BarkBase https://pubmed.ncbi.nlm.nih.gov/31181663/ .');
table(s,[['Role','Source and status'],['Normal / developmental RNA','Bgee native TPM; BarkBase 37,149 features × 150 TMM-CPM columns'],['Additional RNA','DoGA and Pfizer RNA: metadata registered, full matrices not integrated'],['Tumor cohorts','GSE304066 provenance blocked; more cohorts tracked in #574'],['Normal MS protein','Pfizer 9,929 reported groups × 38 tissue pools; integration #585'],['Normal IHC / cell atlases','Limited study IHC. Systematic IHC and cell coverage remain gaps']],{widths:[330,822],size:24,h:345});
note(s,'46 Bgee/BarkBase libraries overlap. Pooled MS cannot identify individual donors or cell types.',575,26);
}
{
const s=slide('Neither can have a simpler definition','Proposed policy, not released','This removes the arbitrary intermediate RNA band. It requires replaying RNA and protein with the same anatomy policy, not only relabeling plots. Missingness and uncertainty remain separate reason/status fields.');
table(s,[['Category','Proposed meaning'],['Unsafe','Credible expression above the single vital/protected cutoff, or accepted clinical harm evidence'],['CTA / oncofetal','Tumor + permitted-context association + adequate coverage + acceptable normal expression'],['Neither','No tumor or association, low expression, other-compartment expression, or unresolved evidence']],{widths:[250,902],size:28,h:310});
note(s,'A missing measurement differs from a low measurement. Every layer still evaluates every gene.',555,28);
}
{
const s=slide('Curated oncofetal nominations can mirror CTA curation','Developmental sources','Shaw et al. (2024): https://doi.org/10.1038/s41467-024-47649-y ; TRIM71 study (2024): https://pmc.ncbi.nlm.nih.gov/articles/PMC11388079/ ; Cardoso-Moreira et al. (2019): https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/ . These newly identified sources are not in current selection counts.');
table(s,[['Source','Contribution and scope'],['Shaw et al. (2024)','157 cancer-specific targets, six oncofetal annotations. Retain exon-specific FN1/TNC scope'],['TRIM71 liver study (2024)','Top-200 oncofetal expression signature. Liver-specific comparison'],['Cardoso-Moreira et al. (2019)','Seven-organ developmental trajectories with child/adult measurements'],['Current Tyser route','Measured embryonic expression. Does not itself establish restriction']],{widths:[390,762],size:25,h:345});
note(s,'Published nominations still need independent normal/tumor filtering and explicit postnatal coverage.',580,26);
}
{
const s=slide('A shared audit API across species','API','The same schema retains native units, references and evidence gaps. Existing cta_selection.results and cta_selection.sources remain unchanged. Audit sources count published-list contributions; the existing sources function inventories measurement resources.');
text(s,'audit = cta_selection.audit(species="human")\naudit.sources     audit.overlaps     audit.funnel\naudit.criteria    audit.explain("CTAG2")',64,184,1152,160,31,C.teal,true);
table(s,[['Species','Assessment'],['Human','Released decisions with detailed measured contexts'],['Dog','Native nominations and explicit unassessed stages'],['Mouse','Stable schema with unavailable curation']],{y:390,h:192,widths:[230,922],size:27});
note(s,'Offline, checksum-verified tables. Human RNA cutoffs do not transfer by orthology.',610,24);
}
for(const [i,s] of p.slides.items.entries()){
 await fs.writeFile(path.join(build,'rendered',`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await(await p.export({slide:s,format:'png',scale:1.5})).arrayBuffer()));
}
await fs.writeFile(path.join(build,'manifest.json'),JSON.stringify({policy:data.summary.policy,commit:'report-sha256:'+data.report_sha256,slideSize:{width:1280,height:720},pages},null,2));
const candidate=path.join(build,'candidate.pptx');
await(await PresentationFile.exportPptx(p)).save(candidate);
const finalPath=path.join(audit,'slides',process.env.CTA_SLIDE_FILENAME||'cta-evidence-walkthrough.pptx');
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath,pythonExecutable:process.env.RUNTIME_PYTHON,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tableOwners.flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:tableOwners,requiredNativeChartOwnerSlides:chartOwners,materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,receiptPath:path.join(build,'validation.json')});
console.log(JSON.stringify({slides:pages.length,finalPath,font,tableOwners,chartOwners}));
