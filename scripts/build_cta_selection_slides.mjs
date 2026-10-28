// Create the presentation from the audited counts and standalone scientific plots.
// Run with RUNTIME_NODE_MODULES and ONCOREF_PRESENTATIONS_SKILL pointing to the supplied runtime.
import fs from 'node:fs/promises';
import {readFileSync} from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
import {fileURLToPath, pathToFileURL} from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const audit = path.join(root, 'docs/audits/cta-cell-thresholds-20261026');
const build = path.join(root, 'tmp/cta-selection-slides');
const snapshot = JSON.parse(await fs.readFile(path.join(audit,'summary.json'),'utf8'));
if(snapshot.current_selected!==215 || snapshot.retained_after_overlay!==154 || snapshot.paper_union!==2538 || snapshot.papers!==13){
  throw Error('Audit snapshot changed; refresh the slide narrative and counts before rebuilding');
}
const skill = process.env.ONCOREF_PRESENTATIONS_SKILL;
const runtime = process.env.RUNTIME_NODE_MODULES;
if (!skill || !runtime) throw Error('Supply the installed presentations skill and runtime paths');
const req = createRequire(path.join(runtime, '__runtime__.cjs'));
const {Presentation, PresentationFile} = await import(pathToFileURL(req.resolve('@oai/artifact-tool')).href);
const {finalizePresentation, resolvePresentationFont, applyPresentationChartFont} = await import(pathToFileURL(path.join(skill, 'container_tools/artifact_tool_utils.mjs')).href);
const font = resolvePresentationFont();
await fs.mkdir(build, {recursive:true});
const presentation = Presentation.create({slideSize:{width:1280,height:720}});
const C={ink:'#172B33',teal:'#13776F',gray:'#476168',red:'#B94430',amber:'#B17312',light:'#EAF1EE',white:'#FFFFFF'};
const auditURL='https://github.com/pirl-unc/oncoref/tree/cta-clinical-policy-1.21.0/docs/audits/cta-cell-thresholds-20261026';

function text(s,t,x,y,w,h,size=24,color=C.ink,bold=false){
  const a=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  a.text=t;a.text.style={typeface:font,fontSize:size,color,bold};return a;
}
function slide(title,notes){
  const s=presentation.slides.add();s.background.fill=C.white;
  text(s,title,54,38,1172,68,42,C.ink,true);
  s.speakerNotes.text=notes+'\nAudit and downloadable data: '+auditURL;
  return s;
}
function table(s,values,x,y,w,h,widths,size=23){
  const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,values,columnWidths:widths});
  t.borders.assign({style:'solid',fill:'#D7E2DE',width:.5});
  t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({margins:{top:4,bottom:4,left:8,right:8}});
  for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){
    const cell=t.getCell(r,c);cell.fill=r===0?C.teal:(r%2?C.white:'#F2F6F4');
    cell.text.style={typeface:font,fontSize:size,color:r===0?C.white:C.ink,bold:r===0};
  }
  return t;
}
function figure(s,name,x,y,w,h){s.images.add({blob:new Uint8Array(readFileSync(path.join(audit,'figures',name+'.png'))),contentType:'image/png',fit:'contain',position:{left:x,top:y,width:w,height:h},alt:name.replaceAll('-',' ')});}
function foot(s,t){text(s,t,54,654,1172,48,18,C.gray);}

// Cover is completed before the remaining slides are authored.
{
 const s=presentation.slides.add();s.background.fill='#F6F8F7';
 text(s,'OncoRef',72,64,1100,42,25,C.teal,true);
 text(s,'CTA and oncofetal\nselection',72,176,1100,190,64,C.ink,true);
 text(s,'Sources, filtering evidence and tissue-specific RNA thresholds',76,407,1100,88,30,C.gray);
 text(s,'215 genes in the current research selection\nExploratory cell-RNA calibration shown separately',76,556,1100,80,24,C.teal);
 s.speakerNotes.text='Snapshot 26 October 2026. Current policy is the 215-gene research selection on draft PR #642. No gene-wide clinical safety claim. https://github.com/pirl-unc/oncoref/pull/642';
 await fs.writeFile(path.join(build,'slide-01.png'),new Uint8Array(await (await presentation.export({slide:s,format:'png',scale:1})).arrayBuffer()));
}
{
 const s=slide('Two cutoffs, set by tissue and IHC evidence','Illustrative calibration, not validated limits. Exact boundary equality excludes. Numbers are independent hypotheses by assay, not transfers from bulk RNA. IHC factors 1,1.25,1.5,2 apply to both cutoffs only for exact adult tissue and cell labels, negative staining, a reliable positive profile control, and no unresolved antibody concern. No IHC relaxation for activated immune states or unmatched corneal subtypes. No reproductive maximum rescue. An RNA exclusion is reference-review evidence, not confirmed protein or peptide presentation. https://www.proteinatlas.org/about/disclaimer');
 text(s,'Illustrative baseline cutoffs with absent or unusable IHC',54,112,1160,36,25,C.teal,true);
 table(s,[['Source / tissue tier','Unit','CTA','Reference'],['HPA cluster / vital','nTPM',2,5],['HPA cluster / protected','nTPM',5,10],['HPA cluster / manageable','nTPM',10,20],['HPA immune / vital','nTPM',2,5],['DICE or Monaco / vital','TPM',2,5],['Adult cornea / vital','CP10K',1,2]],54,165,710,330,[350,120,95,145],22);
 text(s,'Below CTA cutoff',810,170,400,34,25,C.teal,true);
 text(s,'Passes this RNA gate only',810,211,400,48,23);
 text(s,'Between the cutoffs',810,292,400,34,25,C.amber,true);
 text(s,'Exclude from CTA selection\nRetain as uncertain self evidence',810,332,400,75,23);
 text(s,'At or above reference cutoff',810,428,400,60,25,C.red,true);
 text(s,'RNA-supported exclusion review',810,487,400,55,23);
 text(s,'Negative IHC: missing/uncertain ×1; Approved ×1.25; Supported ×1.5; Enhanced ×2',54,544,1160,40,22,C.ink,true);
 text(s,'Relax only for matching tissue and cell type with a positive assay control. Reliable positive normal IHC excludes independently.',54,591,1160,52,22);
 foot(s,'Support guards: clusters ≥20 reads; DICE ≥3 reported values; cornea ≥5 positive cells and ≥20 counts in ≥2 distinct groups.');
}
{
 const s=slide('Nomination sources and their overlaps','Counts are deduplicated mapped coding genes per paper. Nested tables within a paper are unioned, not summed. The 13 registered papers yield 2,538 genes from 12 with mapped eligible identities. Boel BAGE remains an unresolved historical identity and contributes zero mapped loci. Trial evidence contributes two extra nominations, KIF20A and CLDN6, giving 2,540. The heatmap gives pairwise intersections including the diagonal. Source DOI links and complete membership are in nomination-sources.csv and nomination-membership.csv.gz.');
 const names=['da Silva 2017','Wang 2016','Chang 2019','Gong 2021','Bruggeman 2018','Jamin 2021','Loriot 2025','Carter 2023','Seager 2024','Bradley 2020','Bai 2016','Garg 2007','Boel 1995'];
 const counts=[1066,994,986,733,692,481,215,103,15,10,1,1,0];
 const selected=[187,140,162,167,7,29,50,30,6,0,0,0,0];
 table(s,[['Paper','Mapped','Selected'],...names.map((n,i)=>[n,counts[i],selected[i]])],54,124,475,430,[255,110,110],20);
 figure(s,'source-overlaps',550,115,675,470);
 text(s,'2,538 paper nominees + 2 trial nominees = 2,540 unique genes',54,583,1160,37,27,C.teal,true);
 foot(s,'CLDN6 is explicitly oncofetal in the trial registry. Broad oncofetal signatures and developmental atlases are not yet ingested.');
}
{
 const s=slide('Evidence used for normal-expression filtering','Coverage counts refer to the 2,540 current nominations in the pinned audit, not the full source databases. Bulk expression comes from HPA v23 consensus. Curated and primary antibody data must be distinguished. DICE summary and original DICE measurements are the same study, never independent replication. The adult corneal source is GSE155683, Collin et al. Donor summaries currently cover CTAG2 heart only, four donors. Clinical evidence qualifies development selection but does not prove normal expression absent. https://v23.proteinatlas.org/about/download https://pmc.ncbi.nlm.nih.gov/articles/PMC6289654/ https://pmc.ncbi.nlm.nih.gov/articles/PMC8343164/');
 table(s,[['Evidence','Coverage / unit','Role'],['HPA bulk RNA','2,534 genes / nTPM','Ceiling, reproductive selectivity, breadth'],['HPA IHC + primary antibodies','Tissue / cell reliability','Protein and assay conflicts'],['HPA tissue / cell clusters','1,095 genes / nTPM','Rare-cell RNA, ≥20 reads'],['HPA and Monaco immune','1,095 each / native nTPM or TPM','Sorted immune-cell RNA'],['Original DICE distributions','2,537 genes / TPM, 15 groups','Recurrence and activated T cells'],['Adult corneal atlas','1,071 genes / CP10K','Ocular cell types and sample groups'],['Donors, literature and clinical','Source and uncertainty records','Variability, protein and clinical context']],54,128,1172,430,[330,375,467],23);
 text(s,'Missing coverage remains unknown. RNA, protein staining and clinical outcomes retain separate meanings.',54,598,1172,48,26,C.teal,true);
 foot(s,'HPA DICE summaries are descriptive only. The current audit is a selected source extract, not complete atlas coverage.');
}
{
 const s=slide('The current filtering funnel','Current disjoint final categories: 2,540 nominations = 1,585 Unsafe + 740 Neither + 215 CTA. Unsafe detail: Protein 845, RNA 603, High RNA 128, Normal 8, Safety 1. Neither detail: Antibody 199, Low 196, Conflict 179, Cellular 70, Selectivity 64, Unassessed 31, Variable 1. These are assigned final reasons, not independent filter counts. Historical counts describe policy versions, not a single nested funnel: 758 pilot, 740 tightened selectivity and CTAG2 donor hold, 589 admission floor, 436 new evidence with one readmission, 399 high-RNA cap, 213 cell/protein review, 215 clinical qualification. The 436 transition removed 154 and readmitted CT55.');
 const chart=s.charts.add('bar',{position:{left:54,top:130,width:695,height:410},categories:['Nominations','After normal exclusions','After unresolved holds'],series:[{name:'Genes remaining',values:[2540,955,215],fill:C.teal}],hasLegend:false,dataLabels:{showValue:true,position:'outEnd'},yAxis:{majorGridlines:{fill:'#DCE4E0',width:1}},xAxis:{}});
 applyPresentationChartFont(chart,{fontFamily:font});
 chart.xAxis.textStyle.fontSize=18; chart.yAxis.textStyle.fontSize=16; chart.dataLabels.textStyle.fontSize=22;
 text(s,'1,585 normal exclusions',800,145,420,40,28,C.red,true);
 text(s,'Protein 845\nRNA 603\nHigh RNA 128\nOther normal / safety 9',800,198,420,148,25);
 text(s,'740 unresolved or ineligible',800,382,420,70,28,C.amber,true);
 text(s,'Assay and literature conflicts,\nlow admission evidence, cellular\nRNA, selectivity or missing data',800,454,420,110,24);
 text(s,'Historical selections: 740, 589, 436, 399, 213, then 215',54,592,1172,44,27,C.teal,true);
 foot(s,'215 = 204 expression passes + 9 supported admissions + 2 clinically qualified genes, PRAME and MAGE-A4. Historical steps are policy revisions.');
}
{
 const s=slide('Where the named CTAs sit in bulk RNA','Source: HPA v23 bulk consensus, same source on both axes. Adult maximum is across vital, protected and manageable bulk compartments. Plotted missing values are omitted and plotted n is shown in the legend. Axes use log10(1+x) with raw nTPM tick labels. This display alone is not the selection rule. PRAME reproductive 102.2/adult 3.6; MAGEA4 21.5/0.6; CTAG1A 8.1/0.1; CTAG1B 7.4/0.2. Clinical qualification retains normal-expression caveats.');
 figure(s,'bulk-rna-scatter',54,110,865,520);
 table(s,[['Gene','Repro.','Adult'],['PRAME',102.2,3.6],['MAGEA4',21.5,.6],['CTAG1A',8.1,.1],['CTAG1B',7.4,.2]],930,174,295,235,[115,90,90],20);
 text(s,'NY-ESO-1 is represented\nby CTAG1A and CTAG1B.\n\nBulk averages can miss\ncell-specific signals.',936,455,285,145,23,C.gray);
 foot(s,'Both axes use bulk nTPM. Cell-cluster nTPM, immune TPM and corneal CP10K are evaluated separately.');
}
{
 const s=slide('The stricter cell gate changes the shortlist','Illustrative overlay only; current 215-gene membership is unchanged. 61 current targets cross the lower boundary: 33 middle-zone and 28 upper-layer exclusions. PRAME: heart fibroblasts 12 nTPM with cutoffs 2/5, colon undifferentiated cells 20.3 with cutoffs 5/10, vascular fibroblasts 2.8 with cutoffs 2/5. No exact matched negative IHC for these cells. MAGEA4 and CTAG1A/B remain below the supported cell-RNA cutoffs. The maximum normalized ratio is a screening display across separately evaluated assays; it is not pooled RNA expression or a validated risk score. All 215 have at least some supported normal-cell measurements, but coverage is incomplete.');
 figure(s,'cell-cutoff-distribution',54,116,820,470);
 table(s,[['Illustrative overlay','N'],['Current selection',215],['Below added exclusion gate',154],['Middle zone',33],['Stronger RNA / protein layer',28]],887,147,338,265,[260,78],22);
 text(s,'PRAME crosses both cutoffs.\nMAGE-A4 and NY-ESO-1\nremain below the added gate.',887,444,338,132,23,C.red,true);
 text(s,'Threshold sensitivity: 123 to 168 of 215 retained when the lower cutoffs vary from 0.5× to 1.5×.',54,589,1160,55,24,C.teal,true);
 foot(s,'Only 2,849 of 314,586 screened cell contexts have exact tissue/cell IHC matches. Stronger IHC quality cannot repair a missing cell-type match.');
}
{
 const s=slide('Oncofetal intake and brachyury follow-up','Current nomination registry has CLDN6 trial evidence explicitly labelled oncofetal but no general oncofetal signature intake. Su et al. 2024 identify a 200-gene oncofetal liver expression program, Table S1; these are nomination leads, not established restricted antigens. https://www.thno.org/v14p4948.htm DOI 10.7150/thno.99633. Developmental data sources: Tyser gastrulation https://www.nature.com/articles/s41586-021-04158-y ; Cao fetal organ atlas https://pubmed.ncbi.nlm.nih.gov/33184181/ ; organogenesis sources and stage metadata in docs/developmental-targets.md. TBXT is not among the 2,540 current nominees and is Neither in canonical results. Adult nucleus-pulposus evidence: https://pubmed.ncbi.nlm.nih.gov/28473691/ . Fetal corneal annotations exist in the registry but are not used as admission measurements. No direct nTPM threshold transfer to counts or CP10K.');
 table(s,[['Source / target','Current status','What it adds'],['CLDN6 clinical trial evidence','Included as a nomination','Explicit oncofetal target, product-specific outcomes'],['Su 2024 oncofetal liver signature','200 reported genes, not imported','A broader fetal / tumor expression program'],['Gastrulation and organogenesis atlases','Assessed as sources, not ingested','Embryonic stage and lineage expression'],['Cao / DESCARTES fetal atlas','Assessed as a source, not ingested','Fetal organs, cell types and coverage'],['TBXT / brachyury','Not selected or currently nominated','Requires fetal admission and adult disc review']],54,134,1172,360,[370,360,442],23);
 text(s,'Developmental expression can expand nomination and admission, once measured in a calibrated source.',54,526,1150,58,28,C.teal,true);
 text(s,'The adult normal-expression gate still applies. A high fetal signal cannot erase adult cell expression.',54,599,1150,50,25);
 foot(s,'Su 2024 DOI: 10.7150/thno.99633. Its 200 genes are expression-program leads; overlap with the current intake is not yet computed.');
}

for(const [i,s] of presentation.slides.items.entries()){
 const blob=await presentation.export({slide:s,format:'png',scale:1});
 await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
const candidatePath=path.join(build,'candidate.pptx');
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);
const finalPath=path.join(audit,'cta-selection.pptx');
const receipt=await finalizePresentation({workspaceDir:root,candidatePath,finalPath,
 pythonExecutable:process.env.RUNTIME_PYTHON,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...[2,3,4,6,7,8].flatMap(n=>['--require-native-table-slide',String(n)])],
 explicitTotalSlideCount:8,requiredNativeTableOwnerSlides:[2,3,4,6,7,8],requiredNativeChartOwnerSlides:[5],
 materializeLiteralChartWorkbooks:true,fontPolicy:{basis:'design',families:[font]},
 verifyArtifactToolImport:true,receiptPath:path.join(build,'validation.json')});
console.log(JSON.stringify({font,slides:presentation.slides.items.length,finalPath,receipt},null,2));
