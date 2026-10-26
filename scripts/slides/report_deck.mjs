/** Reusable native slide helpers for tables, charts, diagrams and cited reports. */
export function createReportDeck(Presentation, metadata) {
const C={ink:'#173B55',muted:'#536675',teal:'#087E8B',amber:'#A76708',rust:'#AC493B',pale:'#F3F6F8',line:'#D9E2E8',white:'#FFFFFF'};
const p=Presentation.create({slideSize:{width:1280,height:720}});p.theme.defaultFont='Arial';
const pages=[];let current;
function text(s,str,x,y,w,h,size=27,color=C.ink,bold=false){
  const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  sh.text=String(str);sh.text.style={fontSize:size,color,bold};sh.text.insets=0;
  current.text.push(String(str));return sh;
}
function line(s,x,y,w,color=C.line){return s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:1},fill:color,line:{fill:'none',width:0}});}
function slide(title,section,source='[1,2] Current target policy and release review'){
  const s=p.slides.add();s.background.fill=C.white;
  current={title,section,text:[],links:[]};pages.push(current);
  text(s,section.toUpperCase(),64,26,1152,25,15,C.teal,true);
  text(s,title,64,62,1152,95,43,C.ink,true);
  line(s,64,654,1152);
  text(s,source,64,669,1060,30,14,C.muted);
  text(s,String(pages.length).padStart(2,'0'),1172,666,44,32,18,C.muted,true);
  s.speakerNotes.text=`${title}\nPolicy: ${metadata.policy}\nRepository: ${metadata.commit}\nSources: ${source}`;
  return s;
}
function note(s,str,y=592,color=C.muted,size=24){text(s,str,64,y,1152,55,size,color);}
function table(s,values,{y=170,h=360,widths,size=25,rowHeights}={}){
  const t=s.tables.add({rows:values.length,columns:values[0].length,left:64,top:y,width:1152,height:h,columnWidths:widths,values});
  t.styleOptions={headerRow:true,bandedRows:false};t.borders.assign({style:'solid',fill:C.line,width:.5});
  values.forEach((row,r)=>{
    t.rows[r].height=rowHeights?.[r]??h/values.length;
    row.forEach((value,c)=>{const cell=t.getCell(r,c);cell.fill=r===0?C.ink:(r%2===0?C.pale:C.white);cell.text.style={fontSize:r===0?size-1:size,color:r===0?C.white:C.ink,bold:r===0||c===0};});
    current.text.push(row.join(' | '));
  });return t;
}
function box(s,str,x,y,w,h,color=C.ink,size=27){
  const b=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:color,line:{fill:'none',width:0}});
  b.text=str;b.text.style={fontSize:size,color:C.white,bold:true};b.text.insets=15;
  current.text.push(str);return b;
}
function connect(s,a,b,side='bottom',target='top'){
  return s.shapes.connect(a,b,{kind:'straight',fromSide:side,toSide:target,line:{fill:C.muted,width:2},tail:{type:'arrow',width:'med',length:'med'}});
}
function bullet(s,heading,body,y,color=C.ink){
  text(s,heading,64,y,320,80,29,color,true);
  text(s,body,410,y,806,95,27,C.ink);
}
function bars(s,categories,values,color,{y=174,h=360}={}){
  const ch=s.charts.add('bar',{position:{left:64,top:y,width:1130,height:h},categories:[...categories].reverse(),
    series:[{name:'Genes',values:[...values].reverse(),fill:color}],barOptions:{direction:'bar',grouping:'clustered',gapWidth:48},
    hasLegend:false,xAxis:{visible:true,majorGridlines:null,textStyle:{fontSize:26,fill:C.ink},line:{fill:'none',width:0}},yAxis:{visible:false,majorGridlines:null},
    dataLabels:{showValue:true,position:'outEnd',numberFormatCode:'#,##0',textStyle:{fontSize:27,bold:true,fill:C.ink}},chartFill:C.white,plotAreaFill:C.white});
  for (const style of [ch.xAxis.textStyle, ch.yAxis.textStyle, ch.dataLabels.textStyle]) style.typeface = 'Arial';
  current.text.push(...categories.map((c,i)=>`${c}: ${values[i]}`));return ch;
}


return {p,pages,C,text,slide,note,table,box,connect,bullet,bars};
}
