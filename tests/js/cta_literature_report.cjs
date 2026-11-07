const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
class Element{constructor(id){this.id=id;this.value='';this.innerHTML='';this.textContent='';this.children=[];this.listeners={};}add(x){this.children.push(x)}addEventListener(t,fn){this.listeners[t]=fn}replaceChildren(){this.children=[]}appendChild(x){this.children.push(x)}scrollIntoView(){}}
const nodes=new Map(),get=id=>{if(!nodes.has(id))nodes.set(id,new Element(id));return nodes.get(id)};
get('collection').value='core';get('coverage-filter').value='missing';
const data=JSON.parse(fs.readFileSync('docs/audits/cta-literature-20261107/review.json','utf8'));
const context=vm.createContext({document:{getElementById:get,createElement:()=>new Element()},Option:class{constructor(t,v){this.text=t;this.value=v}},fetch:async()=>({ok:true,json:async()=>data}),URLSearchParams,location:{search:''},console});
vm.runInContext(fs.readFileSync('scripts/cta_evidence_assets/literature/review.js','utf8'),context);
setImmediate(()=>{
 assert.equal(get('review-count').textContent,'96 of 96 reviewed proteins');assert.equal(get('source-count').textContent,'83 of 150 paper and supporting-source records');assert.equal(get('coverage-count').textContent,'33 of 602 released CTA records');
 const set=(id,value)=>{get(id).value=value;get(id).listeners.input()};
 set('baseline','no_match');assert.equal(get('review-count').textContent,'41 of 96 reviewed proteins');set('baseline','broad_table_only');assert.equal(get('review-count').textContent,'55 of 96 reviewed proteins');set('baseline','');
 set('outcome','normal_review');assert.equal(get('review-count').textContent,'18 of 96 reviewed proteins');set('query','SRY');assert.equal(get('review-count').textContent,'1 of 96 reviewed proteins');assert(get('review-table').innerHTML.includes('22568433'));assert(get('review-table').innerHTML.includes('ensembl.org'));set('query','');set('outcome','');
 set('collection','supplementary');assert.equal(get('source-count').textContent,'67 of 150 paper and supporting-source records');set('role','normal_tissue');assert(!get('source-table').innerHTML.includes('No matching sources'));
 set('coverage-filter','all');assert.equal(get('coverage-count').textContent,'602 of 602 released CTA records');set('coverage-filter','unlisted');assert.equal(get('coverage-count').textContent,'122 of 602 released CTA records');
 get('review-pages').children[1].onclick();assert.equal(get('review-pages').children[2].textContent,'Page 2 of 7');
 console.log('DOM rendering, baseline/outcome/source/coverage filters, gene/paper links and pagination passed.');
});
const explorerContext=vm.createContext({document:{getElementById:get},fetch:async()=>({ok:true,json:async()=>data})});
const explorerPrelude=fs.readFileSync('scripts/cta_evidence_assets/explorer.js','utf8').split('function renderComparison()')[0];
vm.runInContext(explorerPrelude,explorerContext);
setImmediate(()=>{
 const html=vm.runInContext(`findingsHTML([{url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC280673/',proteins:['SP7'],finding:'outdated adult-positive claim',category:'Normal review',title:'SP7'}])`,explorerContext);
 assert(html.includes('adult femoral osteoblast RNA was not detected'));
 assert(!html.includes('outdated adult-positive claim'));
 console.log('Updated source findings replace earlier explorer text without altering frozen measurements.');
});
