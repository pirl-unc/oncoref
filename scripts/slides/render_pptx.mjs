// Render a finalized PPTX for visual review and its companion PDF.
// RUNTIME_NODE_MODULES must point to the supplied artifact-tool runtime.
import fs from 'node:fs/promises';
import path from 'node:path';
import {createRequire} from 'node:module';
import {pathToFileURL} from 'node:url';
const [input,output]=process.argv.slice(2);
if(!input||!output||!process.env.RUNTIME_NODE_MODULES) throw Error('Usage: render_pptx.mjs input.pptx render_directory');
const req=createRequire(path.join(process.env.RUNTIME_NODE_MODULES,'__runtime__.cjs'));
const {FileBlob,PresentationFile}=await import(pathToFileURL(req.resolve('@oai/artifact-tool')).href);
const p=await PresentationFile.importPptx(await FileBlob.load(input));
await fs.mkdir(output,{recursive:true});
for(const [i,s] of p.slides.items.entries()){
 const blob=await p.export({slide:s,format:'png',scale:1.5});
 await fs.writeFile(path.join(output,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
console.log(`Rendered ${p.slides.items.length} final slides`);
