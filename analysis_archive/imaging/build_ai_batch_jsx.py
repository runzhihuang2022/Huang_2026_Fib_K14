from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(r"data/external/local_1f040df3bed6")
AI = ROOT / "07_AI"
AI.mkdir(parents=True, exist_ok=True)

pdfs = sorted(p for p in ROOT.rglob("*.pdf") if AI not in p.parents)
jobs = []
for p in pdfs:
    rel = p.relative_to(ROOT)
    safe = "__".join(rel.with_suffix("").parts)
    jobs.append([str(rel).replace("\\", "/"), p.as_posix(), (AI / f"{safe}.ai").as_posix()])

def js(s: str) -> str:
    return json.dumps(s, ensure_ascii=True)

lines = [
    "app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;",
    "var jobs=[" + ",".join("[" + ",".join(js(x) for x in row) + "]" for row in jobs) + "];",
    f"var out=new File({js((AI / 'ai_validation.csv').as_posix())});",
    'out.encoding="UTF-8"; out.open("w");',
    'out.writeln("Source,AI_Path,Opened,Illustrator_Version,Width_pt,Height_pt,Layers,TextFrames,PathItems,GroupItems,PlacedItems,RasterItems,PDF_Compatible,Error");',
    'function q(s){return \'"\'+String(s).replace(/"/g,\'""\')+\'"\';}',
    'for(var i=0;i<jobs.length;i++){',
    ' var doc=null;',
    ' try{',
    '  doc=app.open(new File(jobs[i][1]));',
    '  if(doc.layers.length>0){doc.layers[0].name="Artwork";}',
    '  var opt=new IllustratorSaveOptions(); opt.pdfCompatible=true; opt.compressed=true; opt.embedLinkedFiles=true; opt.embedICCProfile=false; opt.fontSubsetThreshold=100.0;',
    '  doc.saveAs(new File(jobs[i][2]),opt);',
    '  out.writeln([q(jobs[i][0]),q(jobs[i][2]),"true",q(app.version),doc.width,doc.height,doc.layers.length,doc.textFrames.length,doc.pathItems.length,doc.groupItems.length,doc.placedItems.length,doc.rasterItems.length,"true",q("")].join(","));',
    ' }catch(e){out.writeln([q(jobs[i][0]),q(jobs[i][2]),"false",q(app.version),"","","","","","","","","false",q(e.message)].join(","));}',
    ' if(doc){doc.close(SaveOptions.DONOTSAVECHANGES);}',
    '}',
    'out.close(); "DONE";',
]

jsx = AI / "convert_all_pdfs_to_ai.jsx"
jsx.write_text("\n".join(lines), encoding="utf-8")
(AI / "pdf_to_ai_jobs.json").write_text(json.dumps(jobs, ensure_ascii=False, indent=2), encoding="utf-8")
print(jsx)
print(f"JOBS={len(jobs)}")
