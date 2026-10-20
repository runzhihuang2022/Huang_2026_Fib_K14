from pathlib import Path
import json, shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT=Path(r"data/external/local_de73c97db013")
PAN=ROOT/"KC_Fib_Mode整片统计_PDF_AI_20260824"/"01_KC_vertical_all_slides"
STAT=ROOT/"KC_Spinous_Mig_whole_slide_statistics_native_pixel_v3"
OUT=STAT/"05_delivery_PDF_AI"; PDF=OUT/"PDF"; AI=OUT/"AI"
for d in [OUT,PDF,AI]: d.mkdir(parents=True,exist_ok=True)
items=[("Mode 1",PAN/"03_007_vertical.png"),("Mode 2",PAN/"05_bs1_vertical.png"),("Mode 3",PAN/"07_BW4_vertical.png"),(None,STAT/"03_figures"/"KC_whole_slide_lineplot.png")]
plt.rcParams.update({"font.family":"Arial","font.size":8,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(168/25.4,28/25.4),facecolor="white")
for i,(title,path) in enumerate(items):
    ax=fig.add_axes([i*.25,0,.25,1]); ax.imshow(Image.open(path)); ax.axis("off")
    if title: ax.text(.02,.98,title,transform=ax.transAxes,color="white",fontsize=8,fontweight="bold",va="top",bbox={"facecolor":"black","edgecolor":"none","alpha":.6,"pad":.35})
for ext in ["pdf","svg","png"]: fig.savefig(OUT/f"KC_Mode1_Mode2_Mode3_whole_slide_statistics_168x28mm.{ext}",dpi=600 if ext=="png" else None)
plt.close(fig)
pdfs=[]
for p in sorted(PAN.glob("*.pdf")):
    q=PDF/p.name; shutil.copy2(p,q); pdfs.append(q)
q=PDF/"KC_whole_slide_lineplot.pdf"; shutil.copy2(STAT/"03_figures"/"KC_whole_slide_lineplot.pdf",q); pdfs.append(q)
pdfs.append(OUT/"KC_Mode1_Mode2_Mode3_whole_slide_statistics_168x28mm.pdf")
def js(s): return json.dumps(s,ensure_ascii=True)
jobs=[[p.as_posix(),(AI/(p.stem+".ai")).as_posix()] for p in pdfs]
jsx=AI/"convert_all_kc_pdf_to_ai.jsx"
lines=["app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;","var jobs=["+",".join("["+",".join(js(x) for x in j)+"]" for j in jobs)+"];",f"var out=new File({js((AI/'ai_validation.csv').as_posix())});",'out.encoding="UTF-8";out.open("w");out.writeln("PDF,AI,Opened,Width_pt,Height_pt,PDF_Compatible,Error");','function q(s){return \'"\'+String(s).replace(/"/g,\'""\')+\'"\';}','for(var i=0;i<jobs.length;i++){var d=null;try{d=app.open(new File(jobs[i][0]));var o=new IllustratorSaveOptions();o.pdfCompatible=true;o.compressed=true;o.embedLinkedFiles=true;d.saveAs(new File(jobs[i][1]),o);out.writeln([q(jobs[i][0]),q(jobs[i][1]),"true",d.width,d.height,"true",q("")].join(","));}catch(e){out.writeln([q(jobs[i][0]),q(jobs[i][1]),"false","","","false",q(e.message)].join(","));}if(d){d.close(SaveOptions.DONOTSAVECHANGES);}}','out.close();"DONE";']
jsx.write_text("\n".join(lines),encoding="utf-8"); print(jsx); print(len(jobs))
