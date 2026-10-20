from pathlib import Path
import csv, hashlib, json, shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

BASE=Path(r"data/external/local_de73c97db013")
ROOT=BASE/"KC_Fib_SPTDI_DPTDI_whole_slide_revision_20260825"
OUT=ROOT/"05_delivery_PDF_AI"
KC_PAN=BASE/"KC_Fib_Mode整片统计_PDF_AI_20260824"/"01_KC_vertical_all_slides"
FIB_PAN=BASE/"Fib_K14_output"/"03_roi_handoff"/"KRT14_PDGFRA_DAPI_only"/"05_vertical_review_panels"
PLOTS=ROOT/"03_lineplots"
for d in [OUT/"KC_panels_PDF",OUT/"Fib_panels_PDF",OUT/"statistics_PDF",OUT/"AI"]: d.mkdir(parents=True,exist_ok=True)

for p in sorted(KC_PAN.glob("*.pdf")): shutil.copy2(p,OUT/"KC_panels_PDF"/p.name)
for p in sorted(FIB_PAN.glob("*.pdf")): shutil.copy2(p,OUT/"Fib_panels_PDF"/p.name)
plot_names=["KC_existing_whole_slide_SPTDI_DPTDI_main","KC_SPTDI_DPTDI_whole_slide_sensitivity","Fib_K14_SPTDI_DPTDI_whole_slide"]
for name in plot_names: shutil.copy2(PLOTS/f"{name}.pdf",OUT/"statistics_PDF"/f"{name}.pdf")

plt.rcParams.update({"font.family":"Arial","font.size":8,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
rows={
 "KC":[("Mode 1",KC_PAN/"03_007_vertical.png"),("Mode 2",KC_PAN/"05_bs1_vertical.png"),("Mode 3",KC_PAN/"07_BW4_vertical.png"),(None,PLOTS/"KC_existing_whole_slide_SPTDI_DPTDI_main.png")],
 "Fib_K14":[("Mode 1",FIB_PAN/"007_DAPI_KRT14_PDGFRA_vertical_review.png"),("Mode 2",FIB_PAN/"BS1_DAPI_KRT14_PDGFRA_vertical_review.png"),("Mode 3",FIB_PAN/"BW7_DAPI_KRT14_PDGFRA_vertical_review.png"),(None,PLOTS/"Fib_K14_SPTDI_DPTDI_whole_slide.png")],
}
for family,items in rows.items():
    fig=plt.figure(figsize=(168/25.4,28/25.4),facecolor="white")
    for i,(title,p) in enumerate(items):
        ax=fig.add_axes([i*.25,0,.25,1]); ax.imshow(Image.open(p)); ax.axis("off")
        if title: ax.text(.02,.98,title,transform=ax.transAxes,color="white",fontsize=8,fontweight="bold",va="top",bbox={"facecolor":"black","edgecolor":"none","alpha":.6,"pad":.35})
    for ext in ["pdf","svg","png"]: fig.savefig(OUT/f"{family}_Mode1_Mode2_Mode3_grouped_whole_slide_168x28mm.{ext}",dpi=600 if ext=="png" else None)
    plt.close(fig)

pdfs=sorted(OUT.rglob("*.pdf")); jobs=[]
for p in pdfs:
    rel=p.relative_to(OUT).with_suffix("")
    name="__".join(rel.parts)+".ai"
    jobs.append([p.as_posix(),(OUT/"AI"/name).as_posix()])
def js(s): return json.dumps(s,ensure_ascii=True)
jsx=OUT/"AI"/"convert_all_revision_pdf_to_ai.jsx"
lines=["app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;","var jobs=["+",".join("["+",".join(js(x) for x in j)+"]" for j in jobs)+"];",f"var out=new File({js((OUT/'AI'/'ai_validation.csv').as_posix())});",'out.encoding="UTF-8";out.open("w");out.writeln("PDF,AI,Opened,Width_pt,Height_pt,PDF_Compatible,Error");','function q(s){return \'"\'+String(s).replace(/"/g,\'""\')+\'"\';}','for(var i=0;i<jobs.length;i++){var d=null;try{d=app.open(new File(jobs[i][0]));var o=new IllustratorSaveOptions();o.pdfCompatible=true;o.compressed=true;o.embedLinkedFiles=true;d.saveAs(new File(jobs[i][1]),o);out.writeln([q(jobs[i][0]),q(jobs[i][1]),"true",d.width,d.height,"true",q("")].join(","));}catch(e){out.writeln([q(jobs[i][0]),q(jobs[i][1]),"false","","","false",q(e.message)].join(","));}if(d){d.close(SaveOptions.DONOTSAVECHANGES);}}','out.close();"DONE";']
jsx.write_text("\n".join(lines),encoding="utf-8")
print(jsx); print(f"PDF_JOBS={len(jobs)}")
