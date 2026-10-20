from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\BW7_reprocess_20260820"
)
SOURCE_AI = ROOT / "04_editable_three_panel_AI" / "02_converted_ai" / "Figure3G_BW7_vertical_3panel_editable.ai"
OUTPUT_DIR = ROOT / "04_editable_three_panel_AI" / "03_final_live_ai"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_AI = OUTPUT_DIR / "Figure3G_BW7_vertical_3panel_live_text.ai"
PREVIEW = OUTPUT_DIR / "Figure3G_BW7_vertical_3panel_live_text_preview.png"
JSX_PATH = OUTPUT_DIR / "build_live_text_layer.jsx"

A4_WIDTH_PT = 595.275590551
W = A4_WIDTH_PT / 4.0
SOURCE_PIXEL_WIDTH, SOURCE_PIXEL_HEIGHT = 7115, 2573
IMAGE_H = SOURCE_PIXEL_HEIGHT / (SOURCE_PIXEL_WIDTH / W)
HEADER_H, GAP_H = 24.0, 4.0
H = 3 * (HEADER_H + IMAGE_H) + 2 * GAP_H
SCALE_PT = (1000.0 / 1.0113780498504639) / (SOURCE_PIXEL_WIDTH / W)


def js_string(value: str | Path) -> str:
    return json.dumps(str(value).replace("\\", "/"), ensure_ascii=True)


def main() -> None:
    if not SOURCE_AI.exists():
        raise FileNotFoundError(SOURCE_AI)
    for target in (OUTPUT_AI, PREVIEW, JSX_PATH):
        if target.exists():
            target.unlink()
    panels = [
        ("KRT6B+SOX9+", [("KRT6B", "FF3038", 8), ("SOX9", "20C968", 8), ("DAPI", "5475FF", 8)]),
        ("PDGFRα+KRT14+", [("PDGFRα", "20C968", 8), ("KRT14", "FF3038", 8), ("DAPI", "5475FF", 8)]),
        ("Spatial correspondence", [("KRT6B+SOX9+", "FF3038", 7), ("PDGFRα+KRT14+", "20C968", 7), ("DAPI", "5475FF", 7)]),
    ]
    panel_js = json.dumps(panels, ensure_ascii=True)
    jsx = f'''app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;
var doc = app.open(new File({js_string(SOURCE_AI)}));
var layer = doc.layers.add(); layer.name = "Editable_Labels_ScaleBars";
var W={W:.12f}, H={H:.12f}, imageH={IMAGE_H:.12f}, headerH={HEADER_H:.12f}, gapH={GAP_H:.12f}, scalePt={SCALE_PT:.12f};
var panels={panel_js};
function color(hex){{var c=new RGBColor();c.red=parseInt(hex.substring(0,2),16);c.green=parseInt(hex.substring(2,4),16);c.blue=parseInt(hex.substring(4,6),16);return c;}}
function rect(top,left,width,height,hex,opacity){{var p=layer.pathItems.rectangle(top,left,width,height);p.stroked=false;p.filled=true;p.fillColor=color(hex);p.opacity=opacity;p.zOrder(ZOrderMethod.SENDTOBACK);return p;}}
function text(s,x,y,size,bold,hex){{var t=layer.textFrames.add();t.contents=s;t.position=[x,y];var a=t.textRange.characterAttributes;a.textFont=app.textFonts.getByName(bold?"Arial-BoldMT":"ArialMT");a.size=size;a.fillColor=color(hex);t.zOrder(ZOrderMethod.BRINGTOFRONT);return t;}}
function line(x1,y1,x2,y2,width,hex){{var p=layer.pathItems.add();p.setEntirePath([[x1,y1],[x2,y2]]);p.closed=false;p.filled=false;p.stroked=true;p.strokeWidth=width;p.strokeColor=color(hex);p.zOrder(ZOrderMethod.BRINGTOFRONT);return p;}}
var topY=H;
for(var i=0;i<panels.length;i++){{
  rect(topY,0,W,headerH,"000000",100); text(panels[i][0],2.5,topY-10.8,12,true,"FFFFFF");
  var x=2.5; for(var j=0;j<panels[i][1].length;j++){{var s=panels[i][1][j];var t=text(s[0],x,topY-21.0,s[2],true,s[1]);x+=t.width+4.5;}}
  var imageBottom=topY-headerH-imageH, barX=4.0, barY=imageBottom+5.0;
  rect(barY+11.0,barX-1.2,Math.max(scalePt,27.0)+2.4,13.0,"000000",64); line(barX,barY,barX+scalePt,barY,0.75,"FFFFFF"); text("1000 \\u00b5m",barX,barY+9.0,7,false,"FFFFFF");
  topY=imageBottom-(i<panels.length-1?gapH:0);
}}
var opt=new IllustratorSaveOptions();opt.pdfCompatible=true;opt.compressed=true;opt.embedLinkedFiles=true;opt.fontSubsetThreshold=100.0;
doc.saveAs(new File({js_string(OUTPUT_AI)}),opt);
var png=new ExportOptionsPNG24();png.antiAliasing=true;png.artBoardClipping=true;png.transparency=false;png.horizontalScale=300;png.verticalScale=300;
doc.exportFile(new File({js_string(PREVIEW)}),ExportType.PNG24,png);doc.close(SaveOptions.DONOTSAVECHANGES);"DONE";'''
    JSX_PATH.write_text(jsx, encoding="ascii")
    print(JSX_PATH)


if __name__ == "__main__":
    main()
