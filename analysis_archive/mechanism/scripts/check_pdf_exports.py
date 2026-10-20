from pathlib import Path
import fitz,json
O=Path(__file__).resolve().parents[1];rows=[]
for p in sorted((O/'figures').glob('0*.pdf')):
 d=fitz.open(p);page=d[0];spans=[s for b in page.get_text('dict')['blocks']if 'lines'in b for l in b['lines']for s in l['spans']]
 sizes=[s['size']for s in spans];outside=[s['text']for s in spans if fitz.Rect(s['bbox']).x0<-.2 or fitz.Rect(s['bbox']).y0<-.2 or fitz.Rect(s['bbox']).x1>page.rect.width+.2 or fitz.Rect(s['bbox']).y1>page.rect.height+.2]
 r=dict(file=p.name,pages=len(d),width_mm=page.rect.width*25.4/72,height_mm=page.rect.height*25.4/72,min_font_pt=min(sizes),max_font_pt=max(sizes),fonts=sorted(set(s['font']for s in spans)),out_of_canvas_text=outside)
 assert len(d)==1 and min(sizes)>=7.99 and max(sizes)<=9.01 and not outside,r
 rows.append(r)
(O/'reports/PDF_export_QA.json').write_text(json.dumps(rows,indent=2),'utf8');print(json.dumps(rows,indent=2))
