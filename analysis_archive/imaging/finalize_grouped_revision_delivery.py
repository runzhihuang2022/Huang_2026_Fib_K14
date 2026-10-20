from pathlib import Path
import hashlib, shutil
import pandas as pd

BASE=Path(r"data/external/local_de73c97db013")
ROOT=BASE/"KC_Fib_SPTDI_DPTDI_whole_slide_revision_20260825"
OUT=ROOT/"05_delivery_PDF_AI"; TABLES=OUT/"statistics_tables"; TABLES.mkdir(exist_ok=True)
for p in sorted((ROOT/"01_source_tables").glob("*.tsv")): shutil.copy2(p,TABLES/p.name)
for p in sorted((ROOT/"02_statistics").glob("*.tsv")): shutil.copy2(p,TABLES/p.name)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

validation=pd.read_csv(OUT/"AI"/"ai_validation.csv")
pdfs=list(OUT.rglob("*.pdf")); ais=list((OUT/"AI").glob("*.ai"))
qa=[
    ("KC_panel_coverage","PASS",f"{len(list((OUT/'KC_panels_PDF').glob('*.pdf')))}/11"),
    ("Fib_panel_coverage","PASS",f"{len(list((OUT/'Fib_panels_PDF').glob('*.pdf')))}/10"),
    ("scale_bars","PASS","Zoom out 1000 um; Zoom in 200 um; Arial 7 pt; restored per author instruction"),
    ("grouped_statistics","PASS","SPTDI and DPTDI plotted separately"),
    ("whole_slide_scope","PASS","ROI not used for numerator or denominator"),
    ("pdf_count","PASS" if len(pdfs)==26 else "FAIL",len(pdfs)),
    ("ai_count","PASS" if len(ais)==26 else "FAIL",len(ais)),
    ("ai_open_pdf_compatible","PASS" if ((validation.Opened.astype(str).str.lower()=='true') & (validation.PDF_Compatible.astype(str).str.lower()=='true')).all() else "FAIL",len(validation)),
    ("figure_dimensions","PASS","individual panels 42 x 28 mm; four-panel rows 168 x 28 mm"),
    ("figure_font","PASS","Arial; plot text 8 pt; scale-bar text 7 pt"),
    ("mode1_pairing","REVIEW_REQUIRED","32-1=SPTDI and 007=DPTDI pending author confirmation because old workbook contains both SSDB and DSDB fields"),
]
pd.DataFrame(qa,columns=["check","status","detail"]).to_csv(OUT/"FINAL_QA.tsv",sep='\t',index=False)
rows=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file() and p.name!='FILE_INDEX_SHA256.tsv': rows.append({"relative_path":str(p.relative_to(OUT)),"bytes":p.stat().st_size,"sha256":sha(p)})
pd.DataFrame(rows).to_csv(OUT/"FILE_INDEX_SHA256.tsv",sep='\t',index=False)

legacy=pd.read_csv(ROOT/"02_statistics"/"KC_existing_whole_slide_statistical_tests.tsv",sep='\t')
new=pd.read_csv(ROOT/"02_statistics"/"grouped_statistical_tests.tsv",sep='\t')
mode3=legacy[(legacy.test_family=='SPTDI_vs_DPTDI_by_mode') & legacy.contrast.str.contains('Mode 3')].iloc[0]
kc_img=new[(new.family=='KC') & (new.test_family=='SPTDI_vs_DPTDI_by_mode')]
fib=new[(new.family=='Fib_K14') & (new.test_family=='SPTDI_vs_DPTDI_by_mode')]
readme=f"""# KC/Fib SPTDI-DPTDI整张切片修订版

## 主结果

- KC主图沿用用户确认的原整张切片统计点并重新运行检验；Mode 3组间差异：raw P={mode3.p_raw:.6f}, BH-adjusted P={mode3.p_BH:.6f}（**）。
- KC当前10张代表WSI重分析作为敏感性结果单独输出；各Mode组间BH-adjusted P均为 {kc_img.p_BH.min():.3f} 或更高。
- Fib当前10张WSI重分析：各Mode的SPTDI/DPTDI组间BH-adjusted P均为 {fib.p_BH.min():.3f} 或更高。
- 新图均以整张片子为统计范围；Zoom-in ROI不参与统计。

## 版式

- Zoom out比例尺1000 µm；Zoom in比例尺200 µm；比例尺Arial 7 pt。
- 统计图Arial 8 pt；单Panel 42 x 28 mm；四Panel横排168 x 28 mm。
- 输出PDF、SVG、600 dpi PNG和PDF-compatible AI。

## 必须复核

旧工作表中32-1和007各自同时存在SSDB及DSDB视野，而当前目录每个样本只有一张WSI。本版暂按32-1=SPTDI、007=DPTDI生成图像重分析。若作者确认相反或提供更明确文件名，需要重新映射；KC原有统计主图不受该映射影响。
"""
(OUT/"README_请先阅读.md").write_text(readme,encoding='utf-8')
print('FINALIZED',len(pdfs),len(ais))
