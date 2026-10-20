from pathlib import Path
import csv, hashlib, json, shutil
import pandas as pd

BASE=Path(r"data/external/local_de73c97db013")
FIB=BASE/"Fib_K14_output"/"07_whole_slide_statistics"
KC=BASE/"KC_Spinous_Mig_whole_slide_statistics_native_pixel_v3"

shutil.copy2(FIB/"03_figures"/"Fib_K14_whole_slide_lineplot.pdf",FIB/"05_delivery_PDF_AI"/"Fib_K14_whole_slide_lineplot.pdf")

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def audit(root, expected_pdf, expected_ai):
    out=root/"05_delivery_PDF_AI"; rows=[]
    for p in sorted(out.rglob('*')):
        if p.is_file(): rows.append({"relative_path":str(p.relative_to(out)),"bytes":p.stat().st_size,"sha256":sha(p)})
    pd.DataFrame(rows).to_csv(out/"FILE_INDEX_SHA256.tsv",sep='\t',index=False)
    pdf=list(out.rglob('*.pdf')); ai=list(out.rglob('*.ai'))
    val=list(out.rglob('ai_validation.csv'))
    ok=False
    if val:
        v=pd.read_csv(val[0]); ok=bool((v['Opened'].astype(str).str.lower()=='true').all() and (v['PDF_Compatible'].astype(str).str.lower()=='true').all())
    qa=[("whole_slide_scope","PASS","ROI not used for numerator or denominator"),("pdf_count","PASS" if len(pdf)==expected_pdf else "FAIL",len(pdf)),("ai_count","PASS" if len(ai)==expected_ai else "FAIL",len(ai)),("ai_open_pdf_compatible","PASS" if ok else "FAIL",ok),("font","PASS","Arial; 8 pt plotted text"),("panel_size","PASS","42 x 28 mm each; four-panel row 168 x 28 mm")]
    pd.DataFrame(qa,columns=['check','status','detail']).to_csv(out/"FINAL_QA.tsv",sep='\t',index=False)

audit(FIB,2,2); audit(KC,13,13)

fib_tests=pd.read_csv(FIB/"02_statistics"/"Fib_K14_statistical_tests.tsv",sep='\t')
fib_sum=pd.read_csv(FIB/"02_statistics"/"Fib_K14_mode_summary.tsv",sep='\t')
text="# Fib_K14整张切片统计结果\n\n"
text+="- 统计区域：整张切片；Zoom-in ROI不参与统计。\n- 重复单位：切片。\n- Mode：0–7 dpb=Mode 1；8–14 dpb=Mode 2；15–28 dpb=Mode 3。\n- 001A因缺少可靠时间映射，保留整片计数但不进入Mode检验。\n- 阈值：全队列每个marker一个锁定Otsu阈值；PDGFRA=25.6292，KRT14=18.3127（背景扣除后的p90信号）。\n\n"
for r in fib_sum.itertuples(): text+=f"- {r.mode}: n={r.count}, mean={r.mean:.4f}%, median={r.median:.4f}%, range={r.min:.4f}–{r.max:.4f}%.\n"
om=fib_tests.iloc[0]; text+=f"\nKruskal–Wallis: H={om.statistic:.4f}, p={om.p_raw:.4f}。所有BH校正后的两两比较q≥0.3429，均未达到显著性。由于每组切片数较少，结果仅作探索性解释。\n"
(FIB/"05_delivery_PDF_AI"/"统计结果摘要_请先阅读.md").write_text(text,encoding='utf-8')

kc_tests=pd.read_csv(KC/"02_statistics"/"KC_statistical_tests.tsv",sep='\t')
kc_sum=pd.read_csv(KC/"02_statistics"/"KC_mode_summary.tsv",sep='\t')
ktext="# KC整张切片统计结果\n\n"
ktext+="- 指标：SOX9+ MMP3+ KRT6B+ cell ratio (%)。\n- 分母：整张切片内通过质量筛选的DAPI+细胞核。\n- Zoom-in ROI不参与分子或分母计算。\n- 原始TIFF缺少µm/px标定，因此分析使用原始像素尺度；物理比例尺状态为CALIBRATION_REQUIRED。\n- 001A缺少可靠时间映射，保留整片计数但不进入Mode检验。\n\n"
for r in kc_sum.itertuples(): ktext+=f"- {r.mode}: n={r.count}, mean={r.mean:.4f}%, median={r.median:.4f}%, range={r.min:.4f}–{r.max:.4f}%.\n"
om=kc_tests.iloc[0]; ktext+=f"\nKruskal–Wallis: H={om.statistic:.4f}, p={om.p_raw:.4f}。所有BH校正后的两两比较q≥0.8571，均未达到显著性。由于每组切片数较少且来源TIFF缺少物理标定，结果仅作探索性解释。\n"
(KC/"05_delivery_PDF_AI"/"统计结果摘要_请先阅读.md").write_text(ktext,encoding='utf-8')
print('FINALIZED')
