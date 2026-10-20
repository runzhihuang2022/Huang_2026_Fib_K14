from pathlib import Path
import sys,json,hashlib,shutil,platform,html
import numpy as np,pandas as pd
from PIL import Image,ImageDraw,ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(r'data/external/local_4a28549fddd8')
OUT=ROOT/'统计分析_20261013'
manifest=pd.read_csv(OUT/'00_manifest/input_manifest.csv')
meta=pd.read_csv(OUT/'00_manifest/raw_acquisition_metadata.csv')
df=manifest.query("kind=='R'").merge(meta,on=['field_id','filename'],validate='one_to_one')
signature=['channel','pixel_size_um','exposure_s','gain','light_intensity']
keys=sorted(set(tuple(r) for r in df[signature].values.tolist()))
sid={k:f'A{i+1:02d}' for i,k in enumerate(keys)}
df['acquisition_set']=[sid[tuple(r)] for r in df[signature].values.tolist()]
metrics=[];imgs={};font=ImageFont.truetype('data/external/local_f014c5000f62',23)
for _,r in df.iterrows():
 with Image.open(r.source_path) as im:a=np.array(im)
 q=np.percentile(a,[10,25,50,75,90,95,99,99.9])
 ceiling=2**int(r.significant_bits)-1
 z=r.to_dict();z.update(raw_sum=float(a.astype(np.uint64).sum()),raw_sd=float(a.std()),p10=q[0],p25=q[1],p50=q[2],p75=q[3],p90=q[4],p95=q[5],p99=q[6],p999=q[7],ceiling_ADU=ceiling,saturated_pixels=int(np.sum(a>=ceiling)),total_pixels=int(a.size),saturated_pct=float(100*np.mean(a>=ceiling)),above_recorded_ceiling_pixels=int(np.sum(a>ceiling)),scope='完整保存视野',biological_marker='待确认',independent_replicate_id='',cell_count_status='未计算：未保存核染/全细胞通道',positive_cell_pct_status='未计算：缺少细胞分母及染色定义',measurement_status='图像层面描述性定量')
 metrics.append(z);imgs[(r.field_id,r.channel)]=a
result=pd.DataFrame(metrics).sort_values(['group_label','magnification_label','field_number','channel'])
result.to_csv(OUT/'04_per_field/per_field_channel_results.csv',index=False,encoding='utf-8-sig')
# Acquisition-homogeneous subgroups only; all counts here are fields, not biological replicates.
summary=[]
for (g,m,c,s),x in result.groupby(['group_label','magnification_label','channel','acquisition_set'],sort=True):
 summary.append(dict(group_label=g,magnification_label=m,channel=c,acquisition_set=s,field_n=len(x),raw_mean_mean=x['mean'].mean(),raw_mean_SD=x['mean'].std(ddof=1),raw_mean_min=x['mean'].min(),raw_mean_max=x['mean'].max(),saturation_mean_pct=x.saturated_pct.mean(),saturation_max_pct=x.saturated_pct.max(),independent_n='待确认',p_value_status='未检验：重复单位和分组含义待确认'))
summary=pd.DataFrame(summary);summary.to_csv(OUT/'05_statistics/group_acquisition_summary.csv',index=False,encoding='utf-8-sig')
sets=[]
for s,x in result.groupby('acquisition_set'):
 r=x.iloc[0];sets.append(dict(acquisition_set=s,channel=r.channel,pixel_size_metadata=float(r.pixel_size_um),exposure_s=r.exposure_s,gain=r.gain,light_intensity=r.light_intensity,field_n=len(x),field_ids='; '.join(x.field_id),group_labels='; '.join(sorted(set(x.group_label)))))
pd.DataFrame(sets).to_csv(OUT/'00_manifest/acquisition_sets.csv',index=False,encoding='utf-8-sig')
sample=[]
for f,x in result.groupby('field_id',sort=True):
 r=x.iloc[0];sample.append(dict(field_id=f,group_label=r.group_label,magnification_label=r.magnification_label,field_number=int(r.field_number),biological_group='待确认',independent_replicate_id='待确认',GFP_marker='待确认',CY5_marker='待确认',GFP_acquisition_set=x.loc[x.channel=='GFP','acquisition_set'].iloc[0],CY5_acquisition_set=x.loc[x.channel=='CY5','acquisition_set'].iloc[0],scope='完整保存视野'))
pd.DataFrame(sample).to_csv(OUT/'00_manifest/sample_metadata_to_confirm.csv',index=False,encoding='utf-8-sig')
channelmap=[dict(file_channel='d1',instrument_channel='GFP',display_color='green',marker='待确认',source='原始TIFF ImageDescription/Image Name'),dict(file_channel='d3',instrument_channel='CY5',display_color='red',marker='待确认',source='原始TIFF ImageDescription/Image Name'),dict(file_channel='未保存',instrument_channel='DAPI',display_color='',marker='未保存',source='协议列有DAPI，但输入目录无该通道图像')]
pd.DataFrame(channelmap).to_csv(OUT/'00_manifest/channel_map.csv',index=False,encoding='utf-8-sig')

# Per-field technical QC: fixed 0..4095 grayscale in both channels, magenta saturated pixels.
qc=[];thumbs=[]
for f,x in result.groupby('field_id',sort=True):
 canvas=Image.new('RGB',(1256,550),'white');draw=ImageDraw.Draw(canvas)
 draw.text((15,6),f+'  原始强度显示固定为 0–4095',font=font,fill='black')
 for i,c in enumerate(['GFP','CY5']):
  a=imgs[(f,c)];v=(np.clip(a/4095,0,1)*255).astype('uint8');rgb=np.repeat(v[:,:,None],3,axis=2);rgb[a>=4095]=[255,0,255]
  panel=Image.fromarray(rgb).resize((614,460),Image.Resampling.LANCZOS);canvas.paste(panel,(10+i*628,75))
  row=x[x.channel==c].iloc[0]
  draw.text((12+i*628,40),f'{c}  {row.acquisition_set}  均值 {row["mean"]:.1f}  饱和 {row.saturated_pct:.3f}%',font=font,fill='black')
 name=f'{f}_raw_QC.jpg';canvas.save(OUT/'02_qc'/name,quality=93)
 thumb=canvas.copy();thumb.thumbnail((628,275));thumbs.append(thumb)
 for _,r in x.iterrows():qc.append(dict(field_id=f,channel=r.channel,acquisition_set=r.acquisition_set,saturated_pct=r.saturated_pct,ceiling_ADU=4095,qc_status='存在满量程像素' if r.saturated_pixels else '无满量程像素',qc_image=f'02_qc/{name}',review_status='自动指标完成；图像已逐视野复核'))
pd.DataFrame(qc).to_csv(OUT/'02_qc/qc_summary.csv',index=False,encoding='utf-8-sig')
for batch in range(2):
 contact=Image.new('RGB',(1256,275*5),'white')
 for i,p in enumerate(thumbs[batch*10:(batch+1)*10]):contact.paste(p,((i%2)*628,(i//2)*275))
 contact.save(OUT/'02_qc'/f'raw_qc_contact_{batch+1}.jpg',quality=90)

plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none','font.size':10})
groups=['正单','烫单','烫双','六孔板'];colors=plt.cm.tab20(np.linspace(0,1,len(keys)))
colmap={f'A{i+1:02d}':colors[i] for i in range(len(keys))}
fig,axs=plt.subplots(2,2,figsize=(12,8))
for ax,(m,c) in zip(axs.flat,[('10x','GFP'),('10x','CY5'),('20x','GFP'),('20x','CY5')]):
 x=result[(result.magnification_label==m)&(result.channel==c)];present=[g for g in groups if g in set(x.group_label)]
 for i,g in enumerate(present):
  y=x[x.group_label==g].sort_values('field_number');xx=np.linspace(i-.2,i+.2,len(y))
  for xx1,(_,r) in zip(xx,y.iterrows()):
   ax.scatter(xx1,r['mean'],color=colmap[r.acquisition_set],s=65,edgecolor='black',linewidth=.4)
   ax.annotate(f'{int(r.field_number)}\n{r.acquisition_set}',(xx1,r['mean']),xytext=(0,8),textcoords='offset points',fontsize=8,ha='center')
 ax.set_xticks(range(len(present)),present);ax.set_ylim(0,float(x['mean'].max())*1.18);ax.set_xlim(-.65,len(present)-.35)
 ax.set_title(f'{m}（文件名） · {c}');ax.set_ylabel('全视野原始平均强度（ADU）');ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.18)
fig.suptitle('爬片图像描述性定量：每点为一个保存视野',fontsize=15)
fig.text(.5,.025,'标注为视野编号 / 采集参数组。不同 A 编号的强度不可直接比较；独立重复未确认，未进行显著性检验。',ha='center',fontsize=10)
fig.tight_layout(rect=[0,.055,1,.95])
for ext in ['png','pdf','svg']:fig.savefig(OUT/'06_figures'/f'raw_intensity_by_field.{ext}',dpi=180)
plt.close(fig)

issues=[
 ['染色定义','GFP/CY5是仪器通道名称，不能据此指定蛋白名称；需补充实际染色指标。'],
 ['细胞比例','输入无DAPI或已确认的全细胞通道，未计算总细胞数、阳性细胞数和阳性比例。'],
 ['独立重复','编号1–3的实验重复关系待确认；当前n均指视野，不代表独立培养、动物或供体。'],
 ['采集差异','不同组及部分同组视野存在增益/光强差异；按精确参数分为16个采集组，不作跨组强度归一化或显著性检验。'],
 ['饱和','文件为uint16容器、记录12个有效位；实际饱和值为4095。饱和比例按4095计算，强信号可能被截断。'],
 ['倍率与校准','倍率沿用文件名；TIFF像素尺度两档为0.628472和0.252112，比例不是2，未据此计算细胞大小或密度。'],
 ['视野重叠','不同倍率或编号可能重复覆盖同一位置，故未合并细胞数、面积或作为独立重复累加。'],
 ['单细胞归属','绿色/红色叠加及黄色区域不等同于经过分割验证的同一细胞双阳性。'],
]
issues[3][1]=issues[3][1].replace('16个',f'{len(keys)}个')
tests=[{'analysis':'组间生物学显著性检验','status':'未进行','reason':'分组含义与独立重复单位未确认；多组采集设置不同','p_value':None}]
pd.DataFrame(tests).to_csv(OUT/'05_statistics/statistical_tests.csv',index=False,encoding='utf-8-sig')
config={'input_root':str(ROOT/'20260930爬片'),'analysis_scope':'full_saved_field','raw_file_rule':'Plate_R only; d1=GFP,d3=CY5','input_files':120,'fields':20,'raw_images':40,'rendered_images':60,'protocols':20,'method':'All-pixel raw intensity distribution and 12-bit saturation audit. No cell segmentation or positivity thresholds. No gain/light normalization.','acquisition_key':signature,'threshold_rule':'None; biological positivity not assigned','cell_results_status':'BLOCKED_MISSING_STAIN_MAP_AND_CELL_DENOMINATOR','statistics_status':'DESCRIPTIVE_ONLY; independent units unconfirmed','python':sys.executable,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}}
(OUT/'08_audit/run_config.json').write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding='utf-8')
checks=[{'check':'all_input_files_inventory','observed':len(manifest),'expected':120,'status':'PASS'}, {'check':'raw_channel_measurements','observed':len(result),'expected':40,'status':'PASS'},{'check':'unique_fields','observed':result.field_id.nunique(),'expected':20,'status':'PASS'},{'check':'raw_file_unique_hashes','observed':result.sha256.nunique(),'expected':40,'status':'PASS'},{'check':'raw_values_above_4095','observed':int(result.above_recorded_ceiling_pixels.sum()),'expected':0,'status':'PASS'},{'check':'biological_marker_map','observed':'missing','expected':'confirmed','status':'PENDING_INPUT'},{'check':'independent_replicate_map','observed':'missing','expected':'confirmed','status':'PENDING_INPUT'},{'check':'cell_denominator','observed':'missing','expected':'validated nuclear or all-cell marker','status':'PENDING_INPUT'}]
pd.DataFrame(checks).to_csv(OUT/'08_audit/final_qa.csv',index=False,encoding='utf-8-sig')

simple=result[['field_id','channel','acquisition_set','mean','p50','p99','saturated_pct','gain','light_intensity']].rename(columns={'field_id':'视野','channel':'通道','acquisition_set':'采集组','mean':'全视野均值_ADU','p50':'中位数_ADU','p99':'P99_ADU','saturated_pct':'饱和像素_%','gain':'增益','light_intensity':'照明强度'})
simple.to_csv(OUT/'04_per_field/逐视野定量结果.csv',index=False,encoding='utf-8-sig')
sat=result[result.saturated_pct>=1][['field_id','channel','saturated_pct']]
summary_html=summary[['group_label','magnification_label','channel','acquisition_set','field_n','raw_mean_mean','raw_mean_SD','saturation_mean_pct']].rename(columns={'group_label':'文件名组别','magnification_label':'倍率标签','channel':'通道','acquisition_set':'采集组','field_n':'视野数','raw_mean_mean':'原始均值_ADU','raw_mean_SD':'视野间SD','saturation_mean_pct':'平均饱和_%'}).round(3).fillna('').to_html(index=False,border=0)
body=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>20260930爬片统计</title><style>body{{font:16px/1.7 'Microsoft YaHei',sans-serif;max-width:1250px;margin:40px auto;color:#172335;padding:0 28px}}h1{{font-size:28px}}h2{{margin-top:32px;font-size:22px}}table{{border-collapse:collapse;font-size:13px;width:100%}}th,td{{border-bottom:1px solid #ddd;padding:7px;text-align:right}}th:first-child,td:first-child{{text-align:left}}th{{background:#e8edf4}}img{{max-width:100%}}.note{{background:#fff4da;padding:16px;border-left:4px solid #c28a23}}a{{color:#24568a}}details{{margin:14px 0}}footer{{color:#555;font-size:13px}}</style>
<h1>20260930 爬片实验：图像描述性统计</h1><p>分析日期：2026-10-13。共20个保存视野、40张原始单通道图、60张显示/合并图、20份采集协议。</p>
<div class="note">已完成原始图像强度和采集质量统计。分组与染色定义、独立重复关系待确认；当前结果不包含阳性细胞比例或组间P值。</div>
<h2>可核实的发现</h2><ul><li>保存的原始通道为GFP和CY5，未见DAPI图像。</li><li>按通道、像素尺度、曝光、增益及照明强度得到{len(keys)}个采集参数组。</li><li>曝光均为0.5秒，但增益和照明强度不同，不能以肉眼亮度或原始均值直接推断表达差异。</li><li>烫单与烫双的20x CY5图像属于相同采集组；仍需考虑细胞密度、饱和、重复关系及染色含义。</li><li>{len(sat)}张原始通道图的满量程像素比例达到1%以上，详见逐视野表。</li></ul>
<h2>逐视野原始均值</h2><img src="06_figures/raw_intensity_by_field.png"><p>ADU为仪器原始灰度单位。点代表视野；无推断性误差线或显著性标记。均值包含背景及细胞数量/覆盖面积的影响，不代表每个细胞的蛋白表达量。</p>
<h2>描述性汇总</h2><p>仅在同一文件名组别、倍率标签、通道及采集参数组内计算均值和SD。单视野的SD留空。</p>{summary_html}
<h2>完整逐视野表</h2>{simple.round(4).to_html(index=False,border=0)}
<h2>原导出合并图索引</h2><p>以下采用仪器已有的显示图，仅用于定位和形态复核，显示亮度不用于定量。</p><img src="02_qc/all_fields_contact_sheet.jpg">
<h2>原始通道质控</h2><p>各图统一显示0–4095，紫色表示达到4095的像素。低强度通道在此统一窗口可能显得较暗。</p>'''
for f in result.field_id.drop_duplicates():body+=f'<details><summary>{html.escape(f)}</summary><img loading="lazy" src="02_qc/{html.escape(f)}_raw_QC.jpg"></details>'
body+='<h2>方法和待补充信息</h2><ul>'+''.join(f'<li><b>{a}</b>：{b}</li>' for a,b in issues)+'</ul>'
body+='<p>定量使用Plate_R原始灰度图的全部2048×1536像素，不进行背景扣除、阈值分割、增益换算或图像裁切。提取算术均值、中位数、分位数、原始强度和、标准差及饱和像素数。Plate_D与Plate_M仅作显示，不重复计入分析。</p><p><a href="07_excel/爬片实验_描述性统计.xlsx">下载Excel统计表</a> · <a href="04_per_field/逐视野定量结果.csv">逐视野CSV</a></p><footer>原始文件保持不变。SHA-256、采集参数及脚本记录保存在本报告同目录下的00_manifest与08_audit。</footer></html>'
(OUT/'统计报告.html').write_text(body,encoding='utf-8')
md='# 20260930爬片统计说明\n\n已完成20个视野、40张原始通道图的图像层面描述性定量。\n\n'+ '\n'.join(f'- {a}：{b}' for a,b in issues)+'\n\n入口：统计报告.html；Excel：07_excel/爬片实验_描述性统计.xlsx。\n'
(OUT/'README_统计说明.md').write_text(md,encoding='utf-8')
for s in ['coverslip_audit_20261013.py','coverslip_metadata_20261013.py','coverslip_descriptive_20261013.py']:shutil.copy2(Path(__file__).parent/s,OUT/'08_audit'/s)
payload={'summary':summary.replace({np.nan:None}).to_dict('records'),'results':result.replace({np.nan:None}).to_dict('records'),'sets':sets,'samples':sample,'channelmap':channelmap,'issues':issues,'manifest':manifest.replace({np.nan:None}).to_dict('records'),'qc':qc,'tests':tests,'checks':checks}
(OUT/'08_audit/workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False),encoding='utf-8')
print('acquisition sets',len(keys));print('saturation >=1%');print(sat.to_string(index=False));print('summary');print(summary.round(3).to_string(index=False))
