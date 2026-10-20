"""Slice-level exploratory comparisons; not a stdlib statistics module."""
from pathlib import Path
import json, csv
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
rows=json.loads((ROOT/'source_data/healing_analysis_input.json').read_text('utf8'))
out=[]
for exp,treat,control in [('KLF4','KLF4-OE','OEC'),('IGFL1','IGFL1','NC')]:
    rr=[r for r in rows if r['experiment']==exp];tests=[]
    for day in [3,5,7,10,12]:
        a=np.array([r['healing_rate_percent'] for r in rr if r['time_dpb']==day and r['group']==treat and r['healing_rate_percent'] is not None])
        b=np.array([r['healing_rate_percent'] for r in rr if r['time_dpb']==day and r['group']==control and r['healing_rate_percent'] is not None])
        t=stats.ttest_ind(a,b,equal_var=False,alternative='two-sided')
        se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
        va,vb=a.var(ddof=1)/len(a),b.var(ddof=1)/len(b)
        df=(va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
        diff=a.mean()-b.mean();ci=stats.t.ppf(.975,df)*se
        tests.append(dict(experiment=exp,time_dpb=day,treatment=treat,control=control,n_treatment_slices=len(a),n_control_slices=len(b),difference_percentage_points=diff,CI95_low=diff-ci,CI95_high=diff+ci,welch_t=t.statistic,welch_df=df,p_two_sided=t.pvalue,analysis_status='EXPLORATORY conditional on slice independence; animal clustering unresolved; no confirmatory significance claim'))
    corrected=multipletests([v['p_two_sided'] for v in tests],method='holm')[1]
    for r,p in zip(tests,corrected):r['p_holm_five_timepoints']=p;r['significance_candidate']='ns' if p>=.05 else '*' if p>=.01 else '**' if p>=.001 else '***'
    out+=tests
with (ROOT/'source_data/healing_Welch_Holm_EXPLORATORY.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,list(out[0]));w.writeheader();w.writerows(out)
(ROOT/'source_data/healing_Welch_Holm_EXPLORATORY.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),'utf8')
for r in out:print(r['experiment'],r['time_dpb'],'difference',round(r['difference_percentage_points'],3),'p_raw',round(r['p_two_sided'],6),'p_holm',round(r['p_holm_five_timepoints'],6),flush=True)
