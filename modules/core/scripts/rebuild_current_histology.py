"""Reproduce the source collection used in the October 16 figure, preserving adjudication flags."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];D=R/'data/processed/histology_current';T=R/'outputs/tables/histology_current';F=R/'outputs/figures/histology_current'
T.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(D/'current_sections.csv');assert d.section_id.is_unique
assert d.healing_rate_percent.dropna().between(0,100).all()
s=d.groupby(['experiment','group','time_dpb']).healing_rate_percent.agg(n_slices='count',mean_percent='mean',sd_percent='std').reset_index();s['sem_percent']=s.sd_percent/np.sqrt(s.n_slices)
expected=pd.read_csv(D/'healing_summary_mean_SD_SEM.csv');z=s.merge(expected,on=['experiment','group','time_dpb'],suffixes=('_new','_expected'),validate='one_to_one')
for col in ['n_slices','mean_percent','sd_percent','sem_percent']:assert np.allclose(z[col+'_new'],z[col+'_expected'],rtol=1e-10,atol=1e-10)
s.to_csv(T/'summary_mean_SD_SEM.csv',index=False)
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
results=[];reports=[]
for exp,groups in [('KLF4',['OEC','KLF4-OE']),('IGFL1',['NC','IGFL1'])]:
 fig,ax=plt.subplots(figsize=(3.7,3.1));rng=np.random.default_rng(20261017);sub=d[d.experiment==exp]
 for j,g in enumerate(groups):
  q=s[(s.experiment==exp)&(s.group==g)].sort_values('time_dpb');c=['#3A9FBF','#E86A5A'][j];x=np.arange(len(q))
  ax.plot(x,q.mean_percent,'-o',c=c,lw=1,label=g,ms=3)
  ax.fill_between(x,q.mean_percent-q.sem_percent,q.mean_percent+q.sem_percent,color=c,alpha=.17)
  for k,day in enumerate(q.time_dpb):
   v=sub[(sub.group==g)&(sub.time_dpb==day)].healing_rate_percent.dropna();ax.scatter(k+rng.normal(0,.035,len(v)),v,s=12,facecolors='none',edgecolors=c,lw=.6)
 for day in [3,5,7,10,12]:
  a=sub[(sub.group==groups[1])&(sub.time_dpb==day)].healing_rate_percent.dropna().to_numpy();b=sub[(sub.group==groups[0])&(sub.time_dpb==day)].healing_rate_percent.dropna().to_numpy()
  tt=stats.ttest_ind(a,b,equal_var=False);results.append(dict(experiment=exp,time_dpb=day,n_treatment=len(a),n_control=len(b),difference=float(a.mean()-b.mean()),p_two_sided=float(tt.pvalue)))
 ax.set(xlabel='Days post-burn',ylabel='H&E re-epithelialization (%)',xticks=range(5),xticklabels=['3','5','7','10','12'],ylim=(-3,108));ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,fontsize=8)
 fig.text(.5,.015,'Sections; mean ± SEM. '+('Provisional D12 calls retained.' if exp=='KLF4' else 'Day-12 annotation conflict retained.'),ha='center',fontsize=7);fig.tight_layout(rect=(0,.045,1,1))
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{exp}_current_sections.{ext}',dpi=600)
 plt.close(fig)
 reports.append(dict(experiment=exp,n_rows=len(sub),n_observed=int(sub.healing_rate_percent.notna().sum()),n_missing=int(sub.healing_rate_percent.isna().sum()),animal_level_inference=False))
tests=pd.DataFrame(results);tests['p_holm_five_timepoints']=np.nan
for exp in tests.experiment.unique():
 idx=tests.index[tests.experiment==exp];p=tests.loc[idx,'p_two_sided'].to_numpy();order=np.argsort(p);adj=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))));out=np.empty_like(adj);out[order]=adj;tests.loc[idx,'p_holm_five_timepoints']=out
ref=pd.read_csv(D/'healing_Welch_Holm_EXPLORATORY.csv');joined=tests.merge(ref,on=['experiment','time_dpb'],suffixes=('_new','_source'))
assert np.allclose(joined.p_holm_five_timepoints_new,joined.p_holm_five_timepoints_source,atol=1e-10)
tests['status']='EXPLORATORY conditional on unverified section independence; not confirmatory animal inference';tests.to_csv(T/'Welch_Holm_exploratory.csv',index=False)
print(d.groupby(['data_origin','adjudication_status'],dropna=False).size().to_string())
report={'summary_equals_source':True,'exploratory_p_equals_source':True,'reports':reports,'all_holm_above_005':bool((tests.p_holm_five_timepoints>=.05).all()),'limitations':['KLF4 D12 includes 11 provisional morphology calls at 100 percent and 3 missing OE slices','IGFL1 D12 one supplied 100 percent value differs from earlier length annotation','Independent-animal/wound nesting unresolved; selection rationale required for current subset'],'source_collection':'October 16 v05, 135 rows; not the earlier 151-row workbook'}
(T/'audit.json').write_text(json.dumps(report,indent=2),'utf8');print(json.dumps(report))
