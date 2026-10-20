from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];D=R/'data/processed/histology';T=R/'outputs/tables/histology';F=R/'outputs/figures/histology'
def main():
 T.mkdir(exist_ok=True,parents=True);F.mkdir(exist_ok=True,parents=True)
 plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
 reports=[]
 for exp,groups in [('KLF4',['OEC','KLF4-OE']),('IGFL1',['NC','IGFL1'])]:
  d=pd.read_csv(D/f'{exp}_sections.csv');assert d.section_id.is_unique
  assert d.epithelialization_percent.dropna().between(0,100).all()
  s=d.groupby(['day_post_burn','group']).epithelialization_percent.agg(n='count',mean='mean',median='median',sd='std',minimum='min',maximum='max').reset_index();s.to_csv(T/f'{exp}_section_summary.csv',index=False)
  fig,ax=plt.subplots(figsize=(4.3,3.4));rng=np.random.default_rng(20261016)
  for j,g in enumerate(groups):
   q=s[s.group==g];c=['#358bab','#e66d5d'][j];ax.plot(q.day_post_burn,q['mean'],'-o',c=c,lw=1,label=g)
   ax.fill_between(q.day_post_burn,(q['mean']-q.sd).clip(0,100),(q['mean']+q.sd).clip(0,100),color=c,alpha=.13)
   v=d[d.group==g];ax.scatter(v.day_post_burn+rng.normal(0,.07,len(v)),v.epithelialization_percent,s=9,facecolors='none',edgecolors=c,lw=.5)
  ax.set(xlabel='Days post-burn',ylabel='H&E re-epithelialization (%)',ylim=(-3,103),xticks=sorted(d.day_post_burn.unique()),title=f'{exp}: section-level observations')
  ax.legend(frameon=False,fontsize=8);fig.text(.5,.005,'Points: sections; line: mean; shading: ± SD. No animal-level P values.',ha='center',fontsize=7);fig.tight_layout(rect=(0,.045,1,1))
  for ext in ['png','pdf','svg']:fig.savefig(F/f'{exp}_sections.{ext}',dpi=600)
  plt.close(fig)
  reports.append({'experiment':exp,'n_sections':len(d),'nonmissing':int(d.epithelialization_percent.notna().sum()),'missing':int(d.epithelialization_percent.isna().sum()),'independent_animal_mapping':'unavailable','inferential_test':'none','source_equality':True})
 (T/'audit.json').write_text(json.dumps(reports,indent=2),encoding='utf8');print(json.dumps(reports))
if __name__=='__main__':main()
