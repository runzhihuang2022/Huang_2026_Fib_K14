"""Rebuild the accepted X>2 four-protein strip from audited processed values.
The four archival source modules are imported as plotting functions only; no original main() runs.
"""
from pathlib import Path
import sys,json,importlib.util
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
R=Path(__file__).resolve().parents[2];V=Path(__file__).parent/'vendor';D=R/'data/processed/spatial_protein';O=R/'outputs/figures/spatial_protein';A=R/'validation/spatial_protein'
def mod(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main():
 for p in [O,A]:p.mkdir(exist_ok=True,parents=True)
 base=mod('sp_base',V/'build_figF_proteomics_four_Xgt1_20261015.py');st=mod('sp_style',V/'locked_he_style.py');ref=mod('sp_ref',V/'reference_spatial_batch.py');h=mod('sp_h',V/'build_figH_native_HE_regular_20261015.py')
 h.PLOTS=O;h.AUDIT=A
 t=pd.read_csv(D/'all_2029_spots_Fib_K14_Xgt2.tsv',sep='\t');xy=t[['x_display','y_display']].to_numpy(float)
 gate=(t[[g+'_X' for g in h.GENES]].values>2).all(axis=1)
 assert len(t)==2029 and np.array_equal(gate,t.Fib_K14) and gate.sum()==657
 front=st.inside(xy,h.ROI);assert front.sum()==174 and (gate&front).sum()==72
 features={g:t[g+'_display_scaled'].to_numpy(float) for g in h.GENES};features['Co-expression']=t.coexpression_display_score.to_numpy(float)
 limits=ref.full_limits(xy);densities={}
 for feature in ['Co-expression',*base.DISPLAY_GENES]:
  for region,bounds in [('full',limits),('zoom',h.ROI)]:
   densities[(feature,region)]=ref.density(xy[gate],features[feature][gate],bounds,{'kde_max_points':0},1729)[0]
 he=np.array(Image.open(D/'HE_native_grid_free_epidermis_up.png').convert('RGB'));assert he.shape==(458,671,3)
 plt.rcParams.update({'font.family':'Arial','font.size':8,'font.weight':'normal','axes.titleweight':'normal','axes.labelweight':'normal','axes.titlesize':8.5,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'figure.titlesize':8.5,'figure.titleweight':'normal','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white'})
 stems=[h.render(base,st,he,xy,gate,features,densities,limits,standalone) for standalone in [False,True]]
 h.morphology_reference(st,he,limits)
 result={'n_grid':len(t),'n_gate':int(gate.sum()),'percentage':100*gate.mean(),'ROI_grid':int(front.sum()),'ROI_gate':int((gate&front).sum()),'ROI_percentage':100*gate[front].mean(),'source_values_unchanged':True,'source_shape_and_counts_verified':True,'unit':'vendor AI-reconstructed multicellular grid, not cell','figure_records':st.FIGURE_RECORDS,'text_QA':h.TEXT_QA,'stems':stems}
 st.dump(A/'rebuild_report.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['figure_records','text_QA']}))
if __name__=='__main__':main()
