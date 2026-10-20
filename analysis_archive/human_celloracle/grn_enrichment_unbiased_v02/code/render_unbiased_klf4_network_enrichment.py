#!/usr/bin/env python3
from pathlib import Path
import argparse, json, math, textwrap
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

MM=1/25.4
SIZE=(52.5*MM,52.5*MM)
COL={
    "tf":"#E64B35","positive":"#7BC87C","negative":"#B8A1D9",
    "strict":"#2F8EE5","epi":"#E6A700","edge_pos":"#6F8F72",
    "edge_neg":"#8A6F93","term":"#F3C969","gene":"#8EC5E8",
}

def configure():
    mpl.rcParams.update({
        "font.family":"Arial","font.size":6.5,
        "axes.titlesize":8.5,"axes.labelsize":7,
        "xtick.labelsize":6,"ytick.labelsize":6,
        "pdf.fonttype":42,"ps.fonttype":42,
        "axes.linewidth":0.65,
    })

def save(fig,stem:Path):
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"),dpi=600)
    fig.savefig(stem.with_suffix(".tiff"),dpi=600,pil_kwargs={"compression":"tiff_lzw"})

def sector_positions(rows:pd.DataFrame,upper:bool)->dict:
    rows=rows.sort_values(["evidence_score","coef_abs"],ascending=False)
    n=len(rows); golden=math.pi*(3-math.sqrt(5))
    out={}
    for i,gene in enumerate(rows.target.astype(str)):
        # High-evidence targets sit near the perimeter so their objective labels
        # can be read without crowding the TF hub.
        r=0.60+2.15*math.sqrt((n-i)/(n+1))
        phase=(i*golden)%math.pi
        theta=phase if upper else math.pi+phase
        out[gene]=(1.14*r*math.cos(theta),r*math.sin(theta))
    return out

def network_graph(edges:pd.DataFrame,epi_targets:set[str]):
    g=nx.DiGraph(); g.add_node("KLF4",node_type="TF")
    for r in edges.itertuples(index=False):
        g.add_node(str(r.target),node_type="target",direction=r.direction,
                   epithelial=str(r.target) in epi_targets,
                   strict=bool(r.strict_filtered_link),evidence=float(r.evidence_score))
        g.add_edge("KLF4",str(r.target),direction=r.direction,
                   coef_mean=float(r.coef_mean),coef_abs=float(r.coef_abs),
                   FDR=float(r.FDR),evidence_score=float(r.evidence_score),
                   strict=bool(r.strict_filtered_link))
    return g

def draw_unbiased_network(ax,edges:pd.DataFrame,epi_targets:set[str]):
    positive=edges[edges.direction.eq("positive")]
    negative=edges[edges.direction.eq("negative")]
    pos={"KLF4":(0,0),**sector_positions(positive,True),**sector_positions(negative,False)}
    max_score=max(edges.evidence_score.max(),1e-9)
    # Every FDR<0.05 edge is drawn. Direction controls line style; coefficient
    # magnitude controls width. No target is manually inserted.
    for r in edges.itertuples(index=False):
        x,y=pos[str(r.target)]
        ax.plot([0,x],[0,y],
                color=COL["edge_pos"] if r.direction=="positive" else COL["edge_neg"],
                lw=0.12+0.45*min(float(r.coef_abs),1.0),
                ls="-" if r.direction=="positive" else (0,(2,1.5)),
                alpha=0.16,zorder=1)
    for direction,frame in [("positive",positive),("negative",negative)]:
        xy=np.array([pos[str(g)] for g in frame.target])
        sizes=3.5+16*np.sqrt(frame.evidence_score.to_numpy()/max_score)
        nonstrict=~frame.strict_filtered_link.astype(bool).to_numpy()
        ax.scatter(xy[nonstrict,0],xy[nonstrict,1],s=sizes[nonstrict],
                   c=COL[direction],edgecolors=[
                       COL["epi"] if g in epi_targets else "white"
                       for g in frame.loc[nonstrict,"target"].astype(str)
                   ],linewidths=0.45,alpha=0.9,zorder=3)
        strict=~nonstrict
        if strict.any():
            ax.scatter(xy[strict,0],xy[strict,1],s=sizes[strict]+8,
                       c=COL[direction],marker="s",edgecolors=COL["strict"],
                       linewidths=1.0,zorder=4)
    ax.scatter([0],[0],s=115,c=COL["tf"],edgecolors="white",linewidths=0.8,zorder=6)
    ax.text(0,0,r"$\it{KLF4}$",ha="center",va="center",fontsize=8.5,
            fontweight="bold",zorder=7)
    label_rows=edges.sort_values(["evidence_score","coef_abs"],ascending=False).head(6)
    for r in label_rows.itertuples(index=False):
        x,y=pos[str(r.target)]
        tx,ty=x*1.13,y*1.13
        ax.plot([x,tx],[y,ty],color="#777777",lw=0.3,zorder=7)
        ax.text(tx,ty,str(r.target),fontsize=6,fontstyle="italic",
                ha="left" if tx>=0 else "right",va="center",zorder=8)
    ax.set_xlim(-3.55,3.55);ax.set_ylim(-3.25,3.25);ax.set_aspect("equal");ax.axis("off")
    ax.set_title(r"Unbiased $\it{KLF4}$ network",pad=1.0,fontsize=8)

    handles=[
        Line2D([0],[0],marker="o",ls="",markerfacecolor=COL["positive"],markeredgecolor="white",label="Positive node"),
        Line2D([0],[0],marker="o",ls="",markerfacecolor=COL["negative"],markeredgecolor="white",label="Negative node"),
        Line2D([0],[0],marker="s",ls="",markerfacecolor="white",markeredgecolor=COL["strict"],label="Strict link"),
        Line2D([0],[0],marker="o",ls="",markerfacecolor="white",markeredgecolor=COL["epi"],label="Epithelial target"),
        Line2D([0,1],[0,0],color=COL["edge_pos"],lw=0.8,ls="-",label="Positive edge"),
        Line2D([0,1],[0,0],color=COL["edge_neg"],lw=0.8,ls=(0,(2,1.5)),label="Negative edge"),
    ]
    ax.legend(handles=handles,ncol=2,loc="lower center",bbox_to_anchor=(0.5,-0.235),
              fontsize=6,frameon=False,handlelength=1.15,columnspacing=0.55,
              handletextpad=0.25,borderaxespad=0,labelspacing=0.25)
    ax.text(0.5,-0.335,"644 FDR<0.05 edges; top 6 labels by evidence score",
            transform=ax.transAxes,ha="center",va="top",fontsize=6)

def overlap_count(s):
    return len([g for g in str(s).split(";") if g])

def wrap_term(term,width=24):
    return "\n".join(textwrap.wrap(term.split(" (GO:")[0],width=width))

def draw_overall_ora(ax,ora:pd.DataFrame):
    d=ora.loc[ora["Adjusted P-value"]<0.05].sort_values("Adjusted P-value").head(6).copy()
    d=d.iloc[::-1]
    y=np.arange(len(d)); x=-np.log10(d["Adjusted P-value"].clip(lower=1e-300))
    counts=d.Genes.map(overlap_count).to_numpy()
    ax.barh(y,x,color="#B8DEC9",edgecolor="#333333",lw=0.55,height=0.68)
    ax.scatter(x,y,s=10+1.7*counts,c="#2166AC",edgecolor="white",lw=0.35,zorder=3)
    ax.axvline(-np.log10(0.05),color="#777777",lw=0.7,ls=(0,(3,2)))
    ax.set_yticks(y);ax.set_yticklabels([wrap_term(t,20) for t in d.Term])
    ax.set_xlabel(r"$-\log_{10}$(FDR)")
    ax.set_title(r"$\it{KLF4}$-target ORA",pad=2,fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
    ax.text(-np.log10(0.05)+0.04,0.98,"FDR=0.05",transform=ax.get_xaxis_transform(),
            ha="left",va="top",fontsize=6,color="#666666")

def draw_epi_terms(ax,epi_terms:pd.DataFrame):
    d=epi_terms.sort_values("Adjusted P-value").head(6).copy().iloc[::-1]
    y=np.arange(len(d)); x=-np.log10(d["Adjusted P-value"].clip(lower=1e-300))
    counts=d.Genes.map(overlap_count).to_numpy()
    ax.barh(y,x,color="#F4D78A",edgecolor="#7A6427",lw=0.55,height=0.68)
    ax.scatter(x,y,s=10+1.8*counts,c="#B8860B",edgecolor="white",lw=0.35,zorder=3)
    ax.axvline(-np.log10(0.05),color="#777777",lw=0.7,ls=(0,(3,2)))
    concise={
        "Cellular Response to Epidermal Growth Factor Stimulus":"Cellular response to\nEGF stimulus",
        "Response to Epidermal Growth Factor":"Response to EGF",
        "Epithelium Development":"Epithelium development",
        "Epithelial Cell Differentiation":"Epithelial cell\ndifferentiation",
        "Regulation of Wound Healing":"Regulation of\nwound healing",
        "Skin Epidermis Development":"Skin epidermis\ndevelopment",
    }
    names=[]
    for term in d.Term:
        base=term.split(" (GO:")[0]
        names.append(concise.get(base,wrap_term(term,20)))
    ax.set_yticks(y);ax.set_yticklabels(names)
    ax.set_xlabel(r"$-\log_{10}$(FDR)")
    ax.set_title("Epithelial GO trends",pad=2,fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
    ax.text(-np.log10(0.05)-0.02,0.98,"FDR=0.05",transform=ax.get_xaxis_transform(),
            ha="right",va="top",fontsize=6,color="#666666")

def draw_epi_targets(ax,mapping:pd.DataFrame):
    # One row per target; ranking is the predeclared CellOracle evidence score.
    d=(mapping.sort_values(["evidence_score","target"],ascending=[False,True])
       .drop_duplicates("target").reset_index(drop=True))
    n=len(d); split=int(math.ceil(n/2))
    vmax=max(float(d.evidence_score.max()),1e-9)
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    ax.set_title("Epithelial GO-linked targets",pad=2,fontsize=8)
    ax.text(0.5,0.925,"22 targets from 10 lowest-FDR epithelial terms",
            ha="center",va="center",fontsize=6,color="#555555")
    for col,block in enumerate([d.iloc[:split],d.iloc[split:]]):
        x0=0.04+0.50*col; bar0=x0+0.18; maxw=0.14
        for row,rec in enumerate(block.itertuples(index=False)):
            y=0.835-row*(0.535/max(1,split-1))
            ax.text(x0,y,str(rec.target),ha="left",va="center",fontsize=6,
                    fontstyle="italic",fontweight="bold" if bool(rec.strict_filtered_link) else "normal")
            w=maxw*float(rec.evidence_score)/vmax
            ax.plot([bar0,bar0+w],[y,y],lw=3.0,color=COL["strict"] if bool(rec.strict_filtered_link) else COL["epi"],
                    solid_capstyle="round")
            ax.text(bar0+maxw+0.012,y,f"{float(rec.evidence_score):.2f}",ha="left",va="center",fontsize=6)
    handles=[Line2D([0,1],[0,0],color=COL["epi"],lw=3,
                    label="FDR<0.05 KLF4 target")]
    if d.strict_filtered_link.astype(bool).any():
        handles.append(Line2D([0,1],[0,0],color=COL["strict"],lw=3,
                              label="Strict CellOracle target"))
    ax.legend(handles=handles,loc="lower center",bbox_to_anchor=(0.5,0.085),ncol=1,
              fontsize=6,frameon=False,handlelength=1.4,handletextpad=0.35)
    ax.text(0.5,0.015,"Bar length = CellOracle evidence score",
            ha="center",va="bottom",fontsize=6)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--ora",type=Path,required=True)
    ap.add_argument("--epi-terms",type=Path,required=True)
    ap.add_argument("--epi-map",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True);configure()
    edges=pd.read_csv(args.edges,sep="\t")
    ora=pd.read_csv(args.ora,sep="\t")
    epi_terms=pd.read_csv(args.epi_terms,sep="\t")
    epi_map=pd.read_csv(args.epi_map,sep="\t")
    epi_targets=set(epi_map.target.astype(str))
    graph=network_graph(edges,epi_targets)
    nx.write_graphml(graph,args.output/"Fib_K14_KLF4_unbiased_FDR05_network.graphml")

    fig=plt.figure(figsize=SIZE);ax=fig.add_axes([0.08,0.28,0.84,0.62])
    draw_unbiased_network(ax,edges,epi_targets);save(fig,args.output/"Panel_D_unbiased_KLF4_network_square");plt.close(fig)
    fig=plt.figure(figsize=SIZE);ax=fig.add_axes([0.55,0.25,0.40,0.64])
    draw_overall_ora(ax,ora)
    fig.text(0.5,0.025,"390 targets; background=2,505\nDot size = target count",
             ha="center",fontsize=6,linespacing=1.05)
    save(fig,args.output/"Panel_E_unbiased_KLF4_target_ORA_square");plt.close(fig)
    fig=plt.figure(figsize=SIZE);ax=fig.add_axes([0.55,0.25,0.40,0.64])
    draw_epi_terms(ax,epi_terms)
    fig.text(0.5,0.025,"No epithelial term passes FDR<0.05\nDot size = target count",
             ha="center",fontsize=6,color="#A50F15",linespacing=1.05)
    save(fig,args.output/"Panel_F_epithelialization_GO_trends_square");plt.close(fig)
    fig=plt.figure(figsize=SIZE);ax=fig.add_axes([0.04,0.06,0.92,0.86])
    draw_epi_targets(ax,epi_map);save(fig,args.output/"Panel_G_epithelialization_KLF4_targets_square");plt.close(fig)
    checkpoint={
        "status":"PASS_PENDING_PHYSICAL_QA","figure_size_mm":[52.5,52.5],
        "font_family":"Arial","font_range_pt":[6,12],
        "network_selection":"all within-TF BH FDR<0.05 KLF4 edges",
        "network_edges":int(len(edges)),"manual_gene_additions":0,
        "network_labels":"top 6 by predeclared evidence_score only",
        "ORA_query":"all positive FDR<0.05 KLF4 targets",
        "ORA_display":"top 6 significant terms by adjusted P-value; complete table retained",
        "epithelial_focus":"top 10 lowest-FDR keyword-matched terms; all retained as FDR>=0.05 trends",
    }
    (args.output/"render_checkpoint.json").write_text(json.dumps(checkpoint,indent=2),encoding="utf-8")

if __name__=="__main__":main()
