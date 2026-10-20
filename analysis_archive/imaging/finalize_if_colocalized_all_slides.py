from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUT = Path(r"data/external/local_5cc3be176607")


def merge_statistics() -> None:
    summaries, tests = [], []
    for family in ["KC", "Fib_K14"]:
        s = pd.read_csv(OUT / "05_statistics" / family / "group_summary.tsv", sep="\t")
        s.insert(0, "marker_family", family)
        summaries.append(s)
        t = pd.read_csv(OUT / "05_statistics" / family / "statistical_tests.tsv", sep="\t")
        t.insert(0, "marker_family", family)
        tests.append(t)
    pd.concat(summaries, ignore_index=True).to_csv(OUT / "05_statistics" / "group_summary.tsv", sep="\t", index=False)
    pd.concat(tests, ignore_index=True).to_csv(OUT / "05_statistics" / "statistical_tests.tsv", sep="\t", index=False)
    notes = {
        "status": "EXPLORATORY",
        "unit_of_analysis": "independent whole slide/subject; CUT_007 cropped derivative excluded from inference",
        "ratio_definition": "100 * co-positive DAPI-segmented cells / all accepted DAPI-segmented cells",
        "summary": "all raw whole-slide observations with group mean and SEM",
        "multiple_testing": "Benjamini-Hochberg within each eligible pairwise family",
        "limitations": [
            "Most group-by-mode cells contain fewer than two independent slides; tests are deliberately not run for those contrasts.",
            "32-1=SPTDI and 007=DPTDI remain author-review mappings.",
            "Automated segmentation and locked positivity thresholds remain exploratory pending biological image review.",
        ],
    }
    (OUT / "05_statistics" / "analysis_notes.json").write_text(json.dumps(notes, ensure_ascii=False, indent=2), encoding="utf-8")


def combined_figure() -> None:
    df = pd.read_csv(OUT / "04_per_slide" / "per_slide.tsv", sep="\t")
    df = df[df["inference_included"].astype(str).str.lower().eq("true")].copy()
    order = ["Mode 1", "Mode 2", "Mode 3"]
    colors = {"SPTDI": "#2D98BC", "DPTDI": "#E56351"}
    labels = {
        "KC": "SOX9+ MMP3+\nKRT6B+ DAPI+\ncell ratio (%)",
        "Fib_K14": "PDGFRA+ KRT14+\nDAPI+ cell ratio (%)",
    }
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(84/25.4, 36/25.4))
    rng = np.random.default_rng(20260825)
    for ax, family in zip(axes, ["KC", "Fib_K14"]):
        sub = df[df.family == family]
        for gi, group in enumerate(["SPTDI", "DPTDI"]):
            means, sems = [], []
            for xi, mode in enumerate(order):
                vals = sub[(sub.group == group) & (sub["mode"] == mode)].ratio_pct.to_numpy(dtype=float)
                means.append(float(np.mean(vals)) if len(vals) else np.nan)
                sems.append(float(pd.Series(vals).sem()) if len(vals) > 1 else 0.0)
                jitter = rng.uniform(-0.045, 0.045, len(vals)) + (-0.06 if group == "SPTDI" else 0.06)
                ax.scatter(xi + jitter, vals, s=15, facecolors="white", edgecolors=colors[group], linewidths=0.8, zorder=3)
            x = np.arange(3)
            ax.plot(x, means, color=colors[group], marker="o", markersize=3.5, linewidth=1.1, label=group)
            ax.fill_between(x, np.asarray(means)-np.asarray(sems), np.asarray(means)+np.asarray(sems), color=colors[group], alpha=0.13, linewidth=0)
        ax.set_xticks(range(3), ["Mode\n1", "Mode\n2", "Mode\n3"])
        ax.set_ylabel(labels[family])
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=0.8, length=3)
        ax.legend(frameon=False, loc="upper right")
    fig.tight_layout(pad=0.8, w_pad=1.3)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT / "06_figures" / f"KC_Fib_colocalized_cell_ratios.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def readme() -> None:
    text = """# IF共染细胞整片统计（2026-08-25）

- KC：11个图像文件均已统计；其中10张完整切片进入推断统计，CUT_007为007的裁剪图，仅作描述性结果。
- Fib_K14：10张KFBF对应的整片细胞表均已统计。
- KC指标：SOX9+ MMP3+ KRT6B+ DAPI+ cell ratio (%)。
- Fib指标：PDGFRA+ KRT14+ DAPI+ cell ratio (%)。
- 分母：声明区域内全部通过筛选的DAPI分割细胞核。
- 统计单位：独立整片/受试者；细胞和裁剪视野不作为生物学重复。
- 结果状态：EXPLORATORY。自动分割和阳性阈值仍需生物学图像复核；32-1=SPTDI、007=DPTDI仍需作者确认。
- 很多group×Mode组合不足2张独立片子，相应比较在统计表中标记为未检验，不生成虚假P值。
"""
    (OUT / "README_请先阅读.md").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    merge_statistics()
    combined_figure()
    readme()
    print(json.dumps({"status": "FINALIZED_TABLES_AND_FIGURES", "output": str(OUT)}, ensure_ascii=False))
