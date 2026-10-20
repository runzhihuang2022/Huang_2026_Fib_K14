from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

BASE = Path(r"data/external/local_de73c97db013")
KC_SOURCE = BASE / "KC_Spinous_Mig_whole_slide_statistics_native_pixel_v3" / "02_statistics" / "KC_whole_slide_per_slide.tsv"
FIB_SOURCE = BASE / "Fib_K14_output" / "07_whole_slide_statistics" / "02_statistics" / "Fib_K14_whole_slide_per_slide.tsv"
OUT = BASE / "KC_Fib_SPTDI_DPTDI_whole_slide_revision_20260825"
for d in [OUT / "01_source_tables", OUT / "02_statistics", OUT / "03_lineplots", OUT / "04_audit"]:
    d.mkdir(parents=True, exist_ok=True)

GROUP_MAP = {
    "001A": ("SPTDI", 10, "001A matched to workbook sample 001; inferred identifier suffix"),
    "003": ("DPTDI", 10, "workbook DSDB"),
    "007": ("DPTDI", 7, "available WSI assigned to DSDB/DPTDI; workbook contains both SSDB and DSDB fields"),
    "32-1": ("SPTDI", 6, "available WSI assigned to SSDB/SPTDI; workbook contains both SSDB and DSDB fields"),
    "BS1": ("DPTDI", 12, "workbook DSDB"),
    "BW3": ("SPTDI", 19, "workbook SSDB"),
    "BW4": ("DPTDI", 19, "workbook DSDB"),
    "BW7": ("DPTDI", 22, "workbook DSDB"),
    "BW8": ("DPTDI", 22, "sample metadata screenshot: deep burn 22 d"),
    "BW9": ("DPTDI", 26, "workbook DSDB"),
    "W14-1": ("SPTDI", 14, "workbook SSDB"),
}


def mode(day: int) -> str:
    return "Mode 1" if day <= 7 else ("Mode 2" if day <= 14 else "Mode 3")


def add_meta(df: pd.DataFrame, family: str, ratio: str, numerator: str) -> pd.DataFrame:
    rows = []
    for r in df.itertuples(index=False):
        sid = str(r.sample)
        group, day, provenance = GROUP_MAP[sid]
        rows.append({"family": family, "sample": sid, "group": group, "day_dpb": day,
                     "mode": mode(day), "numerator_n": int(getattr(r, numerator)),
                     "denominator_n": int(r.nuclei_count), "ratio_pct": float(getattr(r, ratio)),
                     "group_mapping_provenance": provenance, "analysis_region": "whole_slide"})
    return pd.DataFrame(rows)


def analyze(df: pd.DataFrame, family: str) -> pd.DataFrame:
    order = ["Mode 1", "Mode 2", "Mode 3"]
    tests = []
    for group in ["SPTDI", "DPTDI"]:
        vals = [df.loc[(df.group == group) & (df["mode"] == m), "ratio_pct"].to_numpy() for m in order]
        if all(len(v) for v in vals) and len(np.unique(np.concatenate(vals))) > 1:
            h, p = kruskal(*vals)
            tests.append({"family": family, "test_family": f"{group}_within_mode", "contrast": "Mode 1 vs Mode 2 vs Mode 3",
                          "test": "Kruskal-Wallis", "statistic": h, "p_raw": p, "p_BH": p})
        pairs = []
        for i, a in enumerate(order):
            for b in order[i + 1:]:
                x = df.loc[(df.group == group) & (df["mode"] == a), "ratio_pct"].to_numpy()
                y = df.loc[(df.group == group) & (df["mode"] == b), "ratio_pct"].to_numpy()
                if len(x) and len(y):
                    u, p = mannwhitneyu(x, y, alternative="two-sided", method="auto")
                    pairs.append({"family": family, "test_family": f"{group}_within_mode", "contrast": f"{a} vs {b}",
                                  "test": "Mann-Whitney U", "statistic": u, "p_raw": p})
        if pairs:
            q = multipletests([x["p_raw"] for x in pairs], method="fdr_bh")[1]
            for x, adj in zip(pairs, q): x["p_BH"] = adj
            tests.extend(pairs)
    between = []
    for m in order:
        x = df.loc[(df.group == "SPTDI") & (df["mode"] == m), "ratio_pct"].to_numpy()
        y = df.loc[(df.group == "DPTDI") & (df["mode"] == m), "ratio_pct"].to_numpy()
        if len(x) and len(y):
            u, p = mannwhitneyu(x, y, alternative="two-sided", method="auto")
            between.append({"family": family, "test_family": "SPTDI_vs_DPTDI_by_mode", "contrast": f"SPTDI vs DPTDI at {m}",
                            "test": "Mann-Whitney U", "statistic": u, "p_raw": p})
    if between:
        q = multipletests([x["p_raw"] for x in between], method="fdr_bh")[1]
        for x, adj in zip(between, q): x["p_BH"] = adj
        tests.extend(between)
    return pd.DataFrame(tests)


def summary(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["family", "group", "mode"], observed=True)["ratio_pct"].agg(
        independent_n="count", mean="mean", sem="sem", median="median", minimum="min", maximum="max").reset_index()


def plot(df: pd.DataFrame, tests: pd.DataFrame, ylabel: str, stem: Path, legend_loc: str = "upper left"):
    order = ["Mode 1", "Mode 2", "Mode 3"]
    colors = {"SPTDI": "#3298B8", "DPTDI": "#E86858"}
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(42/25.4, 28/25.4))
    rng = np.random.default_rng(20260825)
    for group in ["SPTDI", "DPTDI"]:
        means, sems = [], []
        for i, m in enumerate(order):
            y = df.loc[(df.group == group) & (df["mode"] == m), "ratio_pct"].to_numpy()
            means.append(float(np.mean(y)) if len(y) else np.nan)
            sems.append(float(pd.Series(y).sem()) if len(y) > 1 else 0.0)
            ax.scatter(i + rng.uniform(-.045, .045, len(y)), y, s=13, facecolors="white",
                       edgecolors=colors[group], linewidths=.75, zorder=3)
        means = np.asarray(means); sems = np.asarray(sems)
        ax.plot(range(3), means, "-o", color=colors[group], lw=1.1, ms=3.2, label=group)
        ax.fill_between(range(3), means-sems, means+sems, color=colors[group], alpha=.14, lw=0)
    ax.set_xticks(range(3), ["Mode\n1", "Mode\n2", "Mode\n3"])
    ax.set_ylabel(ylabel, fontsize=8, labelpad=0)
    ax.legend(frameon=False, loc=legend_loc, ncol=1, fontsize=8, handlelength=.65,
              handletextpad=.25, labelspacing=.05, borderaxespad=.05)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(width=.8, length=3, labelsize=8)
    sig = tests[(tests.test_family == "SPTDI_vs_DPTDI_by_mode") & (tests.p_BH < .05)]
    if not sig.empty:
        ymax = max(df.ratio_pct.max(), .01)
        for r in sig.itertuples():
            idx = next(i for i, m in enumerate(order) if m in r.contrast)
            label = "***" if r.p_BH < .001 else ("**" if r.p_BH < .01 else "*")
            ax.text(idx, ymax*1.05, label, ha="center", va="bottom", fontsize=8, fontweight="bold")
        ax.set_ylim(top=ymax*1.18)
    fig.subplots_adjust(left=.40, right=.97, bottom=.28, top=.95)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(stem.with_suffix("."+ext), dpi=600 if ext == "png" else None)
    plt.close(fig)


kc = add_meta(pd.read_csv(KC_SOURCE, sep="\t", dtype={"sample": str}), "KC",
              "triple_positive_ratio_pct", "triple_positive_n")
fib = add_meta(pd.read_csv(FIB_SOURCE, sep="\t", dtype={"sample": str}), "Fib_K14",
               "double_positive_ratio_pct", "double_positive_n")
all_df = pd.concat([kc, fib], ignore_index=True)
all_df.to_csv(OUT / "01_source_tables" / "whole_slide_grouped_per_slide.tsv", sep="\t", index=False)
mapping = all_df[["family", "sample", "group", "day_dpb", "mode", "group_mapping_provenance"]].drop_duplicates()
mapping.to_csv(OUT / "01_source_tables" / "sample_group_mode_mapping.tsv", sep="\t", index=False)

kc_tests = analyze(kc, "KC"); fib_tests = analyze(fib, "Fib_K14")
tests = pd.concat([kc_tests, fib_tests], ignore_index=True)
tests.to_csv(OUT / "02_statistics" / "grouped_statistical_tests.tsv", sep="\t", index=False)
summary(all_df).to_csv(OUT / "02_statistics" / "grouped_summary.tsv", sep="\t", index=False)

plot(kc, kc_tests, "SOX9+ MMP3+\nKRT6B+\ncell ratio (%)", OUT / "03_lineplots" / "KC_SPTDI_DPTDI_whole_slide_sensitivity", "upper right")
plot(fib, fib_tests, "PDGFRA+ KRT14+\ncell ratio (%)", OUT / "03_lineplots" / "Fib_K14_SPTDI_DPTDI_whole_slide", "upper center")

legacy = BASE / "Figure3D_replot_source_package_20260820" / "01_附件折线图_数字化复现" / "source_data" / "digitized_values.csv"
legacy_raw = pd.read_csv(legacy)
legacy_kc = pd.DataFrame({
    "family": "KC_existing_whole_slide",
    "sample": [f"{g}_{int(t)}_{int(r)}" for g, t, r in zip(legacy_raw.group, legacy_raw.time_dpb, legacy_raw.replicate_digitized)],
    "group": legacy_raw.group,
    "day_dpb": legacy_raw.time_dpb.astype(int),
    "mode": legacy_raw.time_dpb.map({7: "Mode 1", 14: "Mode 2", 28: "Mode 3"}),
    "numerator_n": np.nan,
    "denominator_n": np.nan,
    "ratio_pct": legacy_raw.SOX9_MMP3_KRT6B_positive_ratio_percent.astype(float),
    "group_mapping_provenance": "user-confirmed previous KC whole-slide plot; numerical point digitized from supplied raster",
    "analysis_region": "whole_slide_user_confirmed_previous_plot",
})
legacy_kc.to_csv(OUT / "01_source_tables" / "KC_existing_whole_slide_digitized_points.tsv", sep="\t", index=False)
legacy_tests = analyze(legacy_kc, "KC_existing_whole_slide")
legacy_tests.to_csv(OUT / "02_statistics" / "KC_existing_whole_slide_statistical_tests.tsv", sep="\t", index=False)
summary(legacy_kc).to_csv(OUT / "02_statistics" / "KC_existing_whole_slide_summary.tsv", sep="\t", index=False)
plot(legacy_kc, legacy_tests, "SOX9+ MMP3+\nKRT6B+\ncell ratio (%)",
     OUT / "03_lineplots" / "KC_existing_whole_slide_SPTDI_DPTDI_main", "upper left")

(OUT / "04_audit" / "analysis_config.json").write_text(json.dumps({
    "scope": "whole slide", "denominator": "valid DAPI-positive nuclei", "groups": {"SSDB": "SPTDI", "DSDB": "DPTDI"},
    "mode_rule": {"Mode 1": "0-7 dpb", "Mode 2": "8-14 dpb", "Mode 3": "15-28 dpb"},
    "independent_unit": "one point/slide for the user-confirmed existing KC plot; one available whole-slide section for current-image sensitivity reanalysis",
    "legacy_digitized_plot_used_for_main_KC_replot": str(legacy),
    "mapping_warning": "32-1 and 007 each have both SSDB and DSDB fields in old workbook; current WSI pairing uses 32-1=SPTDI and 007=DPTDI pending author confirmation",
    "kc_scale_bar_calibration": "author-directed inherited project calibration 0.25284451246261597 um/px"
}, ensure_ascii=False, indent=2), encoding="utf-8")
print(OUT)
print(tests.to_string(index=False))
