from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BLUE = "#2F98B9"
SALMON = "#ED6253"


def exact_permutation_pvalue(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    """Two-sided exact randomization test for a difference in means."""
    pooled = np.concatenate([a, b]).astype(float)
    n_a = len(a)
    observed = abs(float(np.mean(a) - np.mean(b)))
    total_sum = float(np.sum(pooled))
    n_total = len(pooled)
    extreme = 0
    total = 0
    tolerance = 1e-14
    for chosen in itertools.combinations(range(n_total), n_a):
        sum_a = float(np.sum(pooled[list(chosen)]))
        diff = abs(sum_a / n_a - (total_sum - sum_a) / (n_total - n_a))
        total += 1
        if diff >= observed - tolerance:
            extreme += 1
    return extreme / total, total


def holm_adjust(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * p[idx])
        adjusted[idx] = min(1.0, running)
    return adjusted.tolist()


def star_label(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def load_data(path: Path) -> tuple[pd.DataFrame, str]:
    data = pd.read_csv(path, sep="\t")
    data.columns = [str(c).strip() for c in data.columns]
    required = {"Time", "Grade", "MIF Fib_K14 proportion"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    sample_col = data.columns[0]
    data = data.rename(columns={sample_col: "Sample"})
    data["Sample"] = data["Sample"].astype(str).str.strip()
    data["Grade"] = data["Grade"].astype(str).str.strip()
    data["Time"] = pd.to_numeric(data["Time"], errors="raise").astype(int)
    data["ratio"] = pd.to_numeric(
        data["MIF Fib_K14 proportion"], errors="raise"
    )
    if data["Sample"].duplicated().any():
        duplicates = data.loc[data["Sample"].duplicated(), "Sample"].tolist()
        raise ValueError(f"Duplicate sample identifiers: {duplicates}")
    if data["ratio"].isna().any() or ((data["ratio"] < 0) | (data["ratio"] > 1)).any():
        raise ValueError("Ratio values must be complete and within [0, 1].")
    unexpected = sorted(set(data["Grade"]) - {"SPTDI", "DPTDI"})
    if unexpected:
        raise ValueError(f"Unexpected Grade labels: {unexpected}")
    data["ratio_percent"] = 100 * data["ratio"]
    return data, sample_col


def make_outputs(input_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    data, original_sample_col = load_data(input_path)

    summary = (
        data.groupby(["Time", "Grade"], sort=True)["ratio_percent"]
        .agg(n="count", mean="mean", sd="std", median="median", minimum="min", maximum="max")
        .reset_index()
    )
    summary["sem"] = summary["sd"] / np.sqrt(summary["n"])

    common_times = sorted(
        set(data.loc[data["Grade"] == "SPTDI", "Time"])
        & set(data.loc[data["Grade"] == "DPTDI", "Time"])
    )
    tests: list[dict[str, float | int | str]] = []
    raw_pvalues: list[float] = []
    for time in common_times:
        a = data.loc[
            (data["Time"] == time) & (data["Grade"] == "SPTDI"), "ratio_percent"
        ].to_numpy()
        b = data.loc[
            (data["Time"] == time) & (data["Grade"] == "DPTDI"), "ratio_percent"
        ].to_numpy()
        pvalue, permutations = exact_permutation_pvalue(a, b)
        raw_pvalues.append(pvalue)
        tests.append(
            {
                "Time": time,
                "SPTDI_n": len(a),
                "DPTDI_n": len(b),
                "SPTDI_mean_percent": float(np.mean(a)),
                "DPTDI_mean_percent": float(np.mean(b)),
                "mean_difference_SPTDI_minus_DPTDI_percent_points": float(
                    np.mean(a) - np.mean(b)
                ),
                "exact_permutation_p": pvalue,
                "enumerated_partitions": permutations,
            }
        )
    adjusted = holm_adjust(raw_pvalues)
    for row, p_adj in zip(tests, adjusted):
        row["Holm_adjusted_p"] = p_adj
        row["significance"] = star_label(p_adj)
    tests_df = pd.DataFrame(tests)

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.linewidth": 1.1,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, ax = plt.subplots(figsize=(6.7, 4.8), constrained_layout=False)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.82, bottom=0.29)

    rng = np.random.default_rng(20260825)
    group_specs = [("SPTDI", BLUE, -0.11), ("DPTDI", SALMON, 0.11)]
    all_times = sorted(data["Time"].unique())
    for grade, color, offset in group_specs:
        group_summary = summary.loc[summary["Grade"] == grade].sort_values("Time")
        x = group_summary["Time"].to_numpy(dtype=float)
        mean = group_summary["mean"].to_numpy(dtype=float)
        sem = group_summary["sem"].fillna(0).to_numpy(dtype=float)
        ax.fill_between(x, np.maximum(0, mean - sem), mean + sem, color=color, alpha=0.15, linewidth=0)
        ax.plot(x, mean, color=color, linewidth=2.4, marker="o", markersize=7.0,
                markerfacecolor=color, markeredgecolor=color, label=grade, zorder=3)
        for time in all_times:
            values = data.loc[
                (data["Grade"] == grade) & (data["Time"] == time), "ratio_percent"
            ].to_numpy()
            if len(values) == 0:
                continue
            jitter = rng.uniform(-0.20, 0.20, size=len(values))
            ax.scatter(
                np.full(len(values), time + offset) + jitter,
                values,
                s=42,
                facecolors="white",
                edgecolors=color,
                linewidths=1.35,
                alpha=0.95,
                zorder=4,
            )

    maximum = float(data["ratio_percent"].max())
    bracket_base = maximum + 0.9
    bracket_step = 0.65
    for i, row in tests_df.iterrows():
        x = float(row["Time"])
        local_max = float(data.loc[data["Time"] == int(x), "ratio_percent"].max())
        y = max(local_max + 0.7, bracket_base + (i % 2) * bracket_step)
        ax.plot([x - 0.34, x - 0.34, x + 0.34, x + 0.34],
                [y - 0.12, y, y, y - 0.12], color="black", linewidth=0.8, clip_on=False)
        label = str(row["significance"])
        if label == "ns":
            label = f"ns\nq={float(row['Holm_adjusted_p']):.3g}"
        else:
            label = f"{label}\nq={float(row['Holm_adjusted_p']):.3g}"
        ax.text(x, y + 0.08, label, ha="center", va="bottom", fontsize=6.5, linespacing=0.9)

    y_top = max(bracket_base + 2.0, maximum + 3.1)
    ax.set_ylim(0, y_top)
    ax.set_xlim(min(all_times) - 1.0, max(all_times) + 1.0)
    ax.set_xticks(all_times)
    ax.set_xlabel("Time", fontsize=11, labelpad=7)
    ax.set_ylabel(r"PDGFRA$^{+}$ KRT14$^{+}$ cell ratio (%)", fontsize=11, labelpad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.01, 1.08), ncol=2,
              fontsize=9, handlelength=2.0, columnspacing=1.2)

    n_lookup = summary.pivot(index="Grade", columns="Time", values="n")
    for row_index, (grade, color, _) in enumerate(group_specs):
        y_axes = -0.18 - row_index * 0.065
        ax.text(-0.03, y_axes, f"{grade} n", transform=ax.transAxes, color=color,
                ha="right", va="center", fontsize=7.2, fontweight="bold", clip_on=False)
        x_min, x_max = ax.get_xlim()
        for time in all_times:
            n_value = n_lookup.loc[grade, time] if time in n_lookup.columns else np.nan
            label = "–" if pd.isna(n_value) else str(int(n_value))
            x_axes = (time - x_min) / (x_max - x_min)
            ax.text(x_axes, y_axes, label, transform=ax.transAxes, color=color,
                    ha="center", va="center", fontsize=7.2, clip_on=False)

    ax.text(0.5, -0.315, "n = independent samples", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.8, color="#444444", clip_on=False)

    stem = output_dir / "PDGFRA_KRT14_cell_ratio_lineplot"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    data.rename(columns={"Sample": original_sample_col}).to_csv(
        output_dir / "source_data_with_percent.tsv", sep="\t", index=False
    )
    summary.to_csv(output_dir / "group_summary.tsv", sep="\t", index=False, float_format="%.8g")
    tests_df.to_csv(output_dir / "statistics_exact_permutation_Holm.tsv", sep="\t", index=False, float_format="%.8g")

    caption = (
        "PDGFRA⁺ KRT14⁺ cell ratio across time in SPTDI and DPTDI samples. "
        "Open circles represent independent samples; filled circles and connecting lines show group means; "
        "shaded bands indicate mean ± SEM. Sample sizes are shown below the x-axis for every group and time point. "
        "At time points containing both groups, two-sided exact permutation tests of the difference in means were used, "
        "with Holm correction across the five displayed comparisons. q denotes the Holm-adjusted p value; "
        "ns, q ≥ 0.05; *q < 0.05; **q < 0.01; ***q < 0.001. No observations were excluded."
    )
    (output_dir / "figure_legend.txt").write_text(caption + "\n", encoding="utf-8")

    qa = [
        f"Input: {input_path}",
        f"Rows: {len(data)}",
        f"Unique samples: {data['Sample'].nunique()}",
        f"Groups: {', '.join(sorted(data['Grade'].unique()))}",
        f"Times: {', '.join(map(str, all_times))}",
        f"Observed ratio range: {data['ratio_percent'].min():.4f}% to {data['ratio_percent'].max():.4f}%",
        f"Shared times tested: {', '.join(map(str, common_times))}",
        "Analysis unit: one row / unique sample identifier.",
        "Missing group-time combinations are displayed as en dashes in the n table and were not tested.",
        "Random jitter uses fixed seed 20260825 and does not change values.",
    ]
    (output_dir / "QA.txt").write_text("\n".join(qa) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    make_outputs(args.input, args.output_dir)


if __name__ == "__main__":
    main()
