from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


BLUE = "#2F98B9"
SALMON = "#ED6253"
MODE_ORDER = ["Mode1", "Mode2", "Mode3"]


def assign_mode(time: int) -> str:
    if 0 <= time <= 7:
        return "Mode1"
    if 8 <= time <= 14:
        return "Mode2"
    if 15 <= time <= 28:
        return "Mode3"
    raise ValueError(f"Time {time} dpi is outside the requested 0–28 dpi range.")


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


def stars(p: float) -> str:
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
    original_sample_col = data.columns[0]
    data = data.rename(columns={original_sample_col: "Sample"})
    data["Sample"] = data["Sample"].astype(str).str.strip()
    data["Time"] = pd.to_numeric(data["Time"], errors="raise").astype(int)
    data["Grade"] = data["Grade"].astype(str).str.strip()
    data["ratio"] = pd.to_numeric(data["MIF Fib_K14 proportion"], errors="raise")
    if data["Sample"].duplicated().any():
        raise ValueError("Sample identifiers are not unique.")
    if ((data["ratio"] < 0) | (data["ratio"] > 1) | data["ratio"].isna()).any():
        raise ValueError("Ratio values must be complete and within [0, 1].")
    if set(data["Grade"]) != {"SPTDI", "DPTDI"}:
        raise ValueError(f"Unexpected Grade labels: {sorted(data['Grade'].unique())}")
    data["Mode"] = data["Time"].map(assign_mode)
    data["Mode"] = pd.Categorical(data["Mode"], categories=MODE_ORDER, ordered=True)
    data["ratio_percent"] = 100.0 * data["ratio"]
    return data, original_sample_col


def analyze(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = (
        data.groupby(["Mode", "Grade"], observed=True, sort=True)["ratio_percent"]
        .agg(n="count", mean="mean", sd="std", median="median", minimum="min", maximum="max")
        .reset_index()
    )
    summary["sem"] = summary["sd"] / np.sqrt(summary["n"])

    tests: list[dict[str, float | int | str]] = []
    raw_p: list[float] = []
    for mode in MODE_ORDER:
        a = data.loc[(data["Mode"] == mode) & (data["Grade"] == "SPTDI"), "ratio_percent"].to_numpy()
        b = data.loc[(data["Mode"] == mode) & (data["Grade"] == "DPTDI"), "ratio_percent"].to_numpy()
        result = mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
        u = float(result.statistic)
        p = float(result.pvalue)
        raw_p.append(p)
        tests.append(
            {
                "Mode": mode,
                "dpi_range": {"Mode1": "0-7", "Mode2": "8-14", "Mode3": "15-28"}[mode],
                "SPTDI_n": len(a),
                "DPTDI_n": len(b),
                "SPTDI_mean_percent": float(np.mean(a)),
                "DPTDI_mean_percent": float(np.mean(b)),
                "mean_difference_SPTDI_minus_DPTDI_percent_points": float(np.mean(a) - np.mean(b)),
                "Mann_Whitney_U": u,
                "rank_biserial_SPTDI_vs_DPTDI": float(2 * u / (len(a) * len(b)) - 1),
                "raw_p": p,
            }
        )
    adjusted = holm_adjust(raw_p)
    for row, q in zip(tests, adjusted):
        row["Holm_adjusted_p"] = q
        row["significance"] = stars(q)
    return summary, pd.DataFrame(tests)


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "font.style": "normal",
            "axes.labelsize": 8,
            "axes.labelweight": "normal",
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 2.8,
            "ytick.major.size": 2.8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def draw_figure(data: pd.DataFrame, summary: pd.DataFrame, tests: pd.DataFrame, output_stem: Path) -> None:
    configure_matplotlib()
    width_mm = 42.0
    height_mm = 40.0
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    fig.subplots_adjust(left=0.36, right=0.985, top=0.91, bottom=0.23)

    x_positions = np.arange(3, dtype=float)
    rng = np.random.default_rng(20260825)
    specs = [("SPTDI", BLUE, -0.075), ("DPTDI", SALMON, 0.075)]
    for grade, color, offset in specs:
        subset = (
            summary.loc[summary["Grade"] == grade]
            .set_index("Mode")
            .loc[MODE_ORDER]
            .reset_index()
        )
        mean = subset["mean"].to_numpy(dtype=float)
        sem = subset["sem"].to_numpy(dtype=float)
        ax.fill_between(x_positions, np.maximum(0, mean - sem), mean + sem,
                        color=color, alpha=0.15, linewidth=0, zorder=1)
        ax.plot(x_positions, mean, color=color, linewidth=1.35, marker="o", markersize=4.2,
                markerfacecolor=color, markeredgewidth=0.8, markeredgecolor=color,
                label=grade, zorder=3)
        for index, mode in enumerate(MODE_ORDER):
            values = data.loc[(data["Mode"] == mode) & (data["Grade"] == grade), "ratio_percent"].to_numpy()
            jitter = rng.uniform(-0.105, 0.105, len(values))
            ax.scatter(np.full(len(values), index + offset) + jitter, values,
                       s=14, facecolors="white", edgecolors=color, linewidths=0.8,
                       zorder=4, clip_on=True)

    observed_max = float(data["ratio_percent"].max())
    locator = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 4, 5, 8, 10])
    y_ticks = locator.tick_values(0, observed_max * 1.16)
    y_ticks = y_ticks[y_ticks >= 0]
    y_max = float(y_ticks[-1])
    annotation_pad = y_max * 0.045
    for index, row in tests.iterrows():
        if float(row["Holm_adjusted_p"]) >= 0.05:
            continue
        mode = str(row["Mode"])
        local_max = float(data.loc[data["Mode"] == mode, "ratio_percent"].max())
        y = min(local_max + annotation_pad, y_max - 2.1 * annotation_pad)
        ax.plot([index - 0.15, index - 0.15, index + 0.15, index + 0.15],
                [y - 0.35 * annotation_pad, y, y, y - 0.35 * annotation_pad],
                color="black", linewidth=0.65, clip_on=False)
        ax.text(index, y + 0.2 * annotation_pad, str(row["significance"]), ha="center", va="bottom",
                fontsize=10, fontweight="bold", fontstyle="normal")

    ax.set_xticks(x_positions, MODE_ORDER)
    ax.set_xlim(-0.32, 2.32)
    ax.set_ylim(0, y_max)
    ax.set_yticks(y_ticks)
    ax.set_ylabel("PDGFRA+ KRT14+\ncell ratio (%)", labelpad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", pad=2.0)
    ax.tick_params(axis="y", pad=1.2)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(-0.02, 1.015),
              ncol=1, handlelength=1.35, borderaxespad=0, labelspacing=0.25,
              handletextpad=0.35)

    fig.savefig(output_stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(output_stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(output_stem.with_suffix(".eps"), facecolor="white")
    fig.savefig(output_stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(output_stem.with_suffix(".tiff"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def make_a4_proof(panel_pdf: Path, output_pdf: Path) -> None:
    try:
        from pypdf import PdfReader, PdfWriter, PageObject, Transformation
    except ImportError:
        return
    mm_to_pt = 72.0 / 25.4
    a4_width = 210 * mm_to_pt
    a4_height = 297 * mm_to_pt
    x0 = 15 * mm_to_pt
    top_margin = 15 * mm_to_pt
    panel_height = 40 * mm_to_pt
    source_page = PdfReader(str(panel_pdf)).pages[0]
    page = PageObject.create_blank_page(width=a4_width, height=a4_height)
    page.merge_transformed_page(
        source_page,
        Transformation().translate(tx=x0, ty=a4_height - top_margin - panel_height),
    )
    writer = PdfWriter()
    writer.add_page(page)
    with output_pdf.open("wb") as handle:
        writer.write(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data, original_sample_col = load_data(args.input)
    summary, tests = analyze(data)
    stem = args.output_dir / "PDGFRA_KRT14_Mode1-3_42mm"
    draw_figure(data, summary, tests, stem)
    make_a4_proof(stem.with_suffix(".pdf"), args.output_dir / "PDGFRA_KRT14_Mode1-3_A4_placement_proof.pdf")

    data.rename(columns={"Sample": original_sample_col}).to_csv(
        args.output_dir / "source_data_with_Mode.tsv", sep="\t", index=False
    )
    summary.to_csv(args.output_dir / "Mode_group_summary.tsv", sep="\t", index=False, float_format="%.8g")
    tests.to_csv(args.output_dir / "Mode_statistics_MannWhitney_Holm.tsv", sep="\t", index=False, float_format="%.8g")

    caption = (
        "PDGFRA+ KRT14+ cell ratio in SPTDI and DPTDI samples across Mode1 (0–7 dpi), "
        "Mode2 (8–14 dpi), and Mode3 (15–28 dpi). Open circles represent individual samples; "
        "filled circles and connecting lines show group means; shaded bands indicate mean ± SEM. "
        "Sample sizes are provided in Source Data. Two-sided unpaired Mann–Whitney U tests "
        "were used for within-Mode comparisons, with Holm correction across the three comparisons. "
        "Exact adjusted p values and effect sizes are provided in Source Data; ns, adjusted p ≥ 0.05; "
        "*adjusted p < 0.05; **adjusted p < 0.01; ***adjusted p < 0.001. No observations were excluded."
    )
    (args.output_dir / "figure_legend.txt").write_text(caption + "\n", encoding="utf-8")

    qa = [
        f"Input: {args.input}",
        f"Rows / unique sample IDs: {len(data)} / {data['Sample'].nunique()}",
        "Mode mapping: Mode1=0-7 dpi; Mode2=8-14 dpi; Mode3=15-28 dpi.",
        "Standalone PDF media size: 42 mm x 40 mm; width is exactly 1/5 of 210 mm A4 width.",
        "Font family: Arial. Base/tick/axis/legend fonts: 8 pt. Significance symbols: 10 pt bold.",
        "Protein names are upright, consistent with protein nomenclature; no gene-symbol italics were used.",
        "Lines/markers/shading: means, individual samples, and mean ± SEM, respectively.",
        "Analysis unit: one row / unique sample identifier.",
        "Caution: pooled within-Mode comparisons assume rows are independent samples and do not model exact dpi.",
    ]
    (args.output_dir / "QA.txt").write_text("\n".join(qa) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
