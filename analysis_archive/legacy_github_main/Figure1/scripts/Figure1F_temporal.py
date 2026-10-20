from __future__ import annotations

import csv
from collections import defaultdict

from scipy.stats import kruskal

from figure_svg_utils import ROOT, ensure_dirs, render_svg


DATA = ROOT / "source_data/Figure1F_section_level_scores.csv"
OUT = ROOT / "outputs/panels"
STEM = "Figure1F_v01"
WIDTH_MM, HEIGHT_MM = 190, 38.6
TIME = ["Normal", "5 dpb", "12 dpb", "19 dpb", "26 dpb", "2 mph"]
TICK = ["N", "5d", "12d", "19d", "26d", "2m"]
PROCESSES = [
    ("Epithelialization", "Epithelialization", "#0072B2"),
    ("Skin development", "Skin development", "#009E73"),
    ("Inflammatory response", "Inflammatory response", "#D55E00"),
    ("Vascularization", "Vascularization", "#56B4E9"),
    ("Matrix remodeling", "Matrix remodeling", "#CC79A7"),
]


def read_values() -> list[dict[str, str]]:
    with DATA.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def build_svg(rows: list[dict[str, str]]) -> str:
    panel_x = [8, 44.2, 80.4, 116.6, 152.8]
    panel_w = 32.2
    top, bottom = 11.2, 30.8
    y_min, y_max = -2.3, 2.3

    def y_pos(value: float) -> float:
        return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM}mm" height="{HEIGHT_MM}mm" viewBox="0 0 {WIDTH_MM} {HEIGHT_MM}">',
        '<rect width="190" height="38.6" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.panel{font-size:5px;font-weight:700}.head{font-size:3.2px;font-weight:700}.title{font-size:2.25px;font-weight:700}.tick{font-size:2.12px}.axis{stroke:#333;stroke-width:.28}.grid{stroke:#c8c8c8;stroke-width:.2;stroke-dasharray:1,1}</style>',
        '<text class="panel" x="0.8" y="5">F</text>',
        '<text class="head" x="7" y="4.8">Temporal dynamics across registered Stereo-seq sections</text>',
        '<text class="tick" x="2.2" y="21" transform="rotate(-90 2.2 21)" text-anchor="middle">Section-level score (z)</text>',
    ]
    for pidx, (column, title, color) in enumerate(PROCESSES):
        x0 = panel_x[pidx]
        parts.append(f'<text class="title" x="{x0 + panel_w/2:.2f}" y="9.2" text-anchor="middle">{title}</text>')
        parts.append(f'<line class="axis" x1="{x0}" y1="{top}" x2="{x0}" y2="{bottom}"/><line class="axis" x1="{x0}" y1="{bottom}" x2="{x0+panel_w}" y2="{bottom}"/>')
        for tick_y in (-2, 0, 2):
            yp = y_pos(tick_y)
            parts.append(f'<line class="grid" x1="{x0}" y1="{yp:.2f}" x2="{x0+panel_w}" y2="{yp:.2f}"/>')
            if pidx == 0:
                parts.append(f'<text class="tick" x="{x0-1.2}" y="{yp+0.7:.2f}" text-anchor="end">{tick_y}</text>')
        means = []
        grouped_values = []
        for tidx, time in enumerate(TIME):
            xp = x0 + 2.0 + tidx * (panel_w - 4.0) / 5
            vals = [float(row[f"{column} z"]) for row in rows if row["timepoint"] == time]
            grouped_values.append(vals)
            means.append((xp, sum(vals) / len(vals)))
            offsets = [0.0] if len(vals) == 1 else [(-0.55 + 1.1 * i / (len(vals)-1)) for i in range(len(vals))]
            for off, value in zip(offsets, vals):
                parts.append(f'<circle cx="{xp+off:.2f}" cy="{y_pos(value):.2f}" r="0.75" fill="white" stroke="{color}" stroke-width="0.45"/>')
            parts.append(f'<text class="tick" x="{xp:.2f}" y="33.5" text-anchor="middle">{TICK[tidx]}</text>')
        path = " ".join(("M" if i == 0 else "L") + f" {x:.2f} {y_pos(v):.2f}" for i, (x, v) in enumerate(means))
        parts.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="0.65"/>')
        for x, value in means:
            parts.append(f'<circle cx="{x:.2f}" cy="{y_pos(value):.2f}" r="0.52" fill="{color}"/>')
        statistic, p_value = kruskal(*grouped_values)
        stars = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        parts.append(f'<text class="tick" x="{x0+panel_w-0.8:.2f}" y="12.2" text-anchor="end">{stars}</text>')
    parts.extend([
        '<text class="tick" x="95" y="37.5" text-anchor="middle">Healing time point</text>',
        '<text class="tick" x="188.5" y="37.5" text-anchor="end">Kruskal–Wallis: *P&lt;0.05; **P&lt;0.01; ***P&lt;0.001; ns, not significant</text>',
        '</svg>',
    ])
    return "\n".join(parts)


def main() -> None:
    ensure_dirs()
    rows = read_values()
    out = OUT / f"{STEM}.svg"
    out.write_text(build_svg(rows), encoding="utf-8")
    render_svg(out, WIDTH_MM, HEIGHT_MM, STEM)
    with (ROOT / "source_data/Figure1F_Kruskal_Wallis.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["process", "test", "statistic", "p_value", "significance", "observational_unit", "group_n"])
        for column, title, color in PROCESSES:
            groups = [[float(row[f"{column} z"]) for row in rows if row["timepoint"] == time] for time in TIME]
            statistic, p_value = kruskal(*groups)
            stars = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            writer.writerow([title, "two-sided Kruskal-Wallis omnibus", statistic, p_value, stars, "Stereo-seq section", "1;1;4;2;2;1"])
    (ROOT / f"QC/{STEM}_content_proof.txt").write_text(
        "Figure 1F v01\n- Exact registered 11-section subset only.\n"
        "- Section-level means are the statistical units; individual section points are shown.\n"
        "- 6 mph and 9 mph excluded.\n- Within-process z standardization used for display.\n"
        "- Timepoint section counts: Normal 1; 5 dpb 1; 12 dpb 4; 19 dpb 2; 26 dpb 2; 2 mph 1.\n"
        "- Two-sided Kruskal-Wallis omnibus tests were added using sections as observational units; all five processes are ns.\n",
        encoding="utf-8",
    )
    print(out)


if __name__ == "__main__":
    main()
