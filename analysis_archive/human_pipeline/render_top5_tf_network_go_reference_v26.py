#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd


MM = 1 / 25.4
SQUARE = (52.5 * MM, 52.5 * MM)
COLORS = {
    "tf": "#E64B35",
    "positive": "#7BC87C",
    "negative": "#B8A1D9",
    "edge_positive": "#6F8F72",
    "edge_negative": "#8A6F93",
    "go": "#F4D78A",
    "go_edge": "#7A6427",
    "go_dot": "#B8860B",
}


def safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6.5,
            "axes.titlesize": 8.5,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.65,
            "axes.unicode_minus": False,
        }
    )


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=600)


def sector(rows: pd.DataFrame, upper: bool) -> dict[str, tuple[float, float]]:
    rows = rows.sort_values(["edge_weight", "target_symbol"], ascending=[False, True])
    n = len(rows)
    golden = math.pi * (3 - math.sqrt(5))
    output: dict[str, tuple[float, float]] = {}
    for i, gene in enumerate(rows.target_symbol.astype(str)):
        radius = 0.60 + 2.05 * math.sqrt((n - i) / (n + 1))
        phase = (i * golden) % math.pi
        theta = phase if upper else math.pi + phase
        output[gene] = (1.14 * radius * math.cos(theta), radius * math.sin(theta))
    return output


def network_figure(edges: pd.DataFrame, subtype: str, tf: str) -> plt.Figure:
    q = (
        edges.sort_values(["edge_weight", "target_symbol"], ascending=[False, True])
        .drop_duplicates("target_symbol")
        .head(28)
        .copy()
    )
    positive = q[q.direction.eq("positive")]
    negative = q[q.direction.eq("negative")]
    positions = {tf: (0.0, 0.0), **sector(positive, True), **sector(negative, False)}
    vmax = max(float(q.edge_weight.max()) if len(q) else 0.0, 1e-12)

    fig = plt.figure(figsize=SQUARE)
    ax = fig.add_axes([0.07, 0.29, 0.86, 0.56])
    for row in q.itertuples(index=False):
        x, y = positions[str(row.target_symbol)]
        color = COLORS["edge_positive"] if row.direction == "positive" else COLORS["edge_negative"]
        linestyle = "-" if row.direction == "positive" else (0, (2, 1.5))
        ax.annotate(
            "",
            xy=(x, y),
            xytext=(0, 0),
            arrowprops={"arrowstyle": "-|>", "lw": 0.35, "color": color, "linestyle": linestyle, "alpha": 0.55, "mutation_scale": 4.5},
            zorder=1,
        )

    for direction, frame in (("positive", positive), ("negative", negative)):
        if frame.empty:
            continue
        xy = np.array([positions[str(g)] for g in frame.target_symbol])
        sizes = 4.0 + 20.0 * np.sqrt(frame.edge_weight.to_numpy(dtype=float) / vmax)
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=sizes,
            c=COLORS[direction],
            edgecolors="white",
            linewidths=0.45,
            alpha=0.95,
            zorder=3,
        )

    ax.scatter([0], [0], s=105, c=COLORS["tf"], edgecolors="white", linewidths=0.8, zorder=6)
    ax.text(0, 0, tf, ha="center", va="center", fontsize=6.5, fontweight="bold", fontstyle="italic", zorder=7)

    labels = q.sort_values("edge_weight", ascending=False).head(min(12, len(q))).copy()
    labels["x"] = labels.target_symbol.astype(str).map(lambda gene: positions[gene][0])
    labels["y"] = labels.target_symbol.astype(str).map(lambda gene: positions[gene][1])
    labels["upper"] = labels.direction.eq("positive")
    labels["right"] = labels.x.ge(0)
    for (upper, right), block in labels.groupby(["upper", "right"], sort=False):
        block = block.sort_values("y", ascending=False)
        slots = np.linspace(2.45, 0.35, len(block)) if upper else np.linspace(-0.35, -2.45, len(block))
        text_x = 3.20 if right else -3.20
        align = "left" if right else "right"
        for slot, row in zip(slots, block.itertuples(index=False)):
            x, y = positions[str(row.target_symbol)]
            end_x = text_x - 0.08 if right else text_x + 0.08
            ax.plot([x, end_x], [y, slot], color="#777777", lw=0.3, zorder=7)
            ax.text(text_x, slot, str(row.target_symbol), fontsize=6, fontstyle="italic", ha=align, va="center", zorder=8)

    ax.set_xlim(-3.65, 3.65)
    ax.set_ylim(-3.35, 3.35)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.text(0.5, 0.965, subtype, ha="center", va="top", fontsize=8.0, fontweight="bold")
    fig.text(0.5, 0.915, f"{tf} regulon network", ha="center", va="top", fontsize=7.0)
    handles = [
        Line2D([0], [0], marker="o", ls="", markerfacecolor=COLORS["positive"], markeredgecolor="white", label="Positive association"),
        Line2D([0], [0], marker="o", ls="", markerfacecolor=COLORS["negative"], markeredgecolor="white", label="Negative association"),
        Line2D([0, 1], [0, 0], color=COLORS["edge_positive"], lw=0.8, label="Positive edge"),
        Line2D([0, 1], [0, 0], color=COLORS["edge_negative"], lw=0.8, ls=(0, (2, 1.5)), label="Negative edge"),
    ]
    ax.legend(
        handles=handles,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.25),
        frameon=False,
        handlelength=1.05,
        columnspacing=0.6,
        handletextpad=0.25,
        borderaxespad=0,
        labelspacing=0.2,
    )
    return fig


def concise_term(term: str, width: int = 18) -> str:
    clean = re.sub(r"\s*\(GO:\d+\)\s*$", "", str(term))
    return "\n".join(textwrap.wrap(clean, width=width))


def go_figure(enrichment: pd.DataFrame, subtype: str, tf: str) -> tuple[plt.Figure, pd.DataFrame]:
    fig = plt.figure(figsize=SQUARE)
    fig.text(0.5, 0.965, subtype, ha="center", va="top", fontsize=8.0, fontweight="bold")
    fig.text(0.5, 0.915, f"{tf} target GO BP", ha="center", va="top", fontsize=7.0)
    ax = fig.add_axes([0.56, 0.19, 0.39, 0.60])
    if enrichment.empty:
        ax.text(0.5, 0.5, "No GO BP result\navailable", ha="center", va="center", fontsize=7)
        ax.axis("off")
        return fig, enrichment

    selected = enrichment.sort_values(["Adjusted P-value", "P-value", "Term"]).head(4).iloc[::-1].copy()
    adjusted = selected["Adjusted P-value"].clip(lower=1e-300)
    x = -np.log10(adjusted)
    counts = selected.Overlap.astype(str).str.split("/").str[0].astype(float).to_numpy()
    y = np.arange(len(selected))
    ax.barh(y, x, color=COLORS["go"], edgecolor=COLORS["go_edge"], lw=0.55, height=0.68)
    ax.scatter(x, y, s=9 + 1.4 * counts, c=COLORS["go_dot"], edgecolor="white", lw=0.35, zorder=3)
    fdr_x = -np.log10(0.05)
    ax.axvline(fdr_x, color="#777777", lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels([concise_term(term) for term in selected.Term])
    ax.set_xlabel("-log10(adjusted P)", labelpad=2)
    ax.spines[["top", "right"]].set_visible(False)
    xmax = max(float(np.nanmax(x)) * 1.10, fdr_x * 1.08, 0.5)
    ax.set_xlim(0, xmax)
    fig.text(0.56, 0.81, "Dashed line: FDR=0.05", ha="left", va="bottom", fontsize=6, color="#666666")
    return fig, selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    configure()

    table_dir = args.project / "05_tables" / "tf_network_enrichment"
    figure_root = args.project / "04_figures" / "tf_network_enrichment_reference"
    edges = pd.read_csv(table_dir / "Top5_TF_all_signed_target_genelist.tsv.gz", sep="\t")
    enrichment = pd.read_csv(table_dir / "Top5_TF_GO_BP_enrichment_expression_background_all.tsv.gz", sep="\t")
    instances = (
        edges[["family", "display_level", "rank", "tf_symbol", "RSS"]]
        .drop_duplicates()
        .sort_values(["family", "display_level", "rank"])
    )

    records: list[dict[str, object]] = []
    displayed_terms: list[pd.DataFrame] = []
    for (family, subtype), block in instances.groupby(["family", "display_level"], sort=False):
        out = figure_root / str(family) / safe(str(subtype))
        out.mkdir(parents=True, exist_ok=True)
        network_bundle = out / f"{safe(str(subtype))}_Top5_TF_signed_networks_reference_v26.pdf"
        go_bundle = out / f"{safe(str(subtype))}_Top5_TF_GO_BP_barplots_reference_v26.pdf"
        with PdfPages(network_bundle) as network_pdf, PdfPages(go_bundle) as go_pdf:
            for row in block.sort_values("rank").itertuples(index=False):
                tf = str(row.tf_symbol)
                edge_block = edges[
                    edges.family.eq(family)
                    & edges.display_level.eq(subtype)
                    & edges.tf_symbol.eq(tf)
                ].copy()
                go_block = enrichment[enrichment.tf_symbol.eq(tf)].copy()
                stem = f"{safe(str(subtype))}_{safe(tf)}"

                network = network_figure(edge_block, str(subtype), tf)
                save(network, out / f"{stem}_signed_network_reference_v26")
                network_pdf.savefig(network)
                plt.close(network)

                go, shown = go_figure(go_block, str(subtype), tf)
                save(go, out / f"{stem}_GO_BP_barplot_reference_v26")
                go_pdf.savefig(go)
                plt.close(go)
                if not shown.empty:
                    shown = shown.copy()
                    shown.insert(0, "family", family)
                    shown.insert(1, "subtype", subtype)
                    shown.insert(2, "rank", int(row.rank))
                    displayed_terms.append(shown)

                graph = nx.DiGraph()
                graph.add_node(tf, node_type="TF")
                for edge in edge_block.itertuples(index=False):
                    graph.add_edge(
                        tf,
                        str(edge.target_symbol),
                        expression_association=str(edge.direction),
                        edge_weight=float(edge.edge_weight),
                        expression_correlation=float(edge.expression_correlation),
                    )
                nx.write_graphml(graph, out / f"{stem}_signed_network_reference_v26.graphml")
                edge_block.to_csv(out / f"{stem}_all_targets.tsv", sep="\t", index=False)
                records.append(
                    {
                        "family": family,
                        "subtype": subtype,
                        "rank": int(row.rank),
                        "tf_symbol": tf,
                        "n_targets": int(edge_block.target_symbol.nunique()),
                        "positive_targets": int((edge_block.direction == "positive").sum()),
                        "negative_targets": int((edge_block.direction == "negative").sum()),
                        "n_go_terms": int(len(go_block)),
                        "status": "PASS" if len(edge_block) else "NOT_ESTIMABLE",
                    }
                )

    pd.DataFrame(records).to_csv(table_dir / "Top5_TF_network_GO_instance_audit.tsv", sep="\t", index=False)
    if displayed_terms:
        pd.concat(displayed_terms, ignore_index=True).to_csv(table_dir / "Top5_TF_GO_BP_displayed_terms.tsv", sep="\t", index=False)
    qa = {
        "status": "PASS" if len(records) == 225 and all(row["status"] == "PASS" for row in records) else "FAIL",
        "instances": len(records),
        "expected_instances": 225,
        "subtypes": int(instances.display_level.nunique()),
        "network_individual_pdf": len(list(figure_root.rglob("*_signed_network_reference_v26.pdf"))),
        "go_individual_pdf": len(list(figure_root.rglob("*_GO_BP_barplot_reference_v26.pdf"))),
        "network_subtype_bundles": len(list(figure_root.rglob("*_Top5_TF_signed_networks_reference_v26.pdf"))),
        "go_subtype_bundles": len(list(figure_root.rglob("*_Top5_TF_GO_BP_barplots_reference_v26.pdf"))),
        "canvas_mm": [52.5, 52.5],
        "font_family": "Arial",
        "font_range_pt": [6, 8.5],
        "edge_sign": "positive/negative expression association, not proven activation/repression",
        "go_statistic": "-log10(BH-adjusted P); complete tables retained",
        "reference_style": "rat_fibk14_grn_v16_render.py visual grammar parameterized for human Top5 TF instances",
    }
    (args.project / "10_qa" / "TF_NETWORK_GO_REFERENCE_V26_AUDIT.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    print(json.dumps(qa, indent=2))


if __name__ == "__main__":
    main()
