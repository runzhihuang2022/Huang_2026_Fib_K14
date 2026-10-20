#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

MM_IN = 1 / 25.4
MODES = [
    ("Mode1", "D0_7", "Mode 1", "0–7 dpb"),
    ("Mode2", "D15_28", "Mode 2", "15–28 dpb"),
    ("Mode3", "D8_14", "Mode 3", "8–14 dpb"),
    ("Mode4", "M1_2", "Mode 4", "1–2 mph"),
]
INTERFACES = ["SAC", "KC", "Endo", "Immune"]
DISPLAY = {"SAC": "Fib–SAC", "KC": "Fib–KC", "Endo": "Fib–Endo", "Immune": "Fib–Immune"}
COLORS = {"SAC": "#1B9E77", "KC": "#D95F02", "Endo": "#2C7FB8", "Immune": "#7A5195"}
FOCUS = {
    "SAC": ["Positive Regulation Of Epithelial Cell Proliferation", "Regulation Of Cell Migration",
            "Epidermis Development", "Epithelial Cell Proliferation"],
    "KC": ["Keratinocyte Differentiation", "Epidermal Cell Differentiation",
           "Positive Regulation Of Cell Migration", "Epidermis Development"],
    "Endo": ["Positive Regulation Of Angiogenesis", "Angiogenesis", "Blood Vessel Morphogenesis",
             "Endothelial Cell Proliferation"],
    "Immune": ["Positive Regulation Of Immune Response", "Macrophage Activation",
               "Regulation Of Inflammatory Response", "Leukocyte Migration"],
}


def split_complex(values: pd.Series) -> set[str]:
    genes: set[str] = set()
    for value in values.dropna().astype(str):
        genes.update(x for x in value.split("_") if x)
    return genes


def clean_term(term: str) -> str:
    return re.sub(r"\s*\(GO:\d+\)$", "", str(term)).strip()


def stars(q: float) -> str:
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def enrich(genes: set[str], background: set[str], go: dict[str, list[str]]) -> pd.DataFrame:
    genes &= background
    rows = []
    for term, raw_members in go.items():
        members = set(raw_members) & background
        overlap = genes & members
        if overlap:
            rows.append({
                "term": term, "lr_gene_n": len(genes), "background_n": len(background),
                "term_n": len(members), "overlap_n": len(overlap),
                "overlap_genes": ";".join(sorted(overlap)),
                "p_value": hypergeom.sf(len(overlap) - 1, len(background), len(members), len(genes)),
            })
    columns = ["term", "lr_gene_n", "background_n", "term_n", "overlap_n", "overlap_genes", "p_value", "fdr_bh"]
    if not rows:
        return pd.DataFrame(columns=columns)
    out = pd.DataFrame(rows)
    p = out.p_value.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.minimum(ranked, 1.0)
    out["fdr_bh"] = q
    return out.sort_values(["fdr_bh", "p_value", "term"]).reset_index(drop=True)


def select_focus(table: pd.DataFrame, interface: str, n: int = 6) -> pd.DataFrame:
    table = table[table.fdr_bh < 0.05].copy()
    selected, seen = [], set()
    for wanted in FOCUS[interface]:
        cleaned = table.term.map(clean_term).str.casefold()
        hit = table[cleaned == wanted.casefold()]
        if hit.empty:
            hit = table[cleaned.str.startswith(wanted.casefold())]
        if not hit.empty:
            row = hit.iloc[0].copy()
            row["selection_basis"] = "pre-specified repair function"
            selected.append(row); seen.add(str(row.term))
    for _, source in table.iterrows():
        if len(selected) >= n:
            break
        if str(source.term) not in seen:
            row = source.copy(); row["selection_basis"] = "top BH-adjusted result used to fill display"
            selected.append(row); seen.add(str(row.term))
    return pd.DataFrame(selected[:n])


def render(table: pd.DataFrame, mode_label: str, time_label: str, interface: str, out: Path,
           reason: str = "", estimable: bool = True) -> None:
    fig, ax = plt.subplots(figsize=(42 * MM_IN, 84 * MM_IN))
    shown = select_focus(table, interface) if len(table) else table.copy()
    if not estimable:
        fig.text(0.5, 0.57, "NOT_ESTIMABLE", ha="center", va="center",
                 fontsize=9, fontweight="bold", color="#A61E1E")
        if "M1_Macrophage n=5" in reason:
            detail = "M1_Macrophage n=5\nM2_Macrophage n=1"
        elif "M1_Macrophage n=2" in reason:
            detail = "M1_Macrophage n=2\nM2_Macrophage n=0"
        else:
            detail = "Required subtype count <10"
        fig.text(0.5, 0.45, detail, ha="center", va="center", fontsize=7, color="#4A4A4A")
        ax.axis("off")
    elif shown.empty:
        fig.text(0.5, 0.54, "NO SIGNIFICANT\nGO-BP", ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color="#555555")
        fig.text(0.5, 0.41, "BH-FDR q ≥ 0.05", ha="center", va="center", fontsize=7, color="#555555")
        ax.axis("off")
    else:
        shown = shown.iloc[::-1].copy()
        score = -np.log10(shown.fdr_bh.clip(lower=1e-300))
        y = np.arange(len(shown))
        ax.barh(y, score, color=COLORS[interface], edgecolor="#333333", linewidth=0.35, height=0.68)
        labels = ["\n".join(textwrap.wrap(clean_term(x), width=15)[:4]) for x in shown.term]
        ax.set_yticks(y, labels)
        ax.set_xlabel("–log10(q)", fontsize=8)
        xmax = max(float(score.max()), 1.0)
        for yi, value, qval in zip(y, score, shown.fdr_bh):
            ax.text(float(value) + xmax * 0.025, yi, stars(float(qval)), fontsize=9,
                    fontweight="bold", ha="left", va="center")
        ax.axvline(-np.log10(0.05), color="#777777", lw=0.55, ls="--")
        ax.set_xlim(0, xmax * 1.24)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=7)
    fig.suptitle(f"{DISPLAY[interface]}\n{mode_label} | {time_label}", y=0.98, fontsize=9, fontweight="bold")
    fig.subplots_adjust(left=0.62, right=0.95, bottom=0.10, top=0.88)
    for ext, kwargs in [("pdf", {}), ("svg", {}), ("png", {"dpi": 600}),
                        ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
        fig.savefig(out.with_suffix(f".{ext}"), **kwargs)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.titlesize": 9,
                         "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    table_dir = root / "02_analysis" / "required_chord" / "tables_extended5"
    out_dir = root / "05_enrichment" / "D_GO_extended5"
    out_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(table_dir / "REQUIRED_CHORD_INSTANCE_AUDIT.tsv", sep="\t")
    go = json.loads((root / "01_input_adapter" / "GO_Biological_Process_2023.json").read_text(encoding="utf-8-sig"))
    background = set(pd.read_csv(root / "02_analysis" / "targeted_cellchat" / "tables" /
                                 "CELLCHAT_SIGNALING_BACKGROUND_GENES.tsv", sep="\t").iloc[:, 0].dropna().astype(str))
    audit_rows, positive = [], []
    for mode, group, mode_label, time_label in MODES:
        for interface in INTERFACES:
            row = audit[(audit.group == group) & (audit.interface == f"Fib_{interface}")].iloc[0]
            source = table_dir / f"{group}_Fib_{interface}_all_significant_LR.tsv"
            lr = pd.read_csv(source, sep="\t") if source.exists() else pd.DataFrame()
            genes = ((split_complex(lr.ligand) | split_complex(lr.receptor)) & background) if len(lr) else set()
            result = enrich(genes, background, go) if row.status == "PASS" else enrich(set(), background, go)
            result.insert(0, "interface", DISPLAY[interface]); result.insert(0, "time_group", time_label); result.insert(0, "mode", mode_label)
            result.to_csv(out_dir / f"{mode}_{interface}_required_cell_GO_all.tsv", sep="\t", index=False)
            focus = select_focus(result, interface) if len(result) else result.copy()
            if len(focus):
                focus.to_csv(out_dir / f"{mode}_{interface}_focused_GO_selection.tsv", sep="\t", index=False)
            pos = result[result.fdr_bh < 0.05].copy() if len(result) else result.copy()
            if len(pos):
                pos["interface_short"] = interface; positive.append(pos)
            render(result, mode_label, time_label, interface,
                   out_dir / f"{mode}_{interface}_required_cell_GO_42x84mm", str(row.reason), row.status == "PASS")
            audit_rows.append({
                "mode": mode_label, "group": group, "time_group": time_label,
                "interface": f"Fib-{interface}", "status": row.status,
                "significant_lr_n": len(lr), "lr_gene_n": len(genes),
                "go_terms_tested": len(result),
                "go_fdr_lt_0.05": int((result.fdr_bh < 0.05).sum()) if len(result) else 0,
                "gene_set_basis": "all significant LR pairs among the five prespecified subtypes",
            })
    pd.DataFrame(audit_rows).to_csv(out_dir / "MODE_REQUIRED_GO_AUDIT.tsv", sep="\t", index=False)
    positives = pd.concat(positive, ignore_index=True).sort_values(["mode", "interface", "fdr_bh"]) if positive else pd.DataFrame()
    positives.to_csv(out_dir / "MODE_REQUIRED_GO_POSITIVE_RESULTS.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
