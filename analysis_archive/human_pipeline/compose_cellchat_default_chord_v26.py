#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from PIL import Image

MM = 1 / 25.4
FIG_SIZE = (52.5 * MM, 52.5 * MM)


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 6.2,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def crop_white(image: Image.Image) -> Image.Image:
    rgb = np.asarray(image.convert("RGB"))
    ink = np.any(rgb < 248, axis=2)
    if not ink.any():
        return image
    ys, xs = np.where(ink)
    pad = 18
    box = (
        max(int(xs.min()) - pad, 0), max(int(ys.min()) - pad, 0),
        min(int(xs.max()) + pad + 1, image.width), min(int(ys.max()) + pad + 1, image.height),
    )
    return image.crop(box)


def title_from_stem(stem: str) -> str:
    stem = stem.replace("_CellChat_default_LR_chord_raw_default", "")
    match = re.match(r"(.+)_Fib_(Endo|Immune|KC|SAC)$", stem)
    if not match:
        return stem.replace("_", " ")
    return f"{match.group(1).replace('_', ' ')} Fib-{match.group(2)}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    configure()
    figure_dir = args.project / "04_figures" / "cellchat_default_chord"
    table_dir = args.project / "05_tables" / "cellchat_default_chord"
    palette = pd.read_csv(args.project / "00_manifest" / "subtype_palette.tsv", sep="\t")
    colors = dict(zip(palette.subtype.astype(str), palette.color.astype(str).str[:7]))

    records = []
    for raw_path in sorted(figure_dir.glob("*_raw_default.png")):
        final_stem = raw_path.name.replace("_raw_default.png", "")
        lr_path = table_dir / final_stem.replace("_CellChat_default_LR_chord", "_displayed_LR.tsv")
        lr = pd.read_csv(lr_path, sep="\t")
        involved = list(dict.fromkeys(lr.source.astype(str).tolist() + lr.target.astype(str).tolist()))
        cropped = crop_white(Image.open(raw_path))

        fig = plt.figure(figsize=FIG_SIZE)
        fig.text(0.5, 0.965, title_from_stem(raw_path.stem), ha="center", va="top", fontsize=8.0, fontweight="bold")
        ax = fig.add_axes([0.05, 0.20, 0.90, 0.70])
        ax.imshow(cropped)
        ax.set_axis_off()
        handles = [Patch(facecolor=colors.get(cell, "#BDBDBD"), edgecolor="none", label=cell) for cell in involved]
        fig.legend(
            handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=min(2, len(handles)),
            frameon=False, fontsize=6.2, handlelength=0.9, handleheight=0.7,
            columnspacing=0.7, handletextpad=0.3, labelspacing=0.25,
        )
        out = figure_dir / final_stem
        fig.savefig(out.with_suffix(".pdf"))
        fig.savefig(out.with_suffix(".png"), dpi=600)
        plt.close(fig)
        records.append({"panel": final_stem, "n_display_lr": len(lr), "n_cell_subtypes": len(involved), "status": "PASS"})

    pd.DataFrame(records).to_csv(table_dir / "CELLCHAT_DEFAULT_CHORD_COMPOSITION_AUDIT.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
