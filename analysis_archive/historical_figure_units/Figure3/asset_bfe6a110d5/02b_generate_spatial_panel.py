"""Generate Figure S4F without importing Scanpy (legacy h5ad compatibility)."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import gaussian_filter


SAMPLES = [
    ("Normal", "NS_C02847B1", 180),
    ("5 dpb", "BW32_A01597A3_SDSDB_5dpb", 135),
    ("12 dpb DPTDI-1", "BW13_1_B1_DSDB_12dpb", 45),
    ("12 dpb SPTDI-1", "BW14_1_C1_SSDB_12dpb", 135),
    ("12 dpb DPTDI-2", "BW13_A3_DSDB_12dpb", 135),
    ("12 dpb SPTDI-2", "BW14_B3_SSDB_12dpb", 135),
    ("19 dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("19 dpb p1", "BW15D_1_D1_SDSDB_19dpb", 140),
    ("26 dpb p1", "BW81_C02846B6_SDSDB_26dpb_part1", 0),
    ("26 dpb p2", "BW81_C02846B6_SDSDB_26dpb_part2", 0),
    ("2 mph", "2mph_A03699G6.SCT", 45),
]


def rotate(xy: np.ndarray, angle: float) -> np.ndarray:
    theta = np.deg2rad(angle)
    center = np.median(xy, axis=0)
    shifted = xy - center
    c, s = np.cos(theta), np.sin(theta)
    return np.column_stack((
        shifted[:, 0] * c - shifted[:, 1] * s,
        shifted[:, 0] * s + shifted[:, 1] * c,
    ))


def density_grid(xy: np.ndarray, w: np.ndarray, bins: int = 360):
    xmin, ymin = xy.min(0); xmax, ymax = xy.max(0)
    xi = np.clip(((xy[:, 0] - xmin) / max(xmax - xmin, 1e-9) * (bins - 1)).astype(int), 0, bins - 1)
    yi = np.clip(((xy[:, 1] - ymin) / max(ymax - ymin, 1e-9) * (bins - 1)).astype(int), 0, bins - 1)
    a = np.zeros((bins, bins)); n = np.zeros_like(a)
    np.add.at(a, (yi, xi), w); np.add.at(n, (yi, xi), 1)
    a = gaussian_filter(a, 2.8); n = gaussian_filter(n, 2.8)
    z = np.divide(a, n, out=np.zeros_like(a), where=n > 1e-5)
    z[n < np.quantile(n[n > 0], 0.05)] = np.nan
    return z, (xmin, xmax, ymin, ymax)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--package", required=True)
    p.add_argument("--compute-only", action="store_true")
    args = p.parse_args()
    root = Path(args.package)
    figures = root / "figures" / "panels"
    tables = root / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    a = ad.read_h5ad(args.h5ad, backed="r")
    sample_values = a.obs["sample_batch_new"].astype(str)
    spatial_all = np.asarray(a.obsm["spatial"])
    rows = []
    grid_payload = {}
    fig, axes = plt.subplots(3, 4, figsize=(7.08, 4.95))
    axes = axes.ravel()
    for ax, (label, sid, angle) in zip(axes, SAMPLES):
        idx = np.flatnonzero(sample_values.eq(sid).to_numpy())
        if not len(idx):
            raise KeyError(sid)
        xy = rotate(spatial_all[idx], angle)
        x = a[idx, ["PDGFRA", "KRT14"]].X
        x = x.toarray() if sparse.issparse(x) else np.asarray(x)
        pdg, k14 = x[:, 0], x[:, 1]
        pdg = pdg / max(np.quantile(pdg[pdg > 0], .99), 1e-9) if np.any(pdg > 0) else pdg
        k14 = k14 / max(np.quantile(k14[k14 > 0], .99), 1e-9) if np.any(k14 > 0) else k14
        w = np.sqrt(np.clip(pdg, 0, 1) * np.clip(k14, 0, 1))
        z, extent = density_grid(xy, w)
        grid_payload[f"grid_{len(rows)}"] = z
        grid_payload[f"extent_{len(rows)}"] = np.asarray(extent)
        vmax = max(np.nanquantile(z, .995), 1e-6)
        sca = ax.imshow(z, origin="lower", extent=extent, cmap="viridis",
                        interpolation="bilinear", aspect="equal",
                        norm=Normalize(0, vmax), rasterized=True)
        ax.set_title(label, pad=2)
        ax.set_aspect("equal"); ax.axis("off")
        rows.append({
            "section": label, "sample_id": sid, "n_spatial_bins": len(idx),
            "joint_weight_mean": float(np.mean(w)),
            "joint_weight_q95": float(np.quantile(w, .95)),
            "joint_weight_q99": float(np.quantile(w, .99)),
        })
        print(label, len(idx), flush=True)
    for ax in axes[len(SAMPLES):]:
        ax.axis("off")
    np.savez_compressed(tables / "S4F_spatial_density_grids.npz", **grid_payload)
    pd.DataFrame(rows).to_csv(
        tables / "S4F_human_spatial_joint_density_summary.csv", index=False
    )
    if args.compute_only:
        plt.close(fig)
        return
    cax = fig.add_axes([.76, .07, .18, .018])
    cb = fig.colorbar(sca, cax=cax, orientation="horizontal")
    cb.set_label("PDGFRA + KRT14 joint spatial density", fontsize=6)
    cb.set_ticks([])
    fig.subplots_adjust(left=.02, right=.98, top=.96, bottom=.12,
                        wspace=.08, hspace=.15)
    stem = figures / "FigureS4_F_all_human_spatial_joint_density"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
