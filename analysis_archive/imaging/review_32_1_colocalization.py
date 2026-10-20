from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, maximum_filter
import tifffile
import zarr


ROOT = Path(r"data/external/local_de73c97db013")
OME = ROOT / "Fib_K14_output/01_ome_tiff/32-1 PDGFRα 520+P63 570+KLF4 620+KRT14 690.ome.tiff"
CELLS = ROOT / "IF_colocalized_cell_statistics_all_slides_20260825/03_per_cell/Fib_K14/32-1_per_cell.tsv.gz"
OUT = ROOT / "IF_colocalized_cell_statistics_all_slides_20260825/02_qc/32-1_colocalization_review"
OUT.mkdir(parents=True, exist_ok=True)

PIXEL_UM_L2 = 1.0113780498504639
PDGFRA_THRESHOLD = 25.6291640625
KRT14_THRESHOLD = 18.312668554687498
CROP_SIZE_L2 = 900


def norm(x: np.ndarray, gamma: float = 0.72) -> np.ndarray:
    x = x.astype(np.float32)
    nz = x[x > 0]
    if nz.size == 0:
        return np.zeros_like(x, dtype=np.float32)
    lo, hi = np.percentile(nz, [1.0, 99.7])
    if hi <= lo:
        return np.zeros_like(x, dtype=np.float32)
    y = np.clip((x - lo) / (hi - lo), 0, 1)
    return y ** gamma


def rgb(block: np.ndarray) -> np.ndarray:
    dapi = norm(block[0], 0.70)
    pdgfra = norm(block[1], 0.78)
    krt14 = norm(block[2], 0.78)
    out = np.zeros((*dapi.shape, 3), np.float32)
    out[..., 0] = krt14
    out[..., 1] = pdgfra
    out[..., 2] = np.maximum(dapi, 0.20 * krt14)
    return np.clip(out, 0, 1)


cells = pd.read_csv(CELLS, sep="\t")
co = cells[cells["co_positive"].astype(bool)].copy()

store6 = tifffile.imread(OME, series=0, level=6, aszarr=True)
arr6 = zarr.open(store6, mode="r")
whole = np.asarray(arr6[[0, 1, 4], :, :])
h6, w6 = whole.shape[1:]
store6.close()

scale_x = 19008 / w6
scale_y = 18946 / h6
co["x_l6"] = co["centroid_x_level2_px"] / scale_x
co["y_l6"] = co["centroid_y_level2_px"] / scale_y

density = np.zeros((h6, w6), np.float32)
xi = np.clip(np.rint(co["x_l6"]).astype(int), 0, w6 - 1)
yi = np.clip(np.rint(co["y_l6"]).astype(int), 0, h6 - 1)
np.add.at(density, (yi, xi), 1)
density = gaussian_filter(density, sigma=13)
peaks = density == maximum_filter(density, size=121)
candidate_y, candidate_x = np.where(peaks & (density > 0))
order = np.argsort(density[candidate_y, candidate_x])[::-1]
centers = []
for idx in order:
    x, y = float(candidate_x[idx]), float(candidate_y[idx])
    if all((x - xx) ** 2 + (y - yy) ** 2 > 190**2 for xx, yy in centers):
        centers.append((x, y))
    if len(centers) == 3:
        break

store2 = tifffile.imread(OME, series=0, level=2, aszarr=True)
arr2 = zarr.open(store2, mode="r")
h2, w2 = arr2.shape[1:]

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
fig = plt.figure(figsize=(12.0, 7.2), facecolor="black")
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.08, wspace=0.025)
ax0 = fig.add_subplot(gs[0, :])
ax0.imshow(rgb(whole))
ax0.scatter(co["x_l6"], co["y_l6"], s=4, c="#ff3b30", alpha=0.68, linewidths=0)

records = []
letters = ["A", "B", "C"]
crop_half6_x = (CROP_SIZE_L2 / scale_x) / 2
crop_half6_y = (CROP_SIZE_L2 / scale_y) / 2
for letter, (cx6, cy6) in zip(letters, centers):
    ax0.add_patch(Rectangle((cx6 - crop_half6_x, cy6 - crop_half6_y), 2 * crop_half6_x, 2 * crop_half6_y,
                            fill=False, edgecolor="white", linewidth=1.2))
    ax0.text(cx6 - crop_half6_x, cy6 - crop_half6_y - 8, letter, color="white", weight="bold", fontsize=10)
ax0.set_title("32-1 whole slide: algorithm-positive PDGFRA+ KRT14+ DAPI+ cells (red points)", color="white", pad=3)
ax0.axis("off")

for ax, letter, (cx6, cy6) in zip([fig.add_subplot(gs[1, i]) for i in range(3)], letters, centers):
    cx2, cy2 = cx6 * scale_x, cy6 * scale_y
    x0 = int(np.clip(round(cx2 - CROP_SIZE_L2 / 2), 0, w2 - CROP_SIZE_L2))
    y0 = int(np.clip(round(cy2 - CROP_SIZE_L2 / 2), 0, h2 - CROP_SIZE_L2))
    x1, y1 = x0 + CROP_SIZE_L2, y0 + CROP_SIZE_L2
    block = np.asarray(arr2[[0, 1, 4], y0:y1, x0:x1])
    local = co[(co["centroid_x_level2_px"] >= x0) & (co["centroid_x_level2_px"] < x1) &
               (co["centroid_y_level2_px"] >= y0) & (co["centroid_y_level2_px"] < y1)].copy()
    ax.imshow(rgb(block))
    for row in local.itertuples(index=False):
        ax.add_patch(Circle((row.centroid_x_level2_px - x0, row.centroid_y_level2_px - y0),
                            radius=7, fill=False, edgecolor="white", linewidth=0.65))
    bar_um = 200
    bar_px = bar_um / PIXEL_UM_L2
    ax.plot([25, 25 + bar_px], [CROP_SIZE_L2 - 28, CROP_SIZE_L2 - 28], color="white", lw=3)
    ax.text(25, CROP_SIZE_L2 - 42, "200 µm", color="white", fontsize=7, va="bottom")
    ax.set_title(f"{letter}: {len(local)} algorithm-positive cells", color="white", pad=2)
    ax.axis("off")
    local[["cell_id", "centroid_x_level2_px", "centroid_y_level2_px", "area_um2",
           "PDGFRA_p90_bgsub", "KRT14_p90_bgsub"]].assign(region=letter).to_csv(
               OUT / f"region_{letter}_positive_cells.tsv", sep="\t", index=False)
    records.append({"region": letter, "x0_level2_px": x0, "y0_level2_px": y0,
                    "x1_level2_px": x1, "y1_level2_px": y1, "positive_cells_n": len(local)})

store2.close()
fig.text(0.012, 0.012,
         "Display: DAPI blue, PDGFRA green, KRT14 red. White circles are segmented DAPI nuclei whose 5 µm centroid disk exceeded both locked thresholds.",
         color="white", fontsize=7)
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT / f"32-1_colocalization_review.{ext}", dpi=300 if ext == "png" else None,
                facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)

pd.DataFrame(records).to_csv(OUT / "review_regions.tsv", sep="\t", index=False)
co[["cell_id", "centroid_x_level2_px", "centroid_y_level2_px", "area_um2",
    "PDGFRA_p90_bgsub", "KRT14_p90_bgsub"]].to_csv(OUT / "32-1_all_algorithm_positive_cells.tsv", sep="\t", index=False)
(OUT / "method.txt").write_text(
    "DAPI nuclei segmented at pyramid level 2 (1.0113780498504639 um/px); accepted nuclear area 20-300 um2. "
    "After Gaussian background subtraction (sigma 15 um), marker signal was the 90th percentile in a 5 um radius disk around each nuclear centroid. "
    f"PDGFRA positive if > {PDGFRA_THRESHOLD}; KRT14 positive if > {KRT14_THRESHOLD}; co-positive if DAPI valid and both marker thresholds were exceeded. "
    "Thresholds were pooled/locked across the Fib_K14 analysis family. This is an exploratory cell approximation, not membrane-resolved single-cell segmentation.\n",
    encoding="utf-8")
print({"all_nuclei": len(cells), "co_positive_n": len(co), "ratio_pct": 100 * len(co) / len(cells), "centers_l6": centers})
