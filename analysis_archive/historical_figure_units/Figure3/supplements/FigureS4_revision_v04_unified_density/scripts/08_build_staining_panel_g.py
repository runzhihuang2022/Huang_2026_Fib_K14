"""Build Figure S4G from the registered VIM/KRT14/DAPI whole-slide scans.

The statistical unit is one slide/sample. Three non-overlapping fields are
sampled per slide for image measurement and averaged before time-course
visualization. Nuclei are detected from DAPI. VIM and KRT14 positivity
thresholds are derived once from the pooled per-nucleus signal distributions
and then applied unchanged to every slide.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openslide
import pandas as pd
from PIL import Image
from skimage import measure, morphology
from skimage.filters import threshold_otsu


ROOT = Path(
    r"data/external/local_3413b1de7846"
    r"\修回最终版文件\主图和附件\张伟染色结果\20260723染色数据给师兄"
    r"\SH  20230915  张伟 VIM+krt10 krt14 X30  DQ3046\VIM KRT14X15"
)
SAMPLES = {
    "BS1": {"file": "BS1 VIM KRT14.mrxs", "day": 11, "group": "DPTDI"},
    "WB14-1": {"file": "WB14-1 VIM KRT14.mrxs", "day": 12, "group": "DPTDI"},
    "BW4": {"file": "BW4 VIM KRT14.mrxs", "day": 12, "group": "DPTDI"},
    "BW7": {"file": "BW7 VIM KRT14.mrxs", "day": 12, "group": "DPTDI"},
    "BW8": {"file": "BW8 VIM KRT14.mrxs", "day": 12, "group": "DPTDI"},
    "W14-1": {"file": "W14-1 VIM KRT14.mrxs", "day": 14, "group": "SPTDI"},
    "BW3": {"file": "BW3 VIM KRT14.mrxs", "day": 19, "group": "SPTDI"},
    "BW9": {"file": "BW9 VIM KRT14.mrxs", "day": 26, "group": "DPTDI"},
}
GROUP_COLORS = {"SPTDI": "#E6AB02", "DPTDI": "#D62728"}
LEVEL = 2
FIELD_SIZE = 1024


def setup() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def rgba_on_black(image: Image.Image) -> np.ndarray:
    rgba = image.convert("RGBA")
    background = Image.new("RGBA", rgba.size, (0, 0, 0, 255))
    background.alpha_composite(rgba)
    return np.asarray(background.convert("RGB"))


def field_centers(slide: openslide.OpenSlide, n_fields: int = 3) -> list[tuple[int, int]]:
    mask_level = min(6, slide.level_count - 1)
    mask_image = rgba_on_black(
        slide.read_region((0, 0), mask_level, slide.level_dimensions[mask_level])
    )
    intensity = np.max(mask_image, axis=2)
    coordinates = np.column_stack(np.nonzero(intensity > max(8, np.percentile(intensity, 70))))
    if len(coordinates) < 30:
        coordinates = np.column_stack(np.nonzero(intensity > 2))
    centered = coordinates - np.mean(coordinates, axis=0)
    _, _, vectors = np.linalg.svd(centered, full_matrices=False)
    projection = centered @ vectors[0]
    centers = []
    downsample = float(slide.level_downsamples[mask_level])
    for quantile in np.linspace(0.2, 0.8, n_fields):
        target = np.quantile(projection, quantile)
        nearest = coordinates[np.argmin(np.abs(projection - target))]
        centers.append((int(nearest[1] * downsample), int(nearest[0] * downsample)))
    return centers


def read_field(slide: openslide.OpenSlide, center: tuple[int, int]) -> np.ndarray:
    downsample = float(slide.level_downsamples[LEVEL])
    half = int(FIELD_SIZE * downsample / 2)
    origin = (max(center[0] - half, 0), max(center[1] - half, 0))
    return rgba_on_black(slide.read_region(origin, LEVEL, (FIELD_SIZE, FIELD_SIZE)))


def nucleus_signals(rgb: np.ndarray) -> pd.DataFrame:
    red = rgb[:, :, 0].astype(float)
    green = rgb[:, :, 1].astype(float)
    blue = rgb[:, :, 2].astype(float)
    dapi = np.clip(blue - 0.25 * red - 0.25 * green, 0, None)
    positive = dapi[dapi > 1]
    if len(positive) < 50:
        return pd.DataFrame(columns=["VIM_signal", "KRT14_signal"])
    threshold = threshold_otsu(positive)
    nuclei = dapi > threshold
    nuclei = morphology.remove_small_objects(nuclei, min_size=10)
    nuclei = morphology.binary_opening(nuclei, morphology.disk(1))
    labels = measure.label(nuclei)
    rows = []
    for region in measure.regionprops(labels):
        if region.area < 10 or region.area > 1800:
            continue
        y, x = (int(round(value)) for value in region.centroid)
        y0, y1 = max(y - 6, 0), min(y + 7, rgb.shape[0])
        x0, x1 = max(x - 6, 0), min(x + 7, rgb.shape[1])
        rows.append(
            {
                "VIM_signal": float(np.percentile(green[y0:y1, x0:x1], 75)),
                "KRT14_signal": float(np.percentile(red[y0:y1, x0:x1], 75)),
            }
        )
    return pd.DataFrame(rows)


def add_scale_bar(ax: plt.Axes, mpp: float, downsample: float) -> None:
    length_pixels = 200 / (mpp * downsample)
    x0, y0 = 55, FIELD_SIZE - 70
    ax.plot([x0, x0 + length_pixels], [y0, y0], color="white", lw=2)
    ax.text(
        x0 + length_pixels / 2,
        y0 - 22,
        "200 µm",
        color="white",
        ha="center",
        va="bottom",
        fontsize=6.5,
    )


def save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", required=True)
    args = parser.parse_args()
    setup()
    package = Path(args.package)
    field_records: list[dict] = []
    nucleus_tables = []
    representative_images = {}
    representative_counts = {}
    calibration = {}
    for sample, metadata in SAMPLES.items():
        slide = openslide.OpenSlide(str(ROOT / metadata["file"]))
        mpp = float(slide.properties[openslide.PROPERTY_NAME_MPP_X])
        calibration[sample] = (mpp, float(slide.level_downsamples[LEVEL]))
        for field_index, center in enumerate(field_centers(slide)):
            rgb = read_field(slide, center)
            signals = nucleus_signals(rgb)
            if len(signals) > representative_counts.get(sample, -1):
                representative_images[sample] = rgb
                representative_counts[sample] = len(signals)
            signals["sample"] = sample
            signals["field"] = field_index + 1
            nucleus_tables.append(signals)
            field_records.append(
                {
                    "sample": sample,
                    "field": field_index + 1,
                    "day": metadata["day"],
                    "group": metadata["group"],
                    "n_nuclei": len(signals),
                }
            )
        slide.close()
    nuclei = pd.concat(nucleus_tables, ignore_index=True)
    vim_values = nuclei.loc[nuclei["VIM_signal"] > 0, "VIM_signal"].to_numpy()
    krt_values = nuclei.loc[nuclei["KRT14_signal"] > 0, "KRT14_signal"].to_numpy()
    vim_threshold = float(threshold_otsu(vim_values))
    krt_threshold = float(threshold_otsu(krt_values))
    nuclei["VIM_positive"] = nuclei["VIM_signal"] > vim_threshold
    nuclei["KRT14_positive"] = nuclei["KRT14_signal"] > krt_threshold
    nuclei["double_positive"] = nuclei["VIM_positive"] & nuclei["KRT14_positive"]
    field_stats = (
        nuclei.groupby(["sample", "field"], observed=True)["double_positive"]
        .agg(["sum", "count", "mean"])
        .reset_index()
        .rename(
            columns={
                "sum": "n_double_positive",
                "count": "n_nuclei",
                "mean": "double_positive_fraction",
            }
        )
    )
    field_meta = pd.DataFrame(field_records).drop(columns=["n_nuclei"])
    field_stats = field_stats.merge(field_meta, on=["sample", "field"], how="left")
    field_stats = field_stats.loc[field_stats["n_nuclei"] >= 50].copy()
    slide_stats = (
        field_stats.groupby(["sample", "day", "group"], observed=True)
        .agg(
            double_positive_fraction=("double_positive_fraction", "mean"),
            field_sem=("double_positive_fraction", "sem"),
            n_fields=("field", "nunique"),
            n_nuclei=("n_nuclei", "sum"),
            n_double_positive=("n_double_positive", "sum"),
        )
        .reset_index()
    )
    nuclei.to_csv(package / "tables" / "S4G_staining_per_nucleus_signals_v01.csv.gz", index=False)
    field_stats.to_csv(package / "tables" / "S4G_staining_field_statistics_v01.csv", index=False)
    slide_stats.to_csv(package / "tables" / "S4G_staining_slide_statistics_v01.csv", index=False)
    pd.DataFrame(
        [
            {
                "VIM_signal_threshold": vim_threshold,
                "KRT14_signal_threshold": krt_threshold,
                "method": "Pooled Otsu thresholds on per-nucleus 75th-percentile channel signals",
            }
        ]
    ).to_csv(package / "tables" / "S4G_staining_thresholds_v01.csv", index=False)

    representatives = ["BS1", "W14-1", "BW3", "BW9"]
    fig, axes = plt.subplots(1, 5, figsize=(7.08, 1.82), gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.45]})
    for ax, sample in zip(axes[:4], representatives):
        ax.imshow(representative_images[sample])
        mpp, downsample = calibration[sample]
        add_scale_bar(ax, mpp, downsample)
        metadata = SAMPLES[sample]
        ax.set_title(f"{metadata['day']} dpb {metadata['group']}", pad=2)
        ax.axis("off")
    stats_ax = axes[4]
    day_order = [11, 12, 14, 19, 26]
    day_positions = {day: index for index, day in enumerate(day_order)}
    for group in ["SPTDI", "DPTDI"]:
        current = slide_stats[slide_stats["group"].eq(group)].sort_values("day")
        stats_ax.errorbar(
            current["day"].map(day_positions),
            current["double_positive_fraction"] * 100,
            yerr=current["field_sem"].fillna(0) * 100,
            color=GROUP_COLORS[group],
            marker="o",
            lw=1.2,
            ms=3.5,
            capsize=2,
            label=group,
        )
    stats_ax.set(
        xlabel="Healing time (dpb)",
        ylabel="",
        xticks=np.arange(len(day_order)),
        xticklabels=day_order,
    )
    stats_ax.set_title("VIM+KRT14 double-positive nuclei (%)", fontsize=7, pad=2)
    stats_ax.legend(frameon=False, loc="upper left")
    stats_ax.grid(axis="y", color="#E6E6E6", lw=0.5)
    stats_ax.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.02,
        0.99,
        "VIM",
        color="#00A651",
        ha="left",
        va="top",
        fontsize=7,
        fontweight="bold",
    )
    fig.text(0.067, 0.99, "+", color="black", ha="left", va="top", fontsize=7)
    fig.text(
        0.083,
        0.99,
        "KRT14",
        color="#E64B35",
        ha="left",
        va="top",
        fontsize=7,
        fontweight="bold",
    )
    fig.text(0.145, 0.99, "+", color="black", ha="left", va="top", fontsize=7)
    fig.text(
        0.16,
        0.99,
        "DAPI",
        color="#2F55D4",
        ha="left",
        va="top",
        fontsize=7,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.015, right=0.995, top=0.88, bottom=0.22, wspace=0.22)
    save(fig, package / "figures" / "panels" / "FigureS4_G_VIM_KRT14_staining_v01")
    print(package / "figures" / "panels" / "FigureS4_G_VIM_KRT14_staining_v01.pdf")


if __name__ == "__main__":
    main()
