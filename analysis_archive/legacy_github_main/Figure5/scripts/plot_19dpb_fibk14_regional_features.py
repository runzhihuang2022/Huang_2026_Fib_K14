#!/usr/bin/env python
"""White-background 19-dpb Stereo-seq feature maps and Fib_K14 regional statistics.

The expression matrix is read from the registered combined H5AD. Wound axes and
cell-state annotations are joined by observation name from the audited metadata
table, while wound territories are reassigned from the manual three-label mask
for each section (1=Un-epi, 2=Epi-Front, 3=Newly-epi). The spatial-anchor JSON
and anatomical boundary mask are required alignment provenance. The script
never falls back to fixed wound-x thresholds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

try:
    import anndata as ad
except ModuleNotFoundError:  # plotting-only runtime does not need AnnData
    ad = None
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
from scipy import sparse, stats
from PIL import Image


SAMPLES = ["19dpb", "19dpb_p1"]
SAMPLE_IDS = {
    "19dpb": "BW15D_C6_SDSDB_19dpb",
    "19dpb_p1": "BW15D_1_D1_SDSDB_19dpb",
}
REGIONS = ["Un-epi", "Epi-Front", "Newly-epi"]
REGION_COLORS = {
    "Un-epi": "#D95F5F",
    "Epi-Front": "#D8B53F",
    "Newly-epi": "#4C78A8",
}
REGION_BACKGROUNDS = {
    "Un-epi": "#FBE9E9",
    "Epi-Front": "#FFF6CF",
    "Newly-epi": "#E8F1FA",
}
FEATURE_GROUPS = {
    "iPSC Reprogramming": ["KLF4", "SOX2", "POU5F1", "MYC"],
    "Epithelialization": ["GRHL3", "TP63", "TACSTD2", "IGFL1"],
    "Scar Formation": ["COL1A1", "ACTA2", "FN1", "TGFB1"],
}
GENES = [g for genes in FEATURE_GROUPS.values() for g in genes]
FEATURES = [x for group, genes in FEATURE_GROUPS.items() for x in [*genes, group]]
PAIR_ORDER = list(combinations(REGIONS, 2))
ZONE_LABELS = {1: "Un-epi", 2: "Epi-Front", 3: "Newly-epi"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_zone_mask(mask_dir: Path, sample: str) -> Path:
    sample_id = SAMPLE_IDS[sample]
    candidates = [
        mask_dir / f"wound_zones_mask_{sample_id}.png",
        mask_dir / f"wound_zones_mask_{sample}.png",
        mask_dir / f"{sample}_wound_zones_mask.png",
    ]
    # The original interactive annotator saved the first section with a generic
    # name in pdf_output. This fallback is intentionally restricted to 19dpb.
    if sample == "19dpb":
        candidates.append(mask_dir.parent / "wound_zones_mask.png")
    for path in candidates:
        if path.exists():
            return path
    expected = "\n  - ".join(str(x) for x in candidates)
    raise FileNotFoundError(
        f"Missing three-label wound-zone mask for {sample} ({sample_id}). "
        f"No coordinate-threshold fallback is allowed. Checked:\n  - {expected}"
    )


def assign_mask_regions(joined: pd.DataFrame, mask_dir: Path) -> tuple[pd.DataFrame, dict]:
    joined = joined.copy()
    joined["wound_territory"] = pd.Series(index=joined.index, dtype="object")
    source_audit: dict[str, dict] = {}
    for sample in SAMPLES:
        sample_id = SAMPLE_IDS[sample]
        anchor_path = mask_dir / f"{sample}_spatial_anchors.json"
        boundary_path = mask_dir / f"boundary_mask_{sample_id}.png"
        zone_path = resolve_zone_mask(mask_dir, sample)
        missing = [str(p) for p in (anchor_path, boundary_path) if not p.exists()]
        if missing:
            raise FileNotFoundError(f"Missing alignment source(s) for {sample}: {missing}")

        anchor = json.loads(anchor_path.read_text(encoding="utf-8"))
        if len(anchor.get("leading_edge", [])) < 2:
            raise ValueError(f"Invalid leading_edge in {anchor_path}")

        zone_img = np.asarray(Image.open(zone_path))
        if zone_img.ndim == 3:
            zone_img = zone_img[:, :, 0]
        labels = set(np.unique(zone_img).astype(int).tolist())
        if not set(ZONE_LABELS).issubset(labels):
            raise ValueError(
                f"{zone_path} is not a complete three-zone mask; labels={sorted(labels)}, "
                f"required={sorted(ZONE_LABELS)}"
            )

        take = joined["sample_key"].astype(str).eq(sample)
        sdf = joined.loc[take]
        rx = sdf["rot_x"].to_numpy(float)
        ry = sdf["rot_y"].to_numpy(float)
        if np.ptp(rx) == 0 or np.ptp(ry) == 0:
            raise ValueError(f"Degenerate rotated-coordinate range for {sample}")
        height, width = zone_img.shape
        px = np.clip(((rx - rx.min()) / np.ptp(rx) * (width - 1)).astype(int), 0, width - 1)
        py = np.clip(((ry - ry.min()) / np.ptp(ry) * (height - 1)).astype(int), 0, height - 1)
        numeric = zone_img[py, px].astype(int)
        assigned = pd.Series(numeric, index=sdf.index).map(ZONE_LABELS)
        joined.loc[take, "wound_territory"] = assigned

        boundary_img = np.asarray(Image.open(boundary_path))
        if boundary_img.ndim == 3:
            boundary_img = boundary_img[:, :, 0]
        source_audit[sample] = {
            "sample_id": sample_id,
            "anchor_json": str(anchor_path),
            "anchor_sha256": sha256(anchor_path),
            "boundary_mask": str(boundary_path),
            "boundary_mask_sha256": sha256(boundary_path),
            "boundary_labels": sorted(np.unique(boundary_img).astype(int).tolist()),
            "wound_zone_mask": str(zone_path),
            "wound_zone_mask_sha256": sha256(zone_path),
            "wound_zone_labels": ZONE_LABELS,
            "assigned_bin_counts": assigned.value_counts(dropna=False).to_dict(),
            "unassigned_bins": int(assigned.isna().sum()),
        }

    joined["wound_territory"] = pd.Categorical(joined["wound_territory"], REGIONS, ordered=True)
    return joined, source_audit


def bh_adjust(values: list[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if pv.size == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * pv.size / np.arange(1, pv.size + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    out[valid] = restored
    return out


def stars(q: float) -> str:
    if not np.isfinite(q):
        return "NA"
    if q <= 1e-4:
        return "****"
    if q <= 1e-3:
        return "***"
    if q <= 1e-2:
        return "**"
    if q <= 5e-2:
        return "*"
    return "ns"


def italic_gene(name: str) -> str:
    return rf"$\it{{{name}}}$" if name in GENES else name


def load_analysis(h5ad_path: Path, metadata_path: Path, mask_dir: Path) -> tuple[pd.DataFrame, dict]:
    if ad is None:
        raise RuntimeError("AnnData is required for extraction mode")
    metadata = pd.read_csv(metadata_path, sep="\t", compression="infer")
    metadata = metadata.loc[metadata["sample_key"].isin(SAMPLES)].copy()
    metadata = metadata.set_index("obs_name", verify_integrity=True)

    a = ad.read_h5ad(h5ad_path, backed="r")
    sample_key = "sample_batch_new"
    expected_ids = set(SAMPLE_IDS.values())
    observed_ids = set(a.obs[sample_key].astype(str).unique())
    missing_ids = expected_ids - observed_ids
    if missing_ids:
        raise ValueError(f"Registered 19-dpb sample IDs missing from H5AD: {sorted(missing_ids)}")

    var_lookup = {str(v).upper(): str(v) for v in a.var_names}
    missing_genes = [g for g in GENES if g not in var_lookup]
    if missing_genes:
        raise ValueError(f"Genes missing from H5AD: {missing_genes}")

    mask = a.obs[sample_key].astype(str).isin(expected_ids).to_numpy()
    selected_vars = [var_lookup[g] for g in GENES]
    sub = a[mask, selected_vars].to_memory()
    expression = sub.X.toarray() if sparse.issparse(sub.X) else np.asarray(sub.X)
    expression = pd.DataFrame(expression, index=sub.obs_names, columns=GENES)
    joined = metadata.join(expression, how="inner")
    if len(joined) != len(metadata):
        raise ValueError(f"Metadata-expression join lost observations: {len(metadata)} -> {len(joined)}")

    # Deterministic four-gene module scores: mean z-score across the two sections.
    for group, genes in FEATURE_GROUPS.items():
        vals = joined[genes].to_numpy(float)
        means = np.nanmean(vals, axis=0)
        sds = np.nanstd(vals, axis=0, ddof=0)
        sds[sds == 0] = 1.0
        joined[group] = np.nanmean((vals - means) / sds, axis=1)

    joined["sample_key"] = pd.Categorical(joined["sample_key"], SAMPLES, ordered=True)
    joined, region_sources = assign_mask_regions(joined, mask_dir)
    audit = {
        "combined_h5ad": str(h5ad_path),
        "metadata_table": str(metadata_path),
        "region_assignment": "manual three-label mask after sample-specific rotated-coordinate registration; no fixed wound-x thresholds",
        "region_sources": region_sources,
        "sample_ids": SAMPLE_IDS,
        "n_bins": {k: int((joined["sample_key"] == k).sum()) for k in SAMPLES},
        "n_fib_k14_bins": {
            k: int(((joined["sample_key"] == k) & (joined["sub_labels"] == "Fib_K14")).sum())
            for k in SAMPLES
        },
        "feature_groups": FEATURE_GROUPS,
        "module_score": "mean of gene-wise z-scores across both registered 19-dpb sections",
    }
    return joined, audit


def make_spatial_blocks(fib: pd.DataFrame, x_width_um: float = 300.0, depth_width: float = 0.1) -> pd.DataFrame:
    x_floor = np.floor(fib["wound_x_um"].to_numpy(float) / x_width_um).astype(int)
    y_floor = np.floor(fib["depth_axis"].to_numpy(float) / depth_width).astype(int)
    fib = fib.copy()
    fib["block_x"] = x_floor
    fib["block_y"] = y_floor
    keys = ["sample_key", "wound_territory", "block_x", "block_y"]
    blocks = fib.groupby(keys, observed=True)[FEATURES].mean().reset_index()
    counts = fib.groupby(keys, observed=True).size().rename("n_bins").reset_index()
    blocks = blocks.merge(counts, on=keys, validate="one_to_one")
    return blocks.loc[blocks["n_bins"] >= 5].copy()


def calculate_statistics(blocks: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary_rows, overall_rows, pair_rows = [], [], []
    for sample in SAMPLES:
        sdf = blocks.loc[blocks["sample_key"] == sample]
        for feature in FEATURES:
            arrays = {}
            for region in REGIONS:
                vals = sdf.loc[sdf["wound_territory"] == region, feature].dropna().to_numpy()
                arrays[region] = vals
                summary_rows.append({
                    "sample_key": sample,
                    "feature": feature,
                    "region": region,
                    "mean": float(np.mean(vals)) if len(vals) else np.nan,
                    "sem": float(stats.sem(vals)) if len(vals) > 1 else np.nan,
                    "median": float(np.median(vals)) if len(vals) else np.nan,
                    "n_spatial_blocks": int(len(vals)),
                })
            valid_arrays = [v for v in arrays.values() if len(v)]
            if len(valid_arrays) == 3:
                try:
                    h, p = stats.kruskal(*valid_arrays)
                except ValueError as exc:
                    if "identical" not in str(exc):
                        raise
                    h, p = 0.0, 1.0
            else:
                h, p = np.nan, np.nan
            overall_rows.append({"sample_key": sample, "feature": feature, "H": h, "p_raw": p})
            for r1, r2 in PAIR_ORDER:
                v1, v2 = arrays[r1], arrays[r2]
                if len(v1) >= 3 and len(v2) >= 3:
                    u, p = stats.mannwhitneyu(v1, v2, alternative="two-sided")
                else:
                    u, p = np.nan, np.nan
                pair_rows.append({
                    "sample_key": sample,
                    "feature": feature,
                    "region_1": r1,
                    "region_2": r2,
                    "U": u,
                    "p_raw": p,
                    "n_blocks_1": len(v1),
                    "n_blocks_2": len(v2),
                })

    summary = pd.DataFrame(summary_rows)
    overall = pd.DataFrame(overall_rows)
    pairwise = pd.DataFrame(pair_rows)
    # Correction families are prespecified per section: 15 overall tests and 45 pairwise tests.
    overall["q_bh_within_section"] = overall.groupby("sample_key", observed=True)["p_raw"].transform(
        lambda x: bh_adjust(x.tolist())
    )
    pairwise["q_bh_within_section"] = pairwise.groupby("sample_key", observed=True)["p_raw"].transform(
        lambda x: bh_adjust(x.tolist())
    )
    pairwise["significance"] = pairwise["q_bh_within_section"].map(stars)
    return summary, overall, pairwise


def style_spatial_axis(ax: plt.Axes, sample_df: pd.DataFrame, show_y: bool, show_x: bool) -> None:
    xmin, xmax = sample_df["wound_x_um"].min(), sample_df["wound_x_um"].max()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(1.02, -0.02)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("Epidermis–dermis axis" if show_y else "")
    ax.set_xlabel("Distance to leading edge (μm)" if show_x else "")
    if not show_y:
        ax.set_yticklabels([])
    ax.tick_params(length=2.5, width=0.7)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
        spine.set_color("#333333")


def plot_spatial(sample_df: pd.DataFrame, sample: str, output: Path, display_limits: dict) -> None:
    print(f"Plotting spatial panels: {sample}", flush=True)
    cmap = LinearSegmentedColormap.from_list(
        "expression", ["#17324D", "#2D6F8E", "#75B6B2", "#F1D36B", "#C43C2F"]
    )
    fig, axes = plt.subplots(3, 5, figsize=(15.5, 9.8), constrained_layout=True)
    for row, (group, genes) in enumerate(FEATURE_GROUPS.items()):
        print(f"  spatial row {row + 1}/3: {group}", flush=True)
        row_features = [*genes, group]
        for col, feature in enumerate(row_features):
            print(f"    feature {feature}", flush=True)
            ax = axes[row, col]
            values = sample_df[feature].to_numpy(float)
            order = np.argsort(values)
            x = sample_df["wound_x_um"].to_numpy(float)[order]
            y = sample_df["depth_axis"].to_numpy(float)[order]
            val = values[order]
            style_spatial_axis(ax, sample_df, show_y=(col == 0), show_x=(row == 2))
            # Render the actual mask-derived territories as the pale tissue
            # background. This preserves non-vertical and sample-specific zone
            # geometry instead of implying fixed +/-500 um boundaries.
            for region in REGIONS:
                region_mask = sample_df["wound_territory"].astype(str).eq(region).to_numpy()[order]
                ax.scatter(
                    x[region_mask], y[region_mask], s=0.55,
                    c=REGION_BACKGROUNDS[region], alpha=0.72,
                    linewidths=0, rasterized=True, zorder=1,
                )
            # A deterministic background sample preserves tissue geometry while
            # keeping high-resolution export memory bounded.
            bg_step = max(1, len(x) // 12000)
            ax.scatter(x[::bg_step], y[::bg_step], s=0.35, c="#D7D7D7", alpha=0.28, linewidths=0, rasterized=True, zorder=2)
            lo, hi = display_limits[feature]
            if feature in GENES:
                positive = val > 0
            else:
                positive = np.isfinite(val)
            sc = ax.scatter(
                x[positive], y[positive], c=np.clip(val[positive], lo, hi),
                s=0.75, cmap=cmap, norm=Normalize(lo, hi), alpha=0.88,
                linewidths=0, rasterized=True, zorder=3,
            )
            ax.set_title(italic_gene(feature), fontsize=12, fontweight="bold" if feature not in GENES else "normal", pad=5)
            if row == 0:
                trans = ax.get_xaxis_transform()
                label_box = dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.5)
                for region in REGIONS:
                    region_x = sample_df.loc[
                        sample_df["wound_territory"].astype(str).eq(region), "wound_x_um"
                    ]
                    if len(region_x):
                        color = "#8A6A0A" if region == "Epi-Front" else REGION_COLORS[region]
                        ax.text(float(region_x.median()), 0.985, region, color=color,
                                transform=trans, ha="center", va="top", fontsize=10, bbox=label_box)
            if col == 0:
                ax.text(-0.26, 0.5, group, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=12, fontweight="bold")
            if col == 4:
                cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.025)
                cb.ax.tick_params(labelsize=10, length=2)
                cb.set_label("Expression" if feature in GENES else "Module score", fontsize=10)
    fig.suptitle(f"19 dpb spatial programs — {sample}", fontsize=14, fontweight="bold")
    for ext in ["png", "pdf", "svg"]:
        print(f"  saving {sample} spatial {ext}", flush=True)
        kwargs = {"dpi": 200} if ext == "png" else {}
        fig.savefig(output / f"{sample}_spatial_programs_white_3x5.{ext}", facecolor="white", bbox_inches="tight", **kwargs)
    plt.close(fig)


def add_brackets(ax: plt.Axes, pairs: pd.DataFrame, base: float, step: float) -> int:
    position = {r: i for i, r in enumerate(REGIONS)}
    significant = pairs.loc[pairs["q_bh_within_section"] <= 0.05].copy()
    for level, (_, row) in enumerate(significant.iterrows()):
        x1, x2 = position[row["region_1"]], position[row["region_2"]]
        y = base + level * step
        ax.plot([x1, x1, x2, x2], [y, y + step * 0.18, y + step * 0.18, y], color="black", lw=0.8, clip_on=False)
        ax.text((x1 + x2) / 2, y + step * 0.2, row["significance"], ha="center", va="bottom", fontsize=10)
    return len(significant)


def plot_statistics(
    blocks: pd.DataFrame,
    summary: pd.DataFrame,
    overall: pd.DataFrame,
    pairwise: pd.DataFrame,
    sample: str,
    output: Path,
) -> None:
    print(f"Plotting Fib_K14 statistics: {sample}", flush=True)
    fig, axes = plt.subplots(3, 5, figsize=(15.5, 9.8), constrained_layout=True)
    for row, (group, genes) in enumerate(FEATURE_GROUPS.items()):
        for col, feature in enumerate([*genes, group]):
            ax = axes[row, col]
            sm = summary.loc[(summary["sample_key"] == sample) & (summary["feature"] == feature)].set_index("region").reindex(REGIONS)
            means = sm["mean"].to_numpy(float)
            sems = sm["sem"].fillna(0).to_numpy(float)
            ns = sm["n_spatial_blocks"].astype(int).to_numpy()
            x = np.arange(3)
            ax.bar(x, means, yerr=sems, capsize=3, color=[REGION_COLORS[r] for r in REGIONS], edgecolor="black", linewidth=0.7, width=0.68)
            bdf = blocks.loc[(blocks["sample_key"] == sample), ["wound_territory", feature]].dropna()
            rng = np.random.default_rng(20260722 + row * 10 + col)
            for i, region in enumerate(REGIONS):
                vals = bdf.loc[bdf["wound_territory"] == region, feature].to_numpy(float)
                jitter = rng.uniform(-0.13, 0.13, len(vals))
                ax.scatter(np.full(len(vals), i) + jitter, vals, s=8, facecolors="white", edgecolors="#333333", linewidths=0.35, alpha=0.75, zorder=3)
            ax.set_xticks(x, [f"{r}\n(n={n})" for r, n in zip(REGIONS, ns)], rotation=25, ha="right")
            ax.set_title(italic_gene(feature), fontsize=12, fontweight="bold" if feature not in GENES else "normal")
            ax.set_ylabel("Mean expression per spatial block" if col == 0 else "")
            ax.spines[["top", "right"]].set_visible(False)
            psub = pairwise.loc[(pairwise["sample_key"] == sample) & (pairwise["feature"] == feature)].copy()
            all_vals = bdf[feature].to_numpy(float)
            data_min = float(np.nanmin(all_vals)) if len(all_vals) else 0.0
            data_max = float(np.nanmax(all_vals)) if len(all_vals) else 1.0
            span = max(data_max - data_min, abs(data_max) * 0.2, 0.05)
            base = data_max + 0.08 * span
            step = 0.12 * span
            n_brackets = add_brackets(ax, psub, base, step)
            ax.set_ylim(data_min - 0.08 * span, base + max(0.9, n_brackets + 0.5) * step)
            qkw = overall.loc[(overall["sample_key"] == sample) & (overall["feature"] == feature), "q_bh_within_section"].iloc[0]
            ax.text(0.02, 0.98, f"Kruskal–Wallis q={qkw:.2g}", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    fig.suptitle(
        f"Fib_K14 regional expression — {sample}\n300 μm × 0.1-depth spatial-block means; pairwise Mann–Whitney, BH-FDR across 45 tests",
        fontsize=14, fontweight="bold",
    )
    for ext in ["png", "pdf", "svg"]:
        kwargs = {"dpi": 200} if ext == "png" else {}
        fig.savefig(output / f"{sample}_FibK14_regional_statistics_3x5.{ext}", facecolor="white", bbox_inches="tight", **kwargs)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--mask-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=["extract", "plot", "all"], default="all")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "text.color": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.edgecolor": "black",
    })

    cache_path = args.output / "all_19dpb_spatial_features.tsv.gz"
    if args.mode in {"extract", "all"}:
        if args.h5ad is None or args.metadata is None or args.mask_dir is None:
            parser.error("--h5ad, --metadata and --mask-dir are required in extract/all mode")
        print("Loading registered 19-dpb data...", flush=True)
        data, audit = load_analysis(args.h5ad, args.metadata, args.mask_dir)
        fib = data.loc[data["sub_labels"] == "Fib_K14"].copy()
        blocks = make_spatial_blocks(fib)
        print(f"Fib_K14 bins: {len(fib)}; retained spatial blocks: {len(blocks)}", flush=True)
        summary, overall, pairwise = calculate_statistics(blocks)
        cache_columns = [
            "sample_key", "sample_batch_new", "sub_labels", "wound_x_um",
            "surface_axis", "depth_axis", "wound_territory", *FEATURES,
        ]
        data[cache_columns].to_csv(cache_path, sep="\t", index_label="obs_name", compression="gzip")
        blocks.to_csv(args.output / "Fib_K14_spatial_block_expression.tsv.gz", sep="\t", index=False, compression="gzip")
        summary.to_csv(args.output / "Fib_K14_regional_block_summary.csv", index=False)
        overall.to_csv(args.output / "Fib_K14_regional_Kruskal_Wallis_FDR.csv", index=False)
        pairwise.to_csv(args.output / "Fib_K14_regional_pairwise_MWU_FDR.csv", index=False)
        audit["statistics"] = {
            "unit": "Fib_K14 spatial-block mean",
            "block_definition": "300 um wound-x by 0.1 relative-depth; >=5 Fib_K14 bins",
            "overall": "Kruskal-Wallis; BH across 15 features within each section",
            "pairwise": "two-sided Mann-Whitney U; BH across 45 comparisons within each section",
            "interpretation": "within-section spatial heterogeneity; not independent-patient biological replication",
        }
        (args.output / "analysis_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        if args.mode == "extract":
            print(f"Extraction complete: {cache_path}", flush=True)
            return
    else:
        print("Loading extracted plotting cache...", flush=True)
        data = pd.read_csv(cache_path, sep="\t", compression="gzip", index_col="obs_name")
        blocks = pd.read_csv(args.output / "Fib_K14_spatial_block_expression.tsv.gz", sep="\t", compression="gzip")
        summary = pd.read_csv(args.output / "Fib_K14_regional_block_summary.csv")
        overall = pd.read_csv(args.output / "Fib_K14_regional_Kruskal_Wallis_FDR.csv")
        pairwise = pd.read_csv(args.output / "Fib_K14_regional_pairwise_MWU_FDR.csv")

    display_limits = {}
    for feature in FEATURES:
        vals = data[feature].to_numpy(float)
        if feature in GENES:
            positive = vals[vals > 0]
            lo = float(np.quantile(positive, 0.01)) if len(positive) else 0.0
            hi = float(np.quantile(positive, 0.99)) if len(positive) else 1.0
        else:
            lo, hi = np.quantile(vals[np.isfinite(vals)], [0.01, 0.99]).astype(float)
        if hi <= lo:
            hi = lo + 1.0
        display_limits[feature] = [lo, hi]

    for sample in SAMPLES:
        sample_df = data.loc[data["sample_key"] == sample].copy()
        plot_spatial(sample_df, sample, args.output, display_limits)
        plot_statistics(blocks, summary, overall, pairwise, sample, args.output)

    audit_path = args.output / "analysis_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["display_limits_1st_99th_percentile"] = display_limits
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
