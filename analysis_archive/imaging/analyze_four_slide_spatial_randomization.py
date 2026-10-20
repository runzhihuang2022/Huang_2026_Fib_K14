from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import chi2


BASE = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap"
)
OUTPUT = BASE / "08_four_slide_spatial_statistics"
OUTPUT.mkdir(parents=True, exist_ok=True)
PIXEL_SIZE_UM = 1.0113780498504639
RADII_UM = np.array([10.0, 25.0, 50.0, 100.0, 200.0])
N_PERMUTATIONS = 1999
RNG_SEED = 20260820

SLIDES = [
    ("BW3", BASE / "BW3_reprocess_20260820" / "03_imagej_double_positive_spatial", "pass"),
    ("BW7", BASE / "BW7_reprocess_20260820" / "03_imagej_double_positive_spatial", "pass"),
    ("BW9", BASE / "ImageJ_double_positive_spatial", "pass"),
    ("001A", BASE / "001A_reprocess_20260820" / "03_imagej_double_positive_spatial", "review_required"),
]


def read_cells(path: Path) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    coords = frame[["x_px", "y_px"]].to_numpy(dtype=np.float64)
    positive = frame["double_positive"].astype(str).str.lower().eq("true").to_numpy()
    return coords, coords[positive]


def symmetric_proximity(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, float]:
    if len(a) == 0 or len(b) == 0:
        return np.full(len(RADII_UM), np.nan), np.nan
    da = cKDTree(b).query(a, k=1)[0] * PIXEL_SIZE_UM
    db = cKDTree(a).query(b, k=1)[0] * PIXEL_SIZE_UM
    proximity = np.array([50.0 * ((da <= radius).mean() + (db <= radius).mean()) for radius in RADII_UM])
    median_distance = 0.5 * (np.median(da) + np.median(db))
    return proximity, float(median_distance)


def holm_adjust(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, index in enumerate(order):
        value = min(1.0, (m - rank) * p[index])
        running = max(running, value)
        adjusted[index] = running
    return adjusted.tolist()


def fisher_combine(pvalues: list[float]) -> float:
    values = np.clip(np.asarray(pvalues, dtype=float), 1e-300, 1.0)
    statistic = -2.0 * np.log(values).sum()
    return float(chi2.sf(statistic, 2 * len(values)))


def exact_sign_flip_p(effects: np.ndarray) -> float:
    observed = float(np.mean(effects))
    values = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(effects)):
        values.append(float(np.mean(effects * np.asarray(signs))))
    return float(np.mean(np.asarray(values) >= observed - 1e-12))


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict] = []
    slide_objects = []
    for sample_id, root, registration_status in SLIDES:
        kc_all, kc_dp = read_cells(root / "03_double_positive_cells" / "KC_KRT6B_SOX9_per_cell_measurements.tsv")
        fib_all, fib_dp = read_cells(root / "03_double_positive_cells" / "Fib_PDGFRa_KRT14_per_cell_measurements.tsv")
        observed, observed_median = symmetric_proximity(kc_dp, fib_dp)
        null = np.empty((N_PERMUTATIONS, len(RADII_UM)), dtype=np.float64)
        null_median = np.empty(N_PERMUTATIONS, dtype=np.float64)
        for iteration in range(N_PERMUTATIONS):
            kc_random = kc_all[rng.choice(len(kc_all), size=len(kc_dp), replace=False)]
            fib_random = fib_all[rng.choice(len(fib_all), size=len(fib_dp), replace=False)]
            null[iteration], null_median[iteration] = symmetric_proximity(kc_random, fib_random)
        for index, radius in enumerate(RADII_UM):
            null_mean = float(np.mean(null[:, index]))
            p_value = float((1 + np.sum(null[:, index] >= observed[index] - 1e-12)) / (N_PERMUTATIONS + 1))
            rows.append(
                {
                    "sample_id": sample_id,
                    "registration_status": registration_status,
                    "radius_um": radius,
                    "KC_double_positive_n": len(kc_dp),
                    "Fib_double_positive_n": len(fib_dp),
                    "observed_symmetric_proximity_pct": observed[index],
                    "null_mean_pct": null_mean,
                    "null_95ci_low_pct": float(np.quantile(null[:, index], 0.025)),
                    "null_95ci_high_pct": float(np.quantile(null[:, index], 0.975)),
                    "effect_observed_minus_null_pct_points": observed[index] - null_mean,
                    "fold_enrichment": observed[index] / null_mean if null_mean > 0 else np.nan,
                    "permutation_p_one_sided": p_value,
                    "permutations": N_PERMUTATIONS,
                }
            )
        distance_p = float((1 + np.sum(null_median <= observed_median + 1e-12)) / (N_PERMUTATIONS + 1))
        slide_objects.append(
            {
                "sample_id": sample_id,
                "registration_status": registration_status,
                "observed": observed,
                "null_mean": null.mean(axis=0),
                "observed_median_distance_um": observed_median,
                "null_median_distance_mean_um": float(null_median.mean()),
                "median_distance_permutation_p_one_sided": distance_p,
            }
        )
        print(f"DONE {sample_id}: KC={len(kc_dp)} Fib={len(fib_dp)}", flush=True)

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT / "per_slide_spatial_label_randomization.csv", index=False, encoding="utf-8-sig")
    primary = result[result["radius_um"] == 10.0].copy()
    effects_all = primary["effect_observed_minus_null_pct_points"].to_numpy()
    effects_pass = primary.loc[primary["registration_status"] == "pass", "effect_observed_minus_null_pct_points"].to_numpy()
    rng_boot = np.random.default_rng(RNG_SEED + 1)
    boot_all = np.array([rng_boot.choice(effects_all, size=len(effects_all), replace=True).mean() for _ in range(10000)])
    boot_pass = np.array([rng_boot.choice(effects_pass, size=len(effects_pass), replace=True).mean() for _ in range(10000)])

    combined_rows = []
    for radius in RADII_UM:
        subset = result[result["radius_um"] == radius]
        all_p = fisher_combine(subset["permutation_p_one_sided"].tolist())
        pass_p = fisher_combine(subset.loc[subset["registration_status"] == "pass", "permutation_p_one_sided"].tolist())
        combined_rows.append({"radius_um": radius, "fisher_p_all4": all_p, "fisher_p_pass3": pass_p})
    combined = pd.DataFrame(combined_rows)
    combined["holm_adjusted_fisher_p_all4"] = holm_adjust(combined["fisher_p_all4"].tolist())
    combined["holm_adjusted_fisher_p_pass3"] = holm_adjust(combined["fisher_p_pass3"].tolist())
    combined.to_csv(OUTPUT / "combined_spatial_randomization_pvalues.csv", index=False, encoding="utf-8-sig")

    primary_summary = {
        "primary_radius_um": 10.0,
        "independent_unit": "registered serial-section pair (slide/sample)",
        "independent_n_all": 4,
        "independent_n_registration_pass": 3,
        "mean_effect_all4_percentage_points": float(effects_all.mean()),
        "bootstrap_95ci_all4": [float(np.quantile(boot_all, 0.025)), float(np.quantile(boot_all, 0.975))],
        "exact_sign_flip_p_one_sided_all4": exact_sign_flip_p(effects_all),
        "mean_effect_pass3_percentage_points": float(effects_pass.mean()),
        "bootstrap_95ci_pass3": [float(np.quantile(boot_pass, 0.025)), float(np.quantile(boot_pass, 0.975))],
        "exact_sign_flip_p_one_sided_pass3": exact_sign_flip_p(effects_pass),
        "epi_front_analysis": "not_tested_missing_user_defined_ROI",
        "registration_sensitivity": "001A excluded from pass-only sensitivity because registration_status=review_required",
        "interpretation_rule": "p<0.05 required; no result is upgraded to proof; serial sections do not establish same-cell colocalization",
    }
    (OUTPUT / "primary_slide_level_inference.json").write_text(
        json.dumps(primary_summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    x = np.arange(len(SLIDES))
    width = 0.36
    axes[0].bar(x - width / 2, primary["observed_symmetric_proximity_pct"], width, label="Observed", color="#6A3D9A")
    axes[0].bar(x + width / 2, primary["null_mean_pct"], width, label="Random-label null", color="#BDBDBD")
    axes[0].set_xticks(x, primary["sample_id"])
    axes[0].set_ylabel("Symmetric proximity within 10 µm (%)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title("Whole-tissue primary endpoint", fontsize=10, fontweight="bold")
    for i, row in primary.reset_index(drop=True).iterrows():
        axes[0].text(i, max(row["observed_symmetric_proximity_pct"], row["null_mean_pct"]) + 0.08, f"p={row['permutation_p_one_sided']:.3g}", ha="center", fontsize=7)

    colors = {"BW3": "#1B9E77", "BW7": "#D95F02", "BW9": "#7570B3", "001A": "#666666"}
    for sample_id in primary["sample_id"]:
        subset = result[result["sample_id"] == sample_id]
        axes[1].plot(subset["radius_um"], subset["effect_observed_minus_null_pct_points"], marker="o", label=sample_id, color=colors[sample_id])
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].set_xscale("log")
    axes[1].set_xticks(RADII_UM, [str(int(v)) for v in RADII_UM])
    axes[1].set_xlabel("Radius (µm)")
    axes[1].set_ylabel("Observed − null (percentage points)")
    axes[1].set_title("Spatial enrichment sensitivity", fontsize=10, fontweight="bold")
    axes[1].legend(frameon=False, fontsize=7, ncol=2)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=8)
    fig.savefig(OUTPUT / "whole_tissue_spatial_randomization.png", dpi=300)
    fig.savefig(OUTPUT / "whole_tissue_spatial_randomization.pdf")
    plt.close(fig)

    with (OUTPUT / "epi_front_roi_manifest_template.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "roi_file_or_coordinates", "coordinate_space", "reviewer", "status", "note"])
        writer.writeheader()
        for sample_id, _, registration_status in SLIDES:
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "roi_file_or_coordinates": "",
                    "coordinate_space": "registered_analysis_image",
                    "reviewer": "",
                    "status": "AUTHOR_INPUT_NEEDED",
                    "note": f"User-defined Epi-Front ROI required; registration={registration_status}",
                }
            )

    report = f"""# Four-slide spatial statistics audit

## Scope

- Samples: BW3, BW7, BW9, 001A.
- Independent unit: one registered serial-section pair per sample (`n=4`; registration-pass sensitivity `n=3`).
- Primary endpoint: symmetric percentage of the two double-positive cell sets within 10 µm.
- Null model: independently permute double-positive labels among all DAPI-segmented cells within each registered section, preserving cell positions, cell counts, and tissue-density structure; {N_PERMUTATIONS} permutations per slide.
- Secondary radii: 25, 50, 100 and 200 µm; Fisher combined p values are Holm-adjusted across radii.

## Primary slide-level inference

- All four slides: mean observed-minus-null effect = {primary_summary['mean_effect_all4_percentage_points']:.4f} percentage points; bootstrap 95% CI {primary_summary['bootstrap_95ci_all4'][0]:.4f} to {primary_summary['bootstrap_95ci_all4'][1]:.4f}; exact one-sided sign-flip p = {primary_summary['exact_sign_flip_p_one_sided_all4']:.4g}.
- Registration-pass slides only: mean effect = {primary_summary['mean_effect_pass3_percentage_points']:.4f} percentage points; bootstrap 95% CI {primary_summary['bootstrap_95ci_pass3'][0]:.4f} to {primary_summary['bootstrap_95ci_pass3'][1]:.4f}; exact one-sided sign-flip p = {primary_summary['exact_sign_flip_p_one_sided_pass3']:.4g}.

## Claim status

- Whole-tissue significant colocalization: determined from the tests above; do not use cell count as independent n.
- Epi-Front enrichment: **not tested**, because no user-reviewed Epi-Front ROI or boundary was provided.
- Same-cell molecular colocalization: not established because the markers are on adjacent serial sections.
- 001A: registration review required (Dice 0.720865), therefore separated as a sensitivity condition.

## Required next input

Complete `epi_front_roi_manifest_template.csv` using ImageJ ROI files or registered-image coordinates for all four samples. The Epi-Front analysis must use those locked regions without post-hoc movement based on marker signal.
"""
    (OUTPUT / "STATISTICAL_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(primary_summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
