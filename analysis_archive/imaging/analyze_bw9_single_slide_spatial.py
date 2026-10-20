from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


BASE = Path(
    r"data/external/local_3413b1de7846"
    r"\20260813找曾老师指导\免疫荧光结果统计\Figure 3G overlap"
)
INPUT = BASE / "ImageJ_double_positive_spatial" / "03_double_positive_cells"
OUTPUT = BASE / "09_BW9_single_slide_spatial_statistics"
PIXEL_SIZE_UM = 1.0113780498504639
RADII_UM = np.array([10.0, 25.0, 50.0, 100.0, 200.0])
N_PERMUTATIONS = 9999
RNG_SEED = 20260821


def read_cells(filename: str) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    frame = pd.read_csv(INPUT / filename, sep="\t", encoding="utf-8-sig")
    coords = frame[["x_px", "y_px"]].to_numpy(dtype=np.float64)
    positive = frame["double_positive"].astype(str).str.lower().eq("true").to_numpy()
    return coords, coords[positive], frame


def proximity(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    da = cKDTree(b).query(a, k=1)[0] * PIXEL_SIZE_UM
    db = cKDTree(a).query(b, k=1)[0] * PIXEL_SIZE_UM
    a_pct = np.array([(da <= radius).mean() * 100.0 for radius in RADII_UM])
    b_pct = np.array([(db <= radius).mean() * 100.0 for radius in RADII_UM])
    symmetric_pct = 0.5 * (a_pct + b_pct)
    return a_pct, b_pct, symmetric_pct


def holm_adjust(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues)
    adjusted = np.empty_like(pvalues, dtype=float)
    running = 0.0
    m = len(pvalues)
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvalues[index]))
        adjusted[index] = running
    return adjusted


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    kc_all, kc_dp, kc_frame = read_cells("KC_KRT6B_SOX9_per_cell_measurements.tsv")
    fib_all, fib_dp, fib_frame = read_cells("Fib_PDGFRa_KRT14_per_cell_measurements.tsv")

    obs_kc, obs_fib, obs_sym = proximity(kc_dp, fib_dp)
    observed_kc_dist = cKDTree(fib_dp).query(kc_dp, k=1)[0] * PIXEL_SIZE_UM
    observed_fib_dist = cKDTree(kc_dp).query(fib_dp, k=1)[0] * PIXEL_SIZE_UM
    observed_median = 0.5 * (np.median(observed_kc_dist) + np.median(observed_fib_dist))

    rng = np.random.default_rng(RNG_SEED)
    null_kc = np.empty((N_PERMUTATIONS, len(RADII_UM)))
    null_fib = np.empty_like(null_kc)
    null_sym = np.empty_like(null_kc)
    null_median = np.empty(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        kc_random = kc_all[rng.choice(len(kc_all), size=len(kc_dp), replace=False)]
        fib_random = fib_all[rng.choice(len(fib_all), size=len(fib_dp), replace=False)]
        null_kc[i], null_fib[i], null_sym[i] = proximity(kc_random, fib_random)
        kd = cKDTree(fib_random).query(kc_random, k=1)[0] * PIXEL_SIZE_UM
        fd = cKDTree(kc_random).query(fib_random, k=1)[0] * PIXEL_SIZE_UM
        null_median[i] = 0.5 * (np.median(kd) + np.median(fd))
        if (i + 1) % 1000 == 0:
            print(f"Permutation {i + 1}/{N_PERMUTATIONS}", flush=True)

    p_sym = (1 + np.sum(null_sym >= obs_sym[None, :] - 1e-12, axis=0)) / (N_PERMUTATIONS + 1)
    p_kc = (1 + np.sum(null_kc >= obs_kc[None, :] - 1e-12, axis=0)) / (N_PERMUTATIONS + 1)
    p_fib = (1 + np.sum(null_fib >= obs_fib[None, :] - 1e-12, axis=0)) / (N_PERMUTATIONS + 1)
    secondary_holm = np.full(len(RADII_UM), np.nan)
    secondary_holm[1:] = holm_adjust(p_sym[1:])

    rows = []
    for j, radius in enumerate(RADII_UM):
        null_mean = null_sym[:, j].mean()
        rows.append(
            {
                "radius_um": radius,
                "KC_direction_observed_pct": obs_kc[j],
                "KC_direction_null_mean_pct": null_kc[:, j].mean(),
                "KC_direction_p_one_sided": p_kc[j],
                "Fib_direction_observed_pct": obs_fib[j],
                "Fib_direction_null_mean_pct": null_fib[:, j].mean(),
                "Fib_direction_p_one_sided": p_fib[j],
                "symmetric_observed_pct": obs_sym[j],
                "symmetric_null_mean_pct": null_mean,
                "symmetric_null_95ci_low_pct": np.quantile(null_sym[:, j], 0.025),
                "symmetric_null_95ci_high_pct": np.quantile(null_sym[:, j], 0.975),
                "effect_observed_minus_null_percentage_points": obs_sym[j] - null_mean,
                "fold_enrichment": obs_sym[j] / null_mean if null_mean > 0 else np.nan,
                "permutation_p_one_sided": p_sym[j],
                "secondary_holm_adjusted_p": secondary_holm[j],
                "is_prespecified_primary_radius": radius == 10.0,
                "permutations": N_PERMUTATIONS,
            }
        )
    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT / "BW9_spatial_randomization_by_radius.csv", index=False, encoding="utf-8-sig")

    null_source = pd.DataFrame(
        {
            "permutation_id": np.arange(1, N_PERMUTATIONS + 1),
            "symmetric_proximity_10um_pct": null_sym[:, 0],
            "symmetric_median_nearest_distance_um": null_median,
        }
    )
    null_source.to_csv(OUTPUT / "BW9_null_distribution_source_data.csv", index=False, encoding="utf-8-sig")

    distance_p = (1 + np.sum(null_median <= observed_median + 1e-12)) / (N_PERMUTATIONS + 1)
    summary = {
        "sample_id": "BW9",
        "independent_experimental_unit": "one registered serial-section pair",
        "independent_n": 1,
        "cell_counts_are_not_independent_replicates": True,
        "KC_DAPI_nuclei": int(len(kc_all)),
        "KC_KRT6B_SOX9_double_positive_n": int(len(kc_dp)),
        "KC_double_positive_percent": float(100 * len(kc_dp) / len(kc_all)),
        "Fib_DAPI_nuclei": int(len(fib_all)),
        "Fib_PDGFRa_KRT14_double_positive_n": int(len(fib_dp)),
        "Fib_double_positive_percent": float(100 * len(fib_dp) / len(fib_all)),
        "primary_radius_um": 10.0,
        "primary_observed_symmetric_proximity_pct": float(obs_sym[0]),
        "primary_null_mean_pct": float(null_sym[:, 0].mean()),
        "primary_effect_percentage_points": float(obs_sym[0] - null_sym[:, 0].mean()),
        "primary_fold_enrichment": float(obs_sym[0] / null_sym[:, 0].mean()),
        "primary_permutation_p_one_sided": float(p_sym[0]),
        "observed_symmetric_median_nearest_distance_um": float(observed_median),
        "null_mean_symmetric_median_nearest_distance_um": float(null_median.mean()),
        "median_distance_permutation_p_one_sided": float(distance_p),
        "randomization_model": "independent double-positive label permutation among all DAPI-segmented cells within each registered section",
        "permutations": N_PERMUTATIONS,
        "random_seed": RNG_SEED,
        "pixel_size_um": PIXEL_SIZE_UM,
        "epi_front_status": "not_tested_missing_user_defined_ROI",
        "same_cell_colocalization_status": "not_applicable_adjacent_serial_sections",
    }
    (OUTPUT / "BW9_primary_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    plt.rcParams.update({"font.family": "Arial", "font.size": 8})
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)
    axes[0].hist(null_sym[:, 0], bins=40, color="#BDBDBD", edgecolor="white", linewidth=0.3)
    axes[0].axvline(obs_sym[0], color="#6A3D9A", linewidth=1.5, label=f"Observed: {obs_sym[0]:.2f}%")
    axes[0].axvline(null_sym[:, 0].mean(), color="black", linestyle="--", linewidth=1, label=f"Null mean: {null_sym[:, 0].mean():.2f}%")
    axes[0].set_xlabel("Symmetric proximity within 10 µm (%)")
    axes[0].set_ylabel("Permutations")
    axes[0].set_title("Primary randomization test", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=7)

    low = np.quantile(null_sym, 0.025, axis=0)
    high = np.quantile(null_sym, 0.975, axis=0)
    mean = null_sym.mean(axis=0)
    axes[1].fill_between(RADII_UM, low, high, color="#BDBDBD", alpha=0.45, label="Null 95% interval")
    axes[1].plot(RADII_UM, mean, "--", color="black", marker="o", markersize=3, label="Null mean")
    axes[1].plot(RADII_UM, obs_sym, color="#6A3D9A", marker="o", markersize=3, label="Observed")
    axes[1].set_xscale("log")
    axes[1].set_xticks(RADII_UM, [str(int(x)) for x in RADII_UM])
    axes[1].set_xlabel("Radius (µm)")
    axes[1].set_ylabel("Symmetric proximity (%)")
    axes[1].set_title("Radius sensitivity", fontweight="bold")
    axes[1].legend(frameon=False, fontsize=7)

    effects = obs_sym - mean
    axes[2].axhline(0, color="black", linewidth=0.7)
    axes[2].plot(RADII_UM, effects, color="#1B9E77", marker="o", markersize=3)
    axes[2].set_xscale("log")
    axes[2].set_xticks(RADII_UM, [str(int(x)) for x in RADII_UM])
    axes[2].set_xlabel("Radius (µm)")
    axes[2].set_ylabel("Observed − null (percentage points)")
    axes[2].set_title("Spatial enrichment effect", fontweight="bold")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUTPUT / "BW9_single_slide_spatial_statistics.png", dpi=300)
    fig.savefig(OUTPUT / "BW9_single_slide_spatial_statistics.pdf")
    plt.close(fig)

    primary = results.iloc[0]
    conclusion = (
        "支持 BW9 全切片范围内的图像内空间富集"
        if primary["permutation_p_one_sided"] < 0.05 and primary["effect_observed_minus_null_percentage_points"] > 0
        else "不支持 BW9 全切片范围内的图像内空间富集"
    )
    report = f"""# BW9 单张切片空间统计报告

## 分析范围与单位

- 样本：BW9，一对已经配准的相邻连续切片；独立实验单位 `n=1`。
- 细胞数不是生物学重复，不能用 740 个双阳性细胞替代样本量。
- KC 切片：KRT6B⁺SOX9⁺ {len(kc_dp)}/{len(kc_all)}（{100*len(kc_dp)/len(kc_all):.4f}%）。
- Fib 切片：PDGFRα⁺KRT14⁺ {len(fib_dp)}/{len(fib_all)}（{100*len(fib_dp)/len(fib_all):.4f}%）。

## 方法

- 主要终点：两类双阳性细胞在 10 µm 内的对称邻近比例，即两个方向邻近百分比的平均值。
- 零模型：分别在两张切片的全部 DAPI 分割细胞中随机置换双阳性标签，保持细胞位置、组织密度和双阳性细胞总数不变。
- 随机置换 {N_PERMUTATIONS} 次，单侧检验空间富集；25–200 µm 为敏感性分析，并对四个次要半径使用 Holm 校正。

## 结果

- 10 µm 观察值 {obs_sym[0]:.4f}%，随机零分布均值 {null_sym[:,0].mean():.4f}%（95% 区间 {np.quantile(null_sym[:,0],0.025):.4f}%–{np.quantile(null_sym[:,0],0.975):.4f}%）。
- 绝对富集效应 {obs_sym[0]-null_sym[:,0].mean():.4f} 个百分点，倍数富集 {obs_sym[0]/null_sym[:,0].mean():.3f} 倍，单侧置换 p={p_sym[0]:.4g}。
- 两方向观察值：KRT6B⁺SOX9⁺ 至 PDGFRα⁺KRT14⁺ 为 {obs_kc[0]:.4f}%；反方向为 {obs_fib[0]:.4f}%。
- 对称最近邻中位距离为 {observed_median:.4f} µm，随机零模型均值为 {null_median.mean():.4f} µm，距离缩短的单侧置换 p={distance_p:.4g}。

## 结论边界

- 结论：**{conclusion}**。
- 该结果只支持 BW9 这一张连续切片配对中的空间邻近富集，不构成跨样本总体推断，也不证明同一细胞内分子共定位或因果关系。
- Epi-Front：未检验；缺少用户锁定并审核的 Epi-Front ROI，未根据 marker 信号事后选择区域。
"""
    (OUTPUT / "BW9_STATISTICAL_REPORT_CN.md").write_text(report, encoding="utf-8")

    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
