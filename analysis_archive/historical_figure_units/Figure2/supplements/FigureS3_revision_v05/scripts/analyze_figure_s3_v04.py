from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
QC = ROOT / "QC"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

STEREO_ROOT = Path(
    os.environ.get(
        "FIGS3_STEREO_ROOT",
        r"data/external/local_91ac4b6618ea"
        r"\时空组学\stereoseq\pdf_output\Wound_Healing_Annotation_Output_White",
    )
)

SECTIONS = [
    ("Normal", "NS_C02847B1", "Normal", 0),
    ("5dpb", "BW32_A01597A3_SDSDB_5dpb", "5 dpb", 1),
    ("12dpb_DPTDI1", "BW13_1_B1_DSDB_12dpb", "12 dpb", 2),
    ("12dpb_SPTDI1", "BW14_1_C1_SSDB_12dpb", "12 dpb", 2),
    ("12dpb_DPTDI2", "BW13_A3_DSDB_12dpb", "12 dpb", 2),
    ("12dpb_SPTDI2", "BW14_B3_SSDB_12dpb", "12 dpb", 2),
    ("19dpb", "BW15D_C6_SDSDB_19dpb", "19 dpb", 3),
    ("19dpb_p1", "BW15D_1_D1_SDSDB_19dpb", "19 dpb", 3),
    ("26dpb_p1", "BW81_C02846B6_SDSDB_26dpb_part1", "26 dpb", 4),
    ("26dpb_p2", "BW81_C02846B6_SDSDB_26dpb_part2", "26 dpb", 4),
    ("2mph", "2mph_A03699G6.SCT", "2 mph", 5),
]

TARGET = "Fib_K14"
DISPLAY_STATES = [TARGET, "KC_Spinous_Mig", "SAC_SG_Progenitor"]
ALL_TESTED_STATES = [
    TARGET,
    "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous",
    "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
    "SAC_SG_Progenitor", "SAC_SG_Clear", "SAC_SG_Dark", "SAC_SG_Ductal",
    "SAC_HF_IRS", "SAC_HF_ORS", "SAC_HF_HFSC", "SAC_HF_DP_DS",
    "SAC_HF_Matrix",
]
N_PERMUTATIONS = 999
SEED = 20260724
PIXEL_TO_UM = 0.5
FIXED_RADIUS_UM = 50.0


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def locate_h5ad(sample_id: str) -> Path:
    exact = STEREO_ROOT / f"{sample_id}_Annotated.h5ad"
    if exact.exists():
        return exact
    matches = list(STEREO_ROOT.glob(f"{sample_id}*Annotated.h5ad"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one annotated H5AD for {sample_id}; found {matches}")
    return matches[0]


def read_section(sample_id: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = locate_h5ad(sample_id)
    with h5py.File(path, "r") as handle:
        group = handle["obs"]["sub_labels"]
        categories = np.asarray(group["categories"]).astype("U")
        codes = np.asarray(group["codes"], dtype=int)
        labels = categories[codes]
        coords = np.asarray(handle["obsm"]["spatial"], dtype=float)
    return coords, labels, categories


def score_from_labels(
    graph: csr_matrix,
    degrees: np.ndarray,
    codes: np.ndarray,
    target_code: int,
    candidate_codes: np.ndarray,
) -> np.ndarray:
    target_indicator = (codes == target_code).astype(np.float64)
    neighbors_of_target = graph @ target_indicator
    target_degree = float(degrees[codes == target_code].sum())
    total_degree = float(degrees.sum())
    if target_degree == 0 or total_degree == 0:
        return np.full(len(candidate_codes), np.nan)
    conditional = (
        np.bincount(codes, weights=neighbors_of_target, minlength=int(codes.max()) + 1)
        / target_degree
    )
    marginal = (
        np.bincount(codes, weights=degrees, minlength=int(codes.max()) + 1)
        / total_degree
    )
    return np.divide(
        conditional[candidate_codes],
        marginal[candidate_codes],
        out=np.zeros(len(candidate_codes), dtype=float),
        where=marginal[candidate_codes] > 0,
    )


def fixed_radius_statistics() -> tuple[pd.DataFrame, pd.DataFrame]:
    section_rows: list[dict] = []
    null_by_time: dict[tuple[str, str], list[np.ndarray]] = {}

    for section_index, (section, sample_id, timepoint, time_order) in enumerate(SECTIONS):
        coords, labels, _ = read_section(sample_id)
        threshold = FIXED_RADIUS_UM / PIXEL_TO_UM
        edges = cKDTree(coords).query_pairs(threshold, output_type="ndarray")
        if not len(edges):
            raise RuntimeError(f"No spatial edges found for {section}")
        row = np.r_[edges[:, 0], edges[:, 1]]
        col = np.r_[edges[:, 1], edges[:, 0]]
        graph = csr_matrix(
            (np.ones(len(row), dtype=np.float32), (row, col)),
            shape=(len(labels), len(labels)),
        )
        degrees = np.asarray(graph.sum(axis=1)).ravel()

        state_to_code = {
            "Other": 0,
            **{state: index + 1 for index, state in enumerate(ALL_TESTED_STATES)},
        }
        codes = np.zeros(len(labels), dtype=np.int16)
        for state in ALL_TESTED_STATES:
            codes[labels == state] = state_to_code[state]
        candidate_codes = np.array(
            [state_to_code[state] for state in ALL_TESTED_STATES]
        )
        observed = score_from_labels(graph, degrees, codes, state_to_code[TARGET], candidate_codes)

        rng = np.random.default_rng(SEED + section_index)
        null = np.empty((N_PERMUTATIONS, len(ALL_TESTED_STATES)), dtype=float)
        for permutation in range(N_PERMUTATIONS):
            shuffled = rng.permutation(codes)
            null[permutation] = score_from_labels(
                graph, degrees, shuffled, state_to_code[TARGET], candidate_codes
            )
        p_values = (1 + np.sum(null >= observed[None, :], axis=0)) / (N_PERMUTATIONS + 1)
        q_values = bh_adjust(p_values)

        for state_index, state in enumerate(ALL_TESTED_STATES):
            section_rows.append(
                {
                    "section": section,
                    "sample_id": sample_id,
                    "timepoint": timepoint,
                    "time_order": time_order,
                    "conditioned_on": TARGET,
                    "cell_state": state,
                    "comparison_type": "self" if state == TARGET else "cross-state",
                    "n_bins": len(labels),
                    "n_Fib_K14": int(np.sum(labels == TARGET)),
                    "n_cell_state": int(np.sum(labels == state)),
                    "fixed_radius_um": FIXED_RADIUS_UM,
                    "cooccurrence_score": observed[state_index],
                    "p_enrichment": p_values[state_index],
                    "fdr_within_17_FibK14_KC_SAC_states": q_values[state_index],
                    "n_permutations": N_PERMUTATIONS,
                }
            )
            null_by_time.setdefault((timepoint, state), []).append(null[:, state_index])
        print(f"{section}: n={len(labels)}, edges={len(edges)}", flush=True)

    section_stats = pd.DataFrame(section_rows)
    time_rows: list[dict] = []
    for timepoint, time_order in [
        ("Normal", 0), ("5 dpb", 1), ("12 dpb", 2),
        ("19 dpb", 3), ("26 dpb", 4), ("2 mph", 5),
    ]:
        for state in ALL_TESTED_STATES:
            subset = section_stats[
                (section_stats["timepoint"] == timepoint)
                & (section_stats["cell_state"] == state)
            ]
            observed_mean = float(subset["cooccurrence_score"].mean())
            time_null = np.mean(np.stack(null_by_time[(timepoint, state)]), axis=0)
            p_value = float(
                (1 + np.sum(time_null >= observed_mean)) / (N_PERMUTATIONS + 1)
            )
            time_rows.append(
                {
                    "timepoint": timepoint,
                    "time_order": time_order,
                    "cell_state": state,
                    "comparison_type": "self" if state == TARGET else "cross-state",
                    "n_sections": len(subset),
                    "mean_score": observed_mean,
                    "sd_score": (
                        float(subset["cooccurrence_score"].std(ddof=1))
                        if len(subset) > 1 else np.nan
                    ),
                    "sem_score": (
                        float(subset["cooccurrence_score"].sem(ddof=1))
                        if len(subset) > 1 else np.nan
                    ),
                    "p_enrichment_timepoint_permutation": p_value,
                    "n_permutations": N_PERMUTATIONS,
                }
            )
    time_stats = pd.DataFrame(time_rows)
    time_stats["fdr_across_102_timepoint_state_tests"] = bh_adjust(
        time_stats["p_enrichment_timepoint_permutation"].to_numpy()
    )
    section_stats.to_csv(
        OUT / "FigureS3_v04_prespecified_section_statistics.csv", index=False
    )
    time_stats.to_csv(
        OUT / "FigureS3_v04_prespecified_timepoint_statistics.csv", index=False
    )
    return section_stats, time_stats


def distance_cooccurrence_26dpb_p2() -> tuple[pd.DataFrame, pd.DataFrame]:
    coords, labels, categories = read_section("BW81_C02846B6_SDSDB_26dpb_part2")
    coordinate_sum = coords.sum(axis=1)
    two_min = np.argpartition(coordinate_sum, 2)[:2]
    min_index = two_min[np.argmin(coordinate_sum[two_min])]
    second_min = two_min[1] if two_min[0] == min_index else two_min[0]
    max_index = int(np.argmax(coordinate_sum))
    threshold_min = float(np.linalg.norm(coords[min_index] - coords[second_min]))
    threshold_max = float(np.linalg.norm(coords[min_index] - coords[max_index]) / 2)
    interval = np.linspace(threshold_min, threshold_max, 50)
    radii = interval[1:]

    all_tree = cKDTree(coords)
    target_mask = labels == TARGET
    target_tree = cKDTree(coords[target_mask])
    target_count = int(target_mask.sum())
    all_counts = []
    target_counts = []
    for category in categories:
        category_mask = labels == category
        category_tree = cKDTree(coords[category_mask])
        marginal_counts = (
            all_tree.count_neighbors(category_tree, radii, cumulative=True).astype(float)
            - int(category_mask.sum())
        )
        conditional_counts = target_tree.count_neighbors(
            category_tree, radii, cumulative=True
        ).astype(float)
        if category == TARGET:
            conditional_counts -= target_count
        all_counts.append(marginal_counts)
        target_counts.append(conditional_counts)
    all_counts = np.asarray(all_counts)
    target_counts = np.asarray(target_counts)
    conditional = target_counts / target_counts.sum(axis=0)
    marginal = all_counts / all_counts.sum(axis=0)
    ratio = np.divide(
        conditional,
        marginal,
        out=np.zeros_like(conditional),
        where=marginal > 0,
    )

    curve_rows = []
    for category_index, category in enumerate(categories):
        for distance_index, distance_coordinate_units in enumerate(radii):
            curve_rows.append(
                {
                    "section": "26dpb_p2",
                    "conditioned_on": TARGET,
                    "cell_state": category,
                    "distance_coordinate_units": distance_coordinate_units,
                    "distance_um": distance_coordinate_units * PIXEL_TO_UM,
                    "cooccurrence_probability_ratio": ratio[
                        category_index, distance_index
                    ],
                }
            )
    curves = pd.DataFrame(curve_rows)
    closest = pd.DataFrame(
        {
            "cell_state": categories,
            "closest_distance_um": radii[0] * PIXEL_TO_UM,
            "closest_distance_score": ratio[:, 0],
        }
    ).sort_values("closest_distance_score", ascending=False)
    closest["rank_including_Fib_K14_self"] = np.arange(1, len(closest) + 1)
    closest["is_highlighted"] = closest["cell_state"].isin(
        [TARGET, "KC_Spinous_Mig"]
    )
    curves.to_csv(OUT / "FigureS3_v04_26dpb_p2_distance_curves.csv", index=False)
    closest.to_csv(OUT / "FigureS3_v04_26dpb_p2_closest_scores.csv", index=False)
    return curves, closest


def main() -> None:
    section_stats, time_stats = fixed_radius_statistics()
    curves, closest = distance_cooccurrence_26dpb_p2()
    summary_lines = [
        "Figure S3 v04 statistical audit",
        "================================",
        "",
        "B: 999 spatial-label permutations; 50-um fixed radius.",
        "Time-point P values compare observed mean scores with mean permuted scores.",
        "BH correction was applied across Fib_K14 plus all 16 KC/SAC states x 6 time points.",
        "",
        time_stats.to_string(index=False),
        "",
        "C: distance curves exactly reproduce the Squidpy cumulative-radius ratio",
        "using cKDTree cumulative pair counts and physical conversion 0.5 um/unit.",
        "",
        closest.head(10).to_string(index=False),
    ]
    (QC / "FigureS3_v04_statistical_audit.txt").write_text(
        "\n".join(summary_lines), encoding="utf-8"
    )
    print(section_stats.shape, time_stats.shape, curves.shape, closest.shape)


if __name__ == "__main__":
    main()
