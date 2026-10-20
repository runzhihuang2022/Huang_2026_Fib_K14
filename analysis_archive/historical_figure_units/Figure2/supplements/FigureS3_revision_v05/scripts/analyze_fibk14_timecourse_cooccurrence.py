from __future__ import annotations

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
H5AD = Path(
    r"data/external/local_91ac4b6618ea"
    r"\时空组学\stereoseq\combined_16samples_Harmony_Clustered_Fixed_Annotated_White.h5ad"
)

SAMPLES = {
    "Normal": ("NS_C02847B1", "Normal", 0),
    "5dpb": ("BW32_A01597A3_SDSDB_5dpb", "5 dpb", 1),
    "12dpb_DPTDI1": ("BW13_1_B1_DSDB_12dpb", "12 dpb", 2),
    "12dpb_SPTDI1": ("BW14_1_C1_SSDB_12dpb", "12 dpb", 2),
    "12dpb_DPTDI2": ("BW13_A3_DSDB_12dpb", "12 dpb", 2),
    "12dpb_SPTDI2": ("BW14_B3_SSDB_12dpb", "12 dpb", 2),
    "19dpb": ("BW15D_C6_SDSDB_19dpb", "19 dpb", 3),
    "19dpb_p1": ("BW15D_1_D1_SDSDB_19dpb", "19 dpb", 3),
    "26dpb_p1": ("BW81_C02846B6_SDSDB_26dpb_part1", "26 dpb", 4),
    "26dpb_p2": ("BW81_C02846B6_SDSDB_26dpb_part2", "26 dpb", 4),
    "2mph": ("2mph_A03699G6.SCT", "2 mph", 5),
}
TARGET = "Fib_K14"
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


def closest_threshold(coords: np.ndarray) -> tuple[float, float, float]:
    """Reproduce Squidpy's default first co-occurrence distance threshold."""
    sums = coords[:, 0] + coords[:, 1]
    two_min = np.argpartition(sums, 2)[:2]
    min_idx = two_min[np.argmin(sums[two_min])]
    min_idx2 = two_min[1] if two_min[0] == min_idx else two_min[0]
    max_idx = int(np.argmax(sums))
    threshold_min = float(np.linalg.norm(coords[min_idx] - coords[min_idx2]))
    threshold_max = float(np.linalg.norm(coords[min_idx] - coords[max_idx]) / 2.0)
    interval = np.linspace(threshold_min, threshold_max, 50)
    return float(interval[1]), threshold_min, threshold_max


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
    conditional = np.bincount(codes, weights=neighbors_of_target, minlength=int(codes.max()) + 1) / target_degree
    marginal = np.bincount(codes, weights=degrees, minlength=int(codes.max()) + 1) / total_degree
    return np.divide(
        conditional[candidate_codes],
        marginal[candidate_codes],
        out=np.zeros(len(candidate_codes), dtype=float),
        where=marginal[candidate_codes] > 0,
    )


def analyze_section(
    section: str,
    sample_id: str,
    timepoint: str,
    time_order: int,
    coords: np.ndarray,
    labels: np.ndarray,
    candidates: list[str],
    section_index: int,
) -> list[dict]:
    _, threshold_min, threshold_max = closest_threshold(coords)
    threshold = FIXED_RADIUS_UM / PIXEL_TO_UM
    tree = cKDTree(coords)
    edges = tree.query_pairs(threshold, output_type="ndarray")
    if len(edges) == 0:
        raise RuntimeError(f"No graph edges for {section}")
    row = np.r_[edges[:, 0], edges[:, 1]]
    col = np.r_[edges[:, 1], edges[:, 0]]
    graph = csr_matrix((np.ones(len(row), dtype=np.float32), (row, col)), shape=(len(labels), len(labels)))
    degrees = np.asarray(graph.sum(axis=1)).ravel()

    reduced_states = ["Other", TARGET] + candidates
    state_to_code = {state: index for index, state in enumerate(reduced_states)}
    codes = np.zeros(len(labels), dtype=np.int16)
    codes[labels == TARGET] = state_to_code[TARGET]
    for state in candidates:
        codes[labels == state] = state_to_code[state]
    candidate_codes = np.array([state_to_code[state] for state in candidates], dtype=int)
    observed = score_from_labels(graph, degrees, codes, state_to_code[TARGET], candidate_codes)

    rng = np.random.default_rng(SEED + section_index)
    null = np.empty((N_PERMUTATIONS, len(candidates)), dtype=float)
    for perm_idx in range(N_PERMUTATIONS):
        shuffled = rng.permutation(codes)
        null[perm_idx] = score_from_labels(
            graph, degrees, shuffled, state_to_code[TARGET], candidate_codes
        )
    null_mean = np.nanmean(null, axis=0)
    null_sd = np.nanstd(null, axis=0, ddof=1)
    z_scores = np.divide(
        observed - null_mean,
        null_sd,
        out=np.zeros_like(observed),
        where=null_sd > 0,
    )
    p_values = (1 + np.sum(null >= observed[None, :], axis=0)) / (N_PERMUTATIONS + 1)
    fdr = bh_adjust(p_values)

    rows: list[dict] = []
    for idx, state in enumerate(candidates):
        rows.append(
            {
                "section": section,
                "sample_id": sample_id,
                "timepoint": timepoint,
                "time_order": time_order,
                "conditioned_on": TARGET,
                "cell_state": state,
                "n_bins": len(labels),
                "n_Fib_K14": int(np.sum(labels == TARGET)),
                "n_cell_state": int(np.sum(labels == state)),
                "closest_threshold_coordinate_units": threshold,
                "fixed_radius_um": FIXED_RADIUS_UM,
                "squidpy_threshold_min": threshold_min,
                "squidpy_threshold_max": threshold_max,
                "n_undirected_edges": len(edges),
                "cooccurrence_score": observed[idx],
                "null_mean": null_mean[idx],
                "null_sd": null_sd[idx],
                "z_score": z_scores[idx],
                "p_enrichment": p_values[idx],
                "fdr_enrichment_within_section": fdr[idx],
                "n_permutations": N_PERMUTATIONS,
            }
        )
    print(section, len(labels), len(edges), threshold, flush=True)
    return rows


def main() -> None:
    with h5py.File(H5AD, "r") as handle:
        obs = handle["obs"]
        batch_group = obs["sample_batch_new"]
        label_group = obs["sub_labels"]
        batch_categories = np.asarray(batch_group["categories"]).astype("U")
        label_categories = np.asarray(label_group["categories"]).astype("U")
        all_batches = batch_categories[np.asarray(batch_group["codes"], dtype=int)]
        all_labels = label_categories[np.asarray(label_group["codes"], dtype=int)]
        all_coords = np.asarray(handle["obsm"]["spatial"])

    candidates = sorted(
        state for state in np.unique(all_labels) if state.startswith("KC_") or state.startswith("SAC_")
    )
    rows: list[dict] = []
    for index, (section, (sample_id, timepoint, time_order)) in enumerate(SAMPLES.items()):
        selected = all_batches == sample_id
        rows.extend(
            analyze_section(
                section,
                sample_id,
                timepoint,
                time_order,
                all_coords[selected],
                all_labels[selected],
                candidates,
                index,
            )
        )
    section_stats = pd.DataFrame(rows)
    section_stats.to_csv(OUT / "FigureS3_FibK14_KC_SAC_section_scores.csv", index=False)

    summaries: list[dict] = []
    for (timepoint, time_order, state), group in section_stats.groupby(
        ["timepoint", "time_order", "cell_state"], sort=False
    ):
        scores = group["cooccurrence_score"].to_numpy(dtype=float)
        z_values = group["z_score"].to_numpy(dtype=float)
        meta_z = float(np.sum(z_values) / np.sqrt(len(z_values)))
        summaries.append(
            {
                "timepoint": timepoint,
                "time_order": time_order,
                "cell_state": state,
                "n_sections": len(group),
                "mean_score": float(np.mean(scores)),
                "sd_score": float(np.std(scores, ddof=1)) if len(scores) > 1 else np.nan,
                "sem_score": float(np.std(scores, ddof=1) / np.sqrt(len(scores))) if len(scores) > 1 else np.nan,
                "min_score": float(np.min(scores)),
                "max_score": float(np.max(scores)),
                "meta_z": meta_z,
                "meta_p_enrichment": float(norm.sf(meta_z)),
            }
        )
    time_stats = pd.DataFrame(summaries)
    time_stats["meta_fdr"] = bh_adjust(time_stats["meta_p_enrichment"].to_numpy())
    time_stats["rank_within_time"] = time_stats.groupby("timepoint")["mean_score"].rank(
        method="first", ascending=False
    )
    time_stats.to_csv(OUT / "FigureS3_FibK14_KC_SAC_timepoint_summary.csv", index=False)
    time_stats[time_stats["rank_within_time"] <= 10].sort_values(
        ["time_order", "rank_within_time"]
    ).to_csv(OUT / "FigureS3_FibK14_KC_SAC_timepoint_TOP10.csv", index=False)


if __name__ == "__main__":
    main()
