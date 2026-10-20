from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
BASE = Path(r"data/external/local_8e645d0a9b03")
H5 = BASE / "Wound_Healing_Annotation_Output_White" / "BW15D_1_D1_SDSDB_19dpb_Annotated.h5ad"
ANCHORS = BASE / "mask" / "19dpb_p1_spatial_anchors.json"
REGIONS = ["Un-epi", "Epi-Front", "Newly-epi"]
TARGETS = ["KC_Spinous_Mig", "SAC_SG_Progenitor", "Fib_K14"]


def distance_to_polyline(points: np.ndarray, line: list[list[float]]) -> np.ndarray:
    vertices = np.asarray(line, dtype=float)
    result = np.full(points.shape[0], np.inf)
    for start, end in zip(vertices[:-1], vertices[1:]):
        vector = end - start
        denominator = float(np.dot(vector, vector))
        if denominator <= np.finfo(float).eps:
            continue
        delta = points - start
        projection = np.clip((delta[:, 0] * vector[0] + delta[:, 1] * vector[1]) / denominator, 0, 1)
        nearest = start + projection[:, None] * vector
        result = np.minimum(result, np.sqrt(((points - nearest) ** 2).sum(axis=1)))
    return result


def points_in_polygon(points: np.ndarray, vertices: list[list[float]]) -> np.ndarray:
    polygon = np.asarray(vertices, dtype=float)
    inside = np.zeros(points.shape[0], dtype=bool)
    x, y = points[:, 0], points[:, 1]
    j = len(polygon) - 1
    for i in range(len(polygon)):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        crosses = ((yi > y) != (yj > y)) & (
            x < (xj - xi) * (y - yi) / ((yj - yi) + np.finfo(float).eps) + xi
        )
        inside ^= crosses
        j = i
    return inside


def decode_categorical(group: h5py.Group) -> np.ndarray:
    categories = np.asarray(group["categories"]).astype(str)
    return categories[np.asarray(group["codes"], dtype=int)]


def region_assignment(coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    anchors = json.loads(ANCHORS.read_text(encoding="utf-8"))
    theta = np.radians(140)
    rotated = np.column_stack((
        coordinates[:, 0] * np.cos(theta) - coordinates[:, 1] * np.sin(theta),
        coordinates[:, 0] * np.sin(theta) + coordinates[:, 1] * np.cos(theta),
    ))
    epidermal_distance = distance_to_polyline(rotated, anchors["epi_baseline"])
    dermal_distance = distance_to_polyline(rotated, anchors["der_bottom"])
    depth = np.divide(epidermal_distance, epidermal_distance + dermal_distance,
                      out=np.zeros_like(epidermal_distance),
                      where=(epidermal_distance + dermal_distance) > 0)
    front = np.asarray(anchors["leading_edge"], dtype=float)
    distance = distance_to_polyline(rotated, anchors["leading_edge"])
    right_is_healed = anchors.get("healed_direction") == "right_is_healed"
    positive = ((rotated[:, 0] > front[0, 0]) & right_is_healed) | ((rotated[:, 0] < front[0, 0]) & (not right_is_healed))
    horizontal = distance * np.where(positive, 1, -1) * 0.33
    max_abs = max(abs(horizontal.min()), abs(horizontal.max())) * 1.1
    model_points = np.column_stack((horizontal, depth))
    polygons = {
        "Un-epi": [[-max_abs, 0], [-600, 0], [-300, 1], [-max_abs, 1]],
        "Epi-Front": [[-600, 0], [400, 0], [700, 1], [-300, 1]],
        "Newly-epi": [[400, 0], [max_abs, 0], [max_abs, 1], [700, 1]],
    }
    region = np.full(coordinates.shape[0], "", dtype=object)
    for name, polygon in polygons.items():
        region[points_in_polygon(model_points, polygon)] = name
    return region, rotated * 0.33, model_points


def main() -> None:
    with h5py.File(H5, "r") as handle:
        labels = decode_categorical(handle["obs"]["sub_labels"])
        coords = np.asarray(handle["obsm"]["spatial"], dtype=float)
        scores = {t: np.asarray(handle["obs"][f"{t}_MergeScore"], dtype=float) for t in TARGETS}
    region, coords_um, model_points = region_assignment(coords)
    if np.any(region == ""):
        raise RuntimeError(f"{np.sum(region == '')} bins were not assigned to a region")

    rows = []
    for cell_type in np.unique(labels):
        total_cell = np.sum(labels == cell_type)
        for reg in REGIONS:
            mask = region == reg
            count = int(np.sum(mask & (labels == cell_type)))
            n = int(mask.sum())
            prop = count / n
            global_prop = total_cell / len(labels)
            enrich = np.log2((prop + 0.5 / n) / (global_prop + 0.5 / len(labels)))
            table = [[count, n - count], [total_cell - count, (len(labels) - n) - (total_cell - count)]]
            _, p = fisher_exact(table, alternative="two-sided")
            rows.append([cell_type, reg, count, n, prop, enrich, p])
    enrichment = pd.DataFrame(rows, columns=["cell_type", "region", "count", "region_bins", "proportion", "log2_enrichment", "fisher_p"])
    enrichment["fisher_fdr"] = multipletests(enrichment["fisher_p"], method="fdr_bh")[1]
    enrichment.to_csv(OUT / "Figure2B_S2AB_19dpb_p1_region_enrichment.csv", index=False)

    fib = enrichment[enrichment.cell_type.eq("Fib_K14")].copy()
    fib.to_csv(OUT / "Figure2E_Fib_K14_region_proportions.csv", index=False)
    pairwise = []
    for i, a in enumerate(REGIONS):
        for b in REGIONS[i + 1:]:
            aa = fib[fib.region.eq(a)].iloc[0]
            bb = fib[fib.region.eq(b)].iloc[0]
            table = [[aa["count"], aa["region_bins"] - aa["count"]], [bb["count"], bb["region_bins"] - bb["count"]]]
            odds, p = fisher_exact(table, alternative="two-sided")
            pairwise.append([a, b, odds, p])
    pairwise = pd.DataFrame(pairwise, columns=["region_a", "region_b", "odds_ratio", "fisher_p"])
    pairwise["fisher_fdr"] = multipletests(pairwise["fisher_p"], method="fdr_bh")[1]
    pairwise.to_csv(OUT / "Figure2E_Fib_K14_pairwise_Fisher.csv", index=False)

    front = region == "Epi-Front"
    score_df = pd.DataFrame({"x_um": coords_um[front, 0], "y_um": coords_um[front, 1], **{t: scores[t][front] for t in TARGETS}})
    score_df.to_csv(OUT / "FigureS2CD_EpiFront_mapping_scores.csv", index=False)
    corr_rows = []
    for a, b, method in [("Fib_K14", "KC_Spinous_Mig", "Pearson"), ("Fib_K14", "SAC_SG_Progenitor", "Spearman")]:
        statistic, p = pearsonr(score_df[a], score_df[b]) if method == "Pearson" else spearmanr(score_df[a], score_df[b])
        corr_rows.append([a, b, method, statistic, p, len(score_df)])
    pd.DataFrame(corr_rows, columns=["state_a", "state_b", "method", "statistic", "p_value", "n_bins"]).to_csv(OUT / "FigureS2CD_correlations.csv", index=False)

    rng = np.random.default_rng(20260722)
    pairs = [("Fib_K14", "KC_Spinous_Mig"), ("Fib_K14", "SAC_SG_Progenitor"), ("KC_Spinous_Mig", "SAC_SG_Progenitor")]
    front_indices = np.where(front)[0]
    front_labels = labels[front]
    front_coords = coords_um[front]
    nn_rows = []
    permutation_rows = []
    for a, b in pairs:
        a_coords = front_coords[front_labels == a]
        b_coords = front_coords[front_labels == b]
        observed = float(cKDTree(b_coords).query(a_coords, k=1)[0].mean())
        perm_values = []
        for _ in range(1000):
            shuffled = rng.permutation(front_labels)
            pa = front_coords[shuffled == a]
            pb = front_coords[shuffled == b]
            perm_values.append(float(cKDTree(pb).query(pa, k=1)[0].mean()))
        perm_values = np.asarray(perm_values)
        p = (1 + np.sum(perm_values <= observed)) / (len(perm_values) + 1)
        nn_rows.append([a, b, len(a_coords), len(b_coords), observed, perm_values.mean(), np.quantile(perm_values, .025), np.quantile(perm_values, .975), p])
        permutation_rows.extend([[a, b, i + 1, v] for i, v in enumerate(perm_values)])
    pd.DataFrame(nn_rows, columns=["state_a", "state_b", "n_a", "n_b", "observed_mean_nn_um", "permuted_mean_um", "perm_q025_um", "perm_q975_um", "permutation_p"]).to_csv(OUT / "FigureS2E_nearest_neighbor_summary.csv", index=False)
    pd.DataFrame(permutation_rows, columns=["state_a", "state_b", "permutation", "mean_nn_um"]).to_csv(OUT / "FigureS2E_nearest_neighbor_permutations.csv", index=False)
    pd.DataFrame({"region": region, "model_x_um": model_points[:, 0], "relative_depth": model_points[:, 1], "sub_label": labels}).to_csv(OUT / "19dpb_p1_region_assignments.csv", index=False)


if __name__ == "__main__":
    main()
