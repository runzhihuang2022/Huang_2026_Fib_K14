from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np


BASE = Path(r"data/external/local_8e645d0a9b03")
CONFIG = [
    ("19dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("19dpb_p1", "BW15D_1_D1_SDSDB_19dpb", 140),
]


def distance_to_polyline(points: np.ndarray, line: list[list[float]]) -> np.ndarray:
    vertices = np.asarray(line, dtype=float)
    result = np.full(points.shape[0], np.inf)
    for start, end in zip(vertices[:-1], vertices[1:]):
        vector = end - start
        denominator = float(np.dot(vector, vector))
        delta = points - start
        projection = np.clip((delta[:, 0] * vector[0] + delta[:, 1] * vector[1]) / denominator, 0, 1)
        nearest = start + projection[:, None] * vector
        result = np.minimum(result, np.sqrt(((points - nearest) ** 2).sum(axis=1)))
    return result


def points_in_polygon(points: np.ndarray, vertices: list[list[float]]) -> np.ndarray:
    # Vectorized ray casting; boundary points are immaterial at the present bin scale.
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


def main() -> None:
    for key, sample, angle in CONFIG:
        object_path = BASE / "Wound_Healing_Annotation_Output_White" / f"{sample}_Annotated.h5ad"
        anchor_path = BASE / "mask" / f"{key}_spatial_anchors.json"
        print(f"Opening {sample}", flush=True)
        with h5py.File(object_path, "r") as handle:
            encoded = handle["obs"]["sub_labels"]
            categories = np.asarray(encoded["categories"]).astype(str)
            labels = categories[np.asarray(encoded["codes"], dtype=int)]
            coordinates = np.asarray(handle["obsm"]["spatial"], dtype=float)

        anchors = json.loads(anchor_path.read_text(encoding="utf-8"))
        theta = np.radians(angle)
        rotated = np.column_stack(
            (
                coordinates[:, 0] * np.cos(theta) - coordinates[:, 1] * np.sin(theta),
                coordinates[:, 0] * np.sin(theta) + coordinates[:, 1] * np.cos(theta),
            )
        )
        epidermal_distance = distance_to_polyline(rotated, anchors["epi_baseline"])
        dermal_distance = distance_to_polyline(rotated, anchors["der_bottom"])
        depth = np.divide(
            epidermal_distance,
            epidermal_distance + dermal_distance,
            out=np.zeros_like(epidermal_distance),
            where=(epidermal_distance + dermal_distance) > 0,
        )
        front = np.asarray(anchors["leading_edge"], dtype=float)
        distance = distance_to_polyline(rotated, anchors["leading_edge"])
        right_is_healed = anchors.get("healed_direction") == "right_is_healed"
        positive = ((rotated[:, 0] > front[0, 0]) & right_is_healed) | (
            (rotated[:, 0] < front[0, 0]) & (not right_is_healed)
        )
        horizontal = distance * np.where(positive, 1, -1) * 0.33
        max_abs = max(abs(horizontal.min()), abs(horizontal.max())) * 1.1
        points = np.column_stack((horizontal, depth))
        polygons = {
            "Un-epi": [[-max_abs, 0], [-600, 0], [-300, 1], [-max_abs, 1]],
            "Epi-Front": [[-600, 0], [400, 0], [700, 1], [-300, 1]],
            "Newly-epi": [[400, 0], [max_abs, 0], [max_abs, 1], [700, 1]],
        }
        print(f"\n{key}: {sample}; healed_direction={anchors.get('healed_direction')}")
        for region, polygon in polygons.items():
            selected = points_in_polygon(points, polygon)
            fib = labels[selected] == "Fib_K14"
            percentage = 100 * fib.mean() if selected.any() else np.nan
            print(region, int(selected.sum()), int(fib.sum()), f"{percentage:.3f}%")


if __name__ == "__main__":
    main()
