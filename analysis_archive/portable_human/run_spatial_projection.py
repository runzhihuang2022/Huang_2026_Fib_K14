#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROTATIONS = {
    "BW32_A01597A3_SDSDB_5dpb": 135,
    "BW13_1_B1_DSDB_12dpb": 45,
    "BW14_1_C1_SSDB_12dpb": 135,
    "BW13_A3_DSDB_12dpb": 135,
    "BW14_B3_SSDB_12dpb": 135,
    "BW15D_C6_SDSDB_19dpb": -225,
    "BW15D_1_D1_SDSDB_19dpb": 140,
    "BW81_C02846B6_SDSDB_26dpb_part1": 0,
    "BW81_C02846B6_SDSDB_26dpb_part2": 0,
    "2mph_A03699G6.SCT": 45,
    "NS_C02847B1": 180,
}


def normalize_counts(x: sparse.csr_matrix) -> sparse.csr_matrix:
    library = np.asarray(x.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise RuntimeError("Zero-library scRNA cells in frozen canonical object.")
    result = x.multiply((10000.0 / library)[:, None]).tocsr()
    result.data = np.log1p(result.data).astype(np.float32, copy=False)
    return result


def reference_centroids(
    normalized: sparse.csr_matrix, labels: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    levels = np.asarray(sorted(pd.unique(labels)))
    codes = pd.Categorical(labels, categories=levels).codes
    indicator = sparse.csr_matrix(
        (
            np.ones(labels.size, dtype=np.float32),
            (codes, np.arange(labels.size, dtype=np.int64)),
        ),
        shape=(levels.size, labels.size),
    )
    centroids = (indicator @ normalized).toarray()
    sizes = np.bincount(codes, minlength=levels.size).astype(np.float32)
    centroids /= sizes[:, None]
    return levels, centroids.astype(np.float32)


def choose_features(
    levels: np.ndarray,
    centroids: np.ndarray,
    sc_genes: np.ndarray,
    stereo_gene_set: set[str],
    top_n: int = 50,
) -> np.ndarray:
    feature_index: set[int] = set()
    for index, _ in enumerate(levels):
        other = np.delete(centroids, index, axis=0)
        score = centroids[index] - np.mean(other, axis=0)
        order = np.argsort(score)[::-1]
        selected = [
            int(i) for i in order
            if score[i] > 0 and sc_genes[i] in stereo_gene_set
        ][:top_n]
        feature_index.update(selected)
    return np.asarray(sorted(feature_index), dtype=np.int64)


def project(
    stereo_x: sparse.csr_matrix,
    stereo_var_index: dict[str, int],
    sc_genes: np.ndarray,
    feature_index: np.ndarray,
    levels: np.ndarray,
    centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    genes = sc_genes[feature_index]
    query_columns = np.asarray([stereo_var_index[x] for x in genes], dtype=np.int64)
    query = stereo_x[:, query_columns].tocsr().astype(np.float32)
    reference = centroids[:, feature_index].astype(np.float32)
    reference_norm = np.linalg.norm(reference, axis=1)
    reference /= np.maximum(reference_norm[:, None], 1e-12)
    query_norm = np.sqrt(np.asarray(query.multiply(query).sum(axis=1)).ravel())
    query = query.multiply((1.0 / np.maximum(query_norm, 1e-12))[:, None])
    scores = np.asarray(query @ reference.T, dtype=np.float32)
    best = np.argmax(scores, axis=1)
    best_score = scores[np.arange(scores.shape[0]), best]
    scores[np.arange(scores.shape[0]), best] = -np.inf
    second = np.max(scores, axis=1)
    return levels[best], best_score, best_score - second


def write_frame(frame: pd.DataFrame, path: Path) -> str:
    try:
        target = path.with_suffix(".parquet")
        frame.to_parquet(target, index=False)
    except Exception:
        target = path.with_suffix(".tsv.gz")
        frame.to_csv(target, sep="\t", index=False, compression="gzip")
    return str(target)


def rotated_coordinates(
    sample_ids: np.ndarray, spatial: np.ndarray
) -> pd.DataFrame:
    output = pd.DataFrame({
        "sample_batch_new": sample_ids,
        "spatial_x_source": spatial[:, 0],
        "spatial_y_source": spatial[:, 1],
        "rotation_degrees_declared": [ROTATIONS[x] for x in sample_ids],
        "spatial_x_rotated": np.nan,
        "spatial_y_rotated": np.nan,
    })
    for sample, angle in ROTATIONS.items():
        index = np.flatnonzero(sample_ids == sample)
        center = spatial[index].mean(axis=0)
        radians = np.deg2rad(angle)
        rotation = np.asarray([
            [np.cos(radians), -np.sin(radians)],
            [np.sin(radians), np.cos(radians)],
        ])
        transformed = (spatial[index] - center) @ rotation.T + center
        output.loc[index, "spatial_x_rotated"] = transformed[:, 0]
        output.loc[index, "spatial_y_rotated"] = transformed[:, 1]
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scrna", required=True, type=Path)
    parser.add_argument("--stereo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    scrna = ad.read_h5ad(args.scrna)
    stereo = ad.read_h5ad(args.stereo)
    if scrna.shape != (279305, 52399):
        raise RuntimeError(f"Frozen scRNA shape violated: {scrna.shape}")
    if stereo.shape[0] != 387061:
        raise RuntimeError(f"Frozen Stereo subset bin count violated: {stereo.shape}")
    samples = stereo.obs["sample_batch_new"].astype(str).to_numpy()
    if set(samples) != set(ROTATIONS):
        raise RuntimeError("Stereo sample scope differs from the frozen 11-section set.")
    if "spatial" not in stereo.obsm:
        raise RuntimeError("Stereo obsm['spatial'] missing.")

    sc_x = scrna.X.tocsr()
    if np.any(sc_x.data < 0) or not np.allclose(sc_x.data, np.rint(sc_x.data)):
        raise RuntimeError("scRNA X is not count-equivalent.")
    normalized = normalize_counts(sc_x)
    sc_genes = scrna.var_names.astype(str).to_numpy()
    stereo_genes = stereo.var_names.astype(str).to_numpy()
    stereo_gene_index = {gene: index for index, gene in enumerate(stereo_genes)}
    stereo_gene_set = set(stereo_gene_index)
    stereo_x = stereo.X.tocsr().astype(np.float32)
    result = pd.DataFrame({
        "bin": stereo.obs_names.astype(str),
        "sample_batch_new": samples,
    })
    inventory = {}
    for track in ["bulk_labels", "sub_labels"]:
        labels = scrna.obs[track].astype(str).to_numpy()
        levels, centroids = reference_centroids(normalized, labels)
        features = choose_features(
            levels, centroids, sc_genes, stereo_gene_set, top_n=50
        )
        predicted, score, margin = project(
            stereo_x, stereo_gene_index, sc_genes, features, levels, centroids
        )
        result[f"predicted_{track}"] = predicted
        result[f"{track}_cosine_score"] = score
        result[f"{track}_score_margin"] = margin
        pd.DataFrame({
            "gene": sc_genes[features],
        }).to_csv(out / f"M12_{track}_projection_features.tsv", sep="\t", index=False)
        inventory[track] = {
            "reference_labels": int(levels.size),
            "projection_features": int(features.size),
            "reference_cells": int(scrna.n_obs),
            "query_bins": int(stereo.n_obs),
        }
    coordinates = rotated_coordinates(samples, np.asarray(stereo.obsm["spatial"])[:, :2])
    result = pd.concat([result, coordinates.drop(columns="sample_batch_new")], axis=1)
    projection_path = write_frame(result, out / "M12_stereoseq_projection")
    for track in ["bulk_labels", "sub_labels"]:
        counts = (
            result.groupby(["sample_batch_new", f"predicted_{track}"])
            .size().rename("bins").reset_index()
        )
        totals = counts.groupby("sample_batch_new")["bins"].transform("sum")
        counts["proportion"] = counts["bins"] / totals
        counts.to_csv(
            out / f"M12_{track}_sample_composition.tsv", sep="\t", index=False
        )
    checkpoint = {
        "module": "M12",
        "status": "PASS",
        "evidence_status": "faithful reconstruction",
        "method": "current-project logCP10K reference-centroid cosine projection",
        "old_csc_results_used": False,
        "scrna_shape": list(map(int, scrna.shape)),
        "stereo_shape": list(map(int, stereo.shape)),
        "shared_genes": int(len(set(sc_genes) & stereo_gene_set)),
        "tracks": inventory,
        "projection_table": projection_path,
        "coordinate_source": "obsm['spatial']",
        "declared_rotations_applied_about_each_section_centroid": True,
        "he_alignment_claimed": False,
        "epidermis_up_claimed": False,
        "physical_scale_claimed": False,
        "visium_used": False,
        "prohibited_global_recomputation": False,
    }
    (out / "M12_spatial_projection_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
