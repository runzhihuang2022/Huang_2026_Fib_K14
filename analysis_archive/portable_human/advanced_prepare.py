#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


GROUPS = ["Normal", "SPTDI", "DPTDI"]


def write_sparse_csc(path: Path, matrix: sparse.csc_matrix, genes, cells) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        handle.attrs["format"] = "csc_genes_by_cells_v1"
        handle.create_dataset("shape", data=np.asarray(matrix.shape, dtype=np.int64))
        handle.create_dataset(
            "data", data=matrix.data.astype(np.float32, copy=False),
            compression="gzip", compression_opts=4, shuffle=True,
        )
        handle.create_dataset(
            "indices", data=matrix.indices.astype(np.int32, copy=False),
            compression="gzip", compression_opts=4, shuffle=True,
        )
        handle.create_dataset(
            "indptr", data=matrix.indptr.astype(np.int64, copy=False),
            compression="gzip", compression_opts=4, shuffle=True,
        )
        handle.create_dataset(
            "genes", data=np.asarray(genes, dtype=object), dtype=string_dtype,
        )
        handle.create_dataset(
            "cells", data=np.asarray(cells, dtype=object), dtype=string_dtype,
        )


def export_subset(
    adata: ad.AnnData, index: np.ndarray, path: Path, metadata_path: Path
) -> dict:
    matrix = adata.X[index].tocsr().transpose().tocsc()
    if np.any(matrix.data < 0) or not np.allclose(matrix.data, np.rint(matrix.data)):
        raise RuntimeError("scRNA module input is not nonnegative count-equivalent.")
    write_sparse_csc(
        path, matrix, adata.var_names.astype(str),
        adata.obs_names[index].astype(str),
    )
    columns = [
        "sample_id", "patient", "group", "bulk_labels", "sub_labels",
    ]
    metadata = adata.obs.iloc[index][columns].astype(str).copy()
    metadata.insert(0, "cell", adata.obs_names[index].astype(str))
    metadata.to_csv(metadata_path, sep="\t", index=False)
    return {
        "cells": int(index.size), "genes": int(adata.n_vars),
        "nnz": int(matrix.nnz), "matrix": str(path), "metadata": str(metadata_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.canonical)
    if adata.shape != (279305, 52399):
        raise RuntimeError(f"Frozen canonical shape violated: {adata.shape}")
    inventory = {"status": "PASS", "cellchat": {}, "monocle3": {}}
    groups = adata.obs["group"].astype(str).to_numpy()
    for group in GROUPS:
        index = np.flatnonzero(groups == group)
        inventory["cellchat"][group] = export_subset(
            adata, index, out / f"cellchat_{group}_counts.h5",
            out / f"cellchat_{group}_metadata.tsv",
        )
    bulk = adata.obs["bulk_labels"].astype(str).to_numpy()
    epithelial_index = np.flatnonzero(bulk == "Keratinocyte")
    inventory["monocle3"]["epithelial"] = export_subset(
        adata, epithelial_index, out / "monocle3_epithelial_counts.h5",
        out / "monocle3_epithelial_metadata.tsv",
    )
    umap = np.asarray(adata.obsm["X_umap"])[epithelial_index, :2]
    pd.DataFrame({
        "cell": adata.obs_names[epithelial_index].astype(str),
        "UMAP_1": umap[:, 0], "UMAP_2": umap[:, 1],
    }).to_csv(out / "monocle3_epithelial_umap.tsv", sep="\t", index=False)
    inventory["prohibited_operations_run"] = []
    (out / "advanced_prepare_checkpoint.json").write_text(
        json.dumps(inventory, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
