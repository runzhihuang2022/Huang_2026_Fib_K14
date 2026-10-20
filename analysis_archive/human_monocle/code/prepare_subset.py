#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


FOCUS_GENES = ["KLF4", "KRT5", "KRT14", "TP63", "EPCAM", "TACSTD2", "PDGFRA", "VIM", "COL1A1"]


def write_sparse_h5(matrix, genes, cells, path: Path) -> None:
    # R Matrix::sparseMatrix(repr="C") expects a column-compressed matrix.
    matrix = sparse.csc_matrix(matrix)
    with h5py.File(path, "w") as h5:
        h5.create_dataset("data", data=matrix.data, compression="gzip")
        h5.create_dataset("indices", data=matrix.indices, compression="gzip")
        h5.create_dataset("indptr", data=matrix.indptr, compression="gzip")
        h5.create_dataset("shape", data=np.asarray(matrix.shape, dtype=np.int64))
        h5.create_dataset("genes", data=np.asarray(genes, dtype="S"))
        h5.create_dataset("cells", data=np.asarray(cells, dtype="S"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=20260811)
    ap.add_argument("--hvg-cells", type=int, default=50000)
    ap.add_argument("--n-hvg", type=int, default=3000)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    required_cols = ["bulk_labels", "sub_labels", "patient", "Time_category", "sample_id", "orig.ident"]
    missing_cols = [c for c in required_cols if c not in adata.obs]
    if missing_cols:
        raise RuntimeError(f"Missing metadata columns: {missing_cols}")

    bulk = adata.obs["bulk_labels"].astype(str)
    sub = adata.obs["sub_labels"].astype(str)
    keep = (bulk == "Keratinocyte") | (sub == "Fib_K14")
    subset = adata[keep].copy()
    if subset.n_obs == 0:
        raise RuntimeError("Fib_K14 + Keratinocyte gate returned zero cells")
    if int((subset.obs["sub_labels"].astype(str) == "Fib_K14").sum()) == 0:
        raise RuntimeError("Fib_K14 is absent after gating")

    rng = np.random.default_rng(args.seed)
    labels = subset.obs["sub_labels"].astype(str)
    chosen = []
    for _, positions in pd.Series(np.arange(subset.n_obs), index=subset.obs_names).groupby(labels.values):
        arr = positions.to_numpy(dtype=int)
        quota = max(250, int(round(args.hvg_cells * len(arr) / subset.n_obs)))
        chosen.extend(rng.choice(arr, size=min(quota, len(arr)), replace=False).tolist())
    chosen = np.asarray(sorted(set(chosen)), dtype=int)
    if len(chosen) > args.hvg_cells:
        chosen = np.sort(rng.choice(chosen, size=args.hvg_cells, replace=False))

    hvg = subset[chosen].copy()
    sc.pp.filter_genes(hvg, min_cells=20)
    sc.pp.normalize_total(hvg, target_sum=1e4)
    sc.pp.log1p(hvg)
    sc.pp.highly_variable_genes(hvg, flavor="seurat", n_top_genes=args.n_hvg)
    selected = hvg.var_names[hvg.var["highly_variable"].to_numpy()].astype(str).tolist()
    present_focus = [g for g in FOCUS_GENES if g in subset.var_names]
    selected = list(dict.fromkeys(selected + present_focus))

    counts = sparse.csr_matrix(subset[:, selected].X).T.tocsr()
    cells = subset.obs_names.astype(str).tolist()
    write_sparse_h5(counts, selected, cells, args.output / "fibk14_kc_counts.h5")

    meta = subset.obs[required_cols].copy()
    meta.insert(0, "cell", cells)
    meta.to_csv(args.output / "fibk14_kc_metadata.tsv", sep="\t", index=False)
    subtype_counts = labels.value_counts().sort_index()
    subtype_counts.rename_axis("sub_labels").rename("n_cells").to_csv(
        args.output / "subset_cell_counts.tsv", sep="\t"
    )
    audit = {
        "status": "PASS",
        "species": "human",
        "genome_build": "GRCh38",
        "input": str(args.input),
        "gate": "bulk_labels == Keratinocyte OR sub_labels == Fib_K14",
        "n_cells": int(subset.n_obs),
        "n_keratinocytes": int((subset.obs["bulk_labels"].astype(str) == "Keratinocyte").sum()),
        "n_fib_k14": int((subset.obs["sub_labels"].astype(str) == "Fib_K14").sum()),
        "n_hvg_selection_cells": int(len(chosen)),
        "n_features_exported": int(len(selected)),
        "focus_gene_presence": {g: bool(g in subset.var_names) for g in FOCUS_GENES},
        "all_cells_retained_for_trajectory": True,
        "stored_global_umap_reused": False,
        "random_seed": args.seed,
    }
    (args.output / "subset_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
