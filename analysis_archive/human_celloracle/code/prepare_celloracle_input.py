#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


FOCUS = ["KLF4", "KRT5", "KRT14", "TP63", "EPCAM", "TACSTD2", "PDGFRA", "VIM", "COL1A1"]
PALETTE = {
    "KC_Basal": "#008941", "KC_Basal_Mig": "#FF3030", "KC_Basal_Prolif": "#FFD700",
    "KC_Spinous": "#00FF00", "KC_Spinous_Mig": "#1A0099", "KC_Spinous_Mat": "#EE4C97",
    "KC_Granular": "#882D17", "Fib_K14": "#4682B4",
}


def stratified_indices(labels: pd.Series, total: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    pos = np.arange(len(labels))
    fib = pos[labels.to_numpy() == "Fib_K14"]
    selected = fib.tolist()
    remaining = max(0, total - len(selected))
    kc_labels = labels[labels != "Fib_K14"]
    counts = kc_labels.value_counts()
    for label, count in counts.items():
        idx = pos[labels.to_numpy() == label]
        quota = max(500, int(round(remaining * count / counts.sum())))
        selected.extend(rng.choice(idx, min(quota, len(idx)), replace=False).tolist())
    selected = np.asarray(sorted(set(selected)), dtype=int)
    if len(selected) > total:
        fib_set = set(fib.tolist())
        other = np.asarray([x for x in selected if x not in fib_set], dtype=int)
        other = rng.choice(other, total - len(fib), replace=False)
        selected = np.sort(np.r_[fib, other])
    return selected


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--pseudotime", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--n-cells", type=int, default=30000)
    ap.add_argument("--n-hvg", type=int, default=2500)
    ap.add_argument("--seed", type=int, default=20260811)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    labels = adata.obs["sub_labels"].astype(str)
    keep = (adata.obs["bulk_labels"].astype(str) == "Keratinocyte") | (labels == "Fib_K14")
    subset = adata[keep].copy()
    labels = subset.obs["sub_labels"].astype(str)
    idx = stratified_indices(labels, args.n_cells, args.seed)
    model = subset[idx].copy()
    model.obs["sub_labels"] = model.obs["sub_labels"].astype(str).astype("category")
    model.uns["sub_labels_colors"] = [PALETTE[x] for x in model.obs["sub_labels"].cat.categories]

    pt = pd.read_csv(args.pseudotime, sep="\t", compression="gzip").set_index("cell")
    model.obs["Pseudotime"] = pt.reindex(model.obs_names)["pseudotime"].to_numpy()
    if model.obs["Pseudotime"].isna().any():
        raise RuntimeError("Pseudotime join produced missing values")
    model.obsm["X_original_umap"] = np.asarray(model.obsm["X_umap"]).copy()

    hvg = model.copy()
    sc.pp.filter_genes(hvg, min_cells=20)
    sc.pp.normalize_total(hvg, target_sum=1e4)
    sc.pp.log1p(hvg)
    sc.pp.highly_variable_genes(hvg, flavor="seurat", n_top_genes=args.n_hvg)
    genes = hvg.var_names[hvg.var["highly_variable"]].astype(str).tolist()
    genes = list(dict.fromkeys(genes + [g for g in FOCUS if g in model.var_names]))
    model = model[:, genes].copy()
    model.layers["raw_count"] = model.X.copy()
    model.write_h5ad(args.output / "human_FibK14_allKC_celloracle_input.h5ad", compression="gzip")

    model.obs["sub_labels"].value_counts().rename("n_cells").to_csv(
        args.output / "model_cell_counts.tsv", sep="\t"
    )
    audit = {
        "status": "PASS", "species": "human", "genome_build": "GRCh38",
        "gate": "Fib_K14 plus all bulk Keratinocyte states",
        "full_subset_cells": int(subset.n_obs), "model_cells": int(model.n_obs),
        "model_genes": int(model.n_vars), "all_fib_k14_retained": True,
        "celloracle_downsampling": "all Fib_K14 plus deterministic stratified KC sample",
        "embedding": "original parent X_umap copied to X_original_umap",
        "base_grn_planned": "CellOracle human promoter base GRN",
        "random_seed": args.seed,
    }
    (args.output / "prepare_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
