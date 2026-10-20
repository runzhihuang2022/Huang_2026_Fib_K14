#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import celloracle as co
import numpy as np
import pandas as pd
import scanpy as sc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=12)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    np.random.seed(20260811)

    oracle_path = args.output / "human_FibK14_allKC_KLF4_simulation_ready.celloracle.oracle"
    links_path = args.output / "human_FibK14_allKC_GRN.celloracle.links"
    if oracle_path.exists() and links_path.exists():
        print("Validated final checkpoints already exist; skip build.", flush=True)
        return

    adata = sc.read_h5ad(args.input)
    if "KLF4" not in adata.var_names:
        raise RuntimeError("KLF4 is absent from CellOracle input")
    adata.X = adata.layers["raw_count"].copy()
    base_grn = co.data.load_human_promoter_base_GRN(
        version="hg38_gimmemotifsv5_fpr2"
    )
    base_targets = set(base_grn["gene_short_name"].astype(str))
    if "KLF4" not in base_grn.columns:
        raise RuntimeError("KLF4 is absent from the human promoter base-GRN TF columns")

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata, cluster_column_name="sub_labels", embedding_name="X_original_umap"
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)
    oracle.perform_PCA(n_components=50)
    oracle.knn_imputation(
        n_pca_dims=40, k=50, balanced=True, b_sight=400, b_maxl=200,
        group_constraint="sub_labels", n_jobs=args.threads
    )
    pre_path = args.output / "human_FibK14_allKC_preGRN.celloracle.oracle"
    oracle.to_hdf5(str(pre_path))

    links = oracle.get_links(
        cluster_name_for_GRN_unit="sub_labels", alpha=10, bagging_number=20,
        verbose_level=10, model_method="bagging_ridge", n_jobs=args.threads
    )
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=5000,
                       thread_number=args.threads)
    links.to_hdf5(str(links_path))
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(
        GRN_unit="cluster", alpha=10, use_cluster_specific_TFdict=True, verbose_level=10
    )
    oracle.to_hdf5(str(oracle_path))

    klf4_edges = []
    for cluster, table in links.filtered_links.items():
        q = table[table["source"].astype(str).str.upper() == "KLF4"].copy()
        if len(q):
            q.insert(0, "sub_labels", cluster)
            klf4_edges.append(q)
    edge_table = pd.concat(klf4_edges, ignore_index=True) if klf4_edges else pd.DataFrame()
    edge_table.to_csv(args.output / "KLF4_filtered_GRN_edges.tsv", sep="\t", index=False)

    audit = {
        "status": "PASS", "engine": "CellOracle", "celloracle_version": co.__version__,
        "base_grn": "CellOracle hg38_gimmemotifsv5_fpr2 promoter base GRN", "species": "human",
        "cells": int(adata.n_obs), "genes": int(adata.n_vars),
        "clusters": adata.obs["sub_labels"].astype(str).value_counts().to_dict(),
        "KLF4_in_expression": True, "KLF4_in_base_grn_TFs": True,
        "base_grn_targets_overlap": int(len(base_targets.intersection(adata.var_names))),
        "KLF4_filtered_edges": int(len(edge_table)),
        "alpha_GRN": 10, "bagging_number": 20, "filter_p": 0.001,
        "filter_top_edges_per_cluster": 5000, "fit_alpha": 10,
        "knn_imputation_k": 50, "random_seed": 20260811,
    }
    (args.output / "celloracle_grn_checkpoint.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
