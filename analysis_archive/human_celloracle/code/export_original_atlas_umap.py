#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    adata = ad.read_h5ad(args.input, backed="r")
    xy = np.asarray(adata.obsm["X_umap"])
    out = pd.DataFrame({
        "cell": adata.obs_names.astype(str),
        "sub_labels": adata.obs["sub_labels"].astype(str).to_numpy(),
        "bulk_labels": adata.obs["bulk_labels"].astype(str).to_numpy(),
        "UMAP_1": xy[:,0], "UMAP_2": xy[:,1],
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
