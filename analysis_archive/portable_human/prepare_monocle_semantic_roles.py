#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from advanced_prepare import export_subset


TRAJECTORY_STATES = ["KC_Basal", "KC_Basal_Mig", "KC_Spinous_Mig"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.canonical)
    labels = adata.obs["sub_labels"].astype(str).to_numpy()
    missing = [x for x in TRAJECTORY_STATES if x not in set(labels)]
    if missing:
        raise RuntimeError(f"Missing declared trajectory states: {missing}")
    index = np.flatnonzero(np.isin(labels, TRAJECTORY_STATES))
    inventory = export_subset(
        adata, index, out / "monocle3_epithelial_counts.h5",
        out / "monocle3_epithelial_metadata.tsv",
    )
    umap = np.asarray(adata.obsm["X_umap"])[index, :2]
    pd.DataFrame({
        "cell": adata.obs_names[index].astype(str),
        "UMAP_1": umap[:, 0], "UMAP_2": umap[:, 1],
    }).to_csv(out / "monocle3_epithelial_umap.tsv", sep="\t", index=False)
    state_counts = pd.Series(labels[index]).value_counts().to_dict()
    checkpoint = {
        "status": "PASS",
        "selection": "exact_manifest_semantic_roles",
        "states": TRAJECTORY_STATES,
        "state_counts": {k: int(v) for k, v in state_counts.items()},
        "inventory": inventory,
        "random_sampling": False,
        "global_reclustering": False,
        "global_embedding_recomputed": False,
    }
    (out / "monocle3_semantic_role_prepare_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
