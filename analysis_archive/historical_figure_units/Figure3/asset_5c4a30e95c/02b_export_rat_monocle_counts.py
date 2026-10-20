"""Export the balanced rat trajectory count matrix for Monocle2."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite


STATES = [
    "Fib_K14", "SAC_Progenitor", "KC_Basal", "KC_Basal_Mig",
    "KC_Spinous_Mig", "KC_Spinous",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scrna", required=True)
    parser.add_argument("--package", required=True)
    args = parser.parse_args()
    tables = Path(args.package) / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.scrna, backed="r")
    labels = adata.obs["sub_labels"].astype(str)
    rng = np.random.default_rng(17)
    selected: list[int] = []
    for state in STATES:
        idx = np.flatnonzero(labels.eq(state).to_numpy())
        selected.extend(rng.choice(idx, min(1000, len(idx)), replace=False))
    traj = adata[np.asarray(selected)].to_memory()
    sc.pp.filter_genes(traj, min_cells=20)
    counts = traj.X.tocsr() if sparse.issparse(traj.X) else sparse.csr_matrix(traj.X)
    sparse.save_npz(str(tables / "S5E_monocle2_counts_cells_by_genes.npz"), counts)
    mmwrite(str(tables / "S5E_monocle2_counts_cells_by_genes.mtx"), counts)
    pd.DataFrame({"gene": traj.var_names.astype(str)}).to_csv(
        tables / "S5E_monocle2_genes.csv", index=False
    )
    pheno = traj.obs.copy()
    pheno["cell_id"] = traj.obs_names.astype(str)
    pheno.to_csv(tables / "S5E_monocle2_phenodata.csv", index=False)
    print(f"Exported {counts.shape[0]:,} cells x {counts.shape[1]:,} genes.")


if __name__ == "__main__":
    main()
