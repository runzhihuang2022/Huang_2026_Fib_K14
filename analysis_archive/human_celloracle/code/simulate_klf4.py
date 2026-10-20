#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import celloracle as co
import numpy as np
import pandas as pd
from scipy import sparse

# CellOracle 0.18.0 still accesses ``csr_matrix.A`` while recent SciPy
# releases removed that convenience property.  Restore the exact historical
# read-only behavior without changing the numerical calculation.
if not hasattr(sparse.csr_matrix, "A"):
    sparse.csr_matrix.A = property(lambda self: self.toarray())
if not hasattr(sparse.csc_matrix, "A"):
    sparse.csc_matrix.A = property(lambda self: self.toarray())


def gene_vector(oracle, gene: str) -> np.ndarray:
    idx = int(np.flatnonzero(oracle.adata.var_names.to_numpy() == gene)[0])
    x = oracle.adata.layers["imputed_count"][:, idx]
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def run_one(oracle_path: Path, condition: str, value: float, output: Path, threads: int):
    oracle = co.load_hdf5(str(oracle_path))
    oracle.simulate_shift(
        perturb_condition={"KLF4": float(value)}, GRN_unit="cluster",
        n_propagation=3, ignore_warning=False, clip_delta_X=True
    )
    oracle.estimate_transition_prob(
        n_neighbors=200, knn_random=True, sampled_fraction=0.30,
        n_jobs=threads, threads=threads, calculate_randomized=True,
        random_seed=20260811
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200, n_jobs=threads)
    oracle.calculate_mass_filter(min_mass=0.01, plot=False)

    coords = np.asarray(oracle.adata.obsm[oracle.embedding_name])
    delta = np.asarray(oracle.delta_embedding)
    out = pd.DataFrame({
        "cell": oracle.adata.obs_names.astype(str),
        "sub_labels": oracle.adata.obs[oracle.cluster_column_name].astype(str).to_numpy(),
        "pseudotime": oracle.adata.obs["Pseudotime"].to_numpy(),
        "UMAP_1": coords[:,0], "UMAP_2": coords[:,1],
        "delta_UMAP_1": delta[:,0], "delta_UMAP_2": delta[:,1],
    })
    out.to_csv(output / f"KLF4_{condition}_cell_vectors.tsv.gz", sep="\t", index=False)
    np.savez_compressed(
        output / f"KLF4_{condition}_celloracle_grid.npz",
        flow_grid=np.asarray(oracle.flow_grid), flow=np.asarray(oracle.flow),
        mass_filter=np.asarray(oracle.mass_filter), perturbation_value=float(value)
    )
    oracle.to_hdf5(str(output / f"KLF4_{condition}.celloracle.oracle"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oracle", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=12)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    base = co.load_hdf5(str(args.oracle))
    expression = gene_vector(base, "KLF4")
    positive = expression[expression > 0]
    oe_value = float(np.quantile(positive, 0.95)) if len(positive) else float(expression.max())
    natural_max = float(expression.max())
    oe_value = min(oe_value, natural_max)
    del base

    ko = run_one(args.oracle, "KO", 0.0, args.output, args.threads)
    oe = run_one(args.oracle, "OE", oe_value, args.output, args.threads)
    summary = {
        "status":"PASS", "engine":"CellOracle 0.18.0",
        "KO_value":0.0, "OE_value":oe_value, "natural_imputed_max":natural_max,
        "OE_rule":"95th percentile among KLF4-positive imputed values, capped at observed maximum",
        "n_propagation":3, "transition_neighbors":200, "sampled_fraction":0.30,
        "sigma_corr":0.05, "grid_n":40, "min_mass":0.01,
        "KO_cells":int(len(ko)), "OE_cells":int(len(oe)), "random_seed":20260811,
    }
    (args.output / "KLF4_perturbation_checkpoint.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
