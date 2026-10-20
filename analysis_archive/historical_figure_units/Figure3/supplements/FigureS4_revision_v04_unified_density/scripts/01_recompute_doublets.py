"""Recompute doublet scores per human scRNA-seq sample from integer UMI counts.

The registered h5ad stores integer counts in ``X`` but has no ``raw`` object,
counts layer, or precomputed doublet score.  Scrublet is therefore run
independently for each ``sample_id`` to avoid cross-sample synthetic doublets.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scrublet as scr
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--expected-rate", type=float, default=0.06)
    p.add_argument("--seed", type=int, default=17)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy()
    if "sample_id" not in obs or "sub_labels" not in obs:
        raise KeyError("Required obs fields sample_id/sub_labels are absent")

    rows: list[pd.DataFrame] = []
    summaries: list[dict] = []
    for sample in sorted(obs["sample_id"].astype(str).unique()):
        mask = obs["sample_id"].astype(str).eq(sample).to_numpy()
        idx = np.flatnonzero(mask)
        sample_adata = adata[idx, :].to_memory()
        counts = sample_adata.X
        if not sparse.issparse(counts):
            counts = sparse.csr_matrix(counts)
        counts = counts.tocsr()
        nz = counts.data
        integer_fraction = float(np.mean(np.isclose(nz, np.rint(nz)))) if nz.size else 1.0
        if integer_fraction < 0.999:
            raise ValueError(f"{sample}: X is not an integer-count matrix")

        scrub = scr.Scrublet(
            counts,
            expected_doublet_rate=args.expected_rate,
            random_state=args.seed,
        )
        scores, predicted = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
            verbose=False,
        )
        threshold = float(scrub.threshold_)
        cell_obs = obs.iloc[idx]
        rows.append(pd.DataFrame({
            "cell_id": cell_obs.index.astype(str),
            "sample_id": sample,
            "sub_labels": cell_obs["sub_labels"].astype(str).to_numpy(),
            "doublet_score": scores,
            "predicted_doublet": predicted,
            "threshold": threshold,
        }))
        summaries.append({
            "sample_id": sample,
            "n_cells": len(idx),
            "expected_doublet_rate": args.expected_rate,
            "threshold": threshold,
            "predicted_doublets": int(predicted.sum()),
            "predicted_doublet_fraction": float(predicted.mean()),
            "integer_fraction_X_nonzero": integer_fraction,
        })
        print(f"{sample}: n={len(idx):,}, threshold={threshold:.4f}, "
              f"predicted={predicted.mean():.3%}", flush=True)

    all_scores = pd.concat(rows, ignore_index=True)
    all_scores.to_csv(outdir / "S4B_scrublet_cell_scores.csv.gz", index=False)
    pd.DataFrame(summaries).to_csv(outdir / "S4B_scrublet_sample_summary.csv", index=False)


if __name__ == "__main__":
    main()
