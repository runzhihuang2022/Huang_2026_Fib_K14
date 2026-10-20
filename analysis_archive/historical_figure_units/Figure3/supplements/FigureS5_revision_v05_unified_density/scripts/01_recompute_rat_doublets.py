"""Recompute rat scRNA-seq doublet scores independently within each sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scrublet as scr
from scipy import sparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--expected-rate", type=float, default=0.06)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy()
    sample_series = obs["sample_id"].astype(str)
    rows: list[pd.DataFrame] = []
    summaries: list[dict] = []

    for sample in sorted(sample_series.unique()):
        idx = np.flatnonzero(sample_series.eq(sample).to_numpy())
        counts = adata[idx].to_memory().X
        counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
        integer_fraction = (
            float(np.mean(np.isclose(counts.data, np.rint(counts.data))))
            if counts.nnz else 1.0
        )
        if integer_fraction < 0.999:
            raise ValueError(f"{sample}: X is not an integer UMI matrix")
        model = scr.Scrublet(
            counts,
            expected_doublet_rate=args.expected_rate,
            random_state=args.seed,
        )
        scores, predicted = model.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
            verbose=False,
        )
        threshold = float(model.threshold_)
        rows.append(pd.DataFrame({
            "cell_id": obs.index[idx].astype(str),
            "sample_id": sample,
            "sub_labels": obs.iloc[idx]["sub_labels"].astype(str).to_numpy(),
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
        print(
            f"{sample}: n={len(idx):,}; threshold={threshold:.4f}; "
            f"predicted={predicted.mean():.2%}",
            flush=True,
        )

    pd.concat(rows, ignore_index=True).to_csv(
        outdir / "S5B_scrublet_cell_scores.csv.gz", index=False
    )
    pd.DataFrame(summaries).to_csv(
        outdir / "S5B_scrublet_sample_summary.csv", index=False
    )


if __name__ == "__main__":
    main()
