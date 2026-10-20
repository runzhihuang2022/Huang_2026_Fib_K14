#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import celloracle as co
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", required=True, type=Path)
    ap.add_argument("--model-h5ad", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    links = co.load_hdf5(str(args.links))
    raw = links.links_dict["Fib_K14"].copy()
    raw = raw[raw["source"].astype(str).str.upper().eq("KLF4")].copy()
    raw["FDR"] = multipletests(raw["p"].to_numpy(), method="fdr_bh")[1]
    raw["direction"] = np.where(raw["coef_mean"] >= 0, "positive", "negative")
    raw["evidence_score"] = raw["coef_abs"] * (-np.log10(raw["p"].clip(lower=1e-300)))

    strict = links.filtered_links["Fib_K14"].copy()
    strict_targets = set(strict.loc[
        strict["source"].astype(str).str.upper().eq("KLF4"), "target"
    ].astype(str))
    raw["strict_filtered_link"] = raw["target"].astype(str).isin(strict_targets)
    raw = raw.sort_values(["evidence_score", "coef_abs"], ascending=False)
    raw.to_csv(args.output / "Fib_K14_KLF4_all_candidate_edges.tsv", sep="\t", index=False)

    significant = raw[raw["FDR"] < 0.05].copy()
    top50 = significant.head(50).copy()
    top50.insert(0, "plot_rank", np.arange(1, len(top50) + 1))
    top50.to_csv(args.output / "Fib_K14_KLF4_top50_signed_network_edges.tsv", sep="\t", index=False)
    strict_rows = raw[raw["strict_filtered_link"]].copy()
    strict_rows.to_csv(args.output / "Fib_K14_KLF4_strict_edges.tsv", sep="\t", index=False)
    pd.Series(sorted(strict_targets), name="gene").to_csv(
        args.output / "Fib_K14_KLF4_strict_targets.txt", index=False, header=False
    )

    positive_targets = top50.loc[top50["direction"].eq("positive"), "target"].astype(str)
    negative_targets = top50.loc[top50["direction"].eq("negative"), "target"].astype(str)
    pd.Series(positive_targets, name="gene").to_csv(
        args.output / "Fib_K14_KLF4_top50_positive_targets.txt", index=False, header=False
    )
    pd.Series(negative_targets, name="gene").to_csv(
        args.output / "Fib_K14_KLF4_top50_negative_targets.txt", index=False, header=False
    )

    model = ad.read_h5ad(args.model_h5ad, backed="r")
    pd.Series(model.var_names.astype(str), name="gene").to_csv(
        args.output / "CellOracle_model_gene_universe.txt", index=False, header=False
    )

    audit = {
        "status": "PASS",
        "species": "human",
        "cluster": "Fib_K14",
        "TF": "KLF4",
        "candidate_edges": int(len(raw)),
        "BH_FDR_lt_0.05": int(len(significant)),
        "significant_positive": int((significant.direction == "positive").sum()),
        "significant_negative": int((significant.direction == "negative").sum()),
        "strict_filtered_edges": int(len(strict_rows)),
        "plot_edges": int(len(top50)),
        "plot_positive_targets": int(len(positive_targets)),
        "plot_negative_targets": int(len(negative_targets)),
        "universe_genes": int(model.n_vars),
        "ranking": "coef_abs * -log10(p), among KLF4 edges with within-TF BH-FDR < 0.05",
        "strict_definition": "CellOracle Links.filter_links(p=0.001, weight=coef_abs, threshold_number=5000)",
        "enrichment_gene_set": "positive targets among the plotted top-50 signed network",
    }
    (args.output / "Fib_K14_KLF4_GRN_checkpoint.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
