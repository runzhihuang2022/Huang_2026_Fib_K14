#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import gseapy as gp
import pandas as pd


def read_lines(path: Path) -> list[str]:
    return [x.strip() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, type=Path)
    ap.add_argument("--strict-genes", required=True, type=Path)
    ap.add_argument("--background", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    genes = read_lines(args.genes)
    strict = read_lines(args.strict_genes)
    background = read_lines(args.background)
    libraries = ["GO_Biological_Process_2025", "KEGG_2021_Human"]
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=libraries,
        organism="human",
        background=background,
        outdir=None,
        cutoff=1.0,
        verbose=True,
    )
    result = enr.results.copy()
    result["query_genes_n"] = len(genes)
    result["background_genes_n"] = len(background)
    result.to_csv(args.output / "Fib_K14_KLF4_positive_target_ORA_all.tsv", sep="\t", index=False)
    result[result["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").to_csv(
        args.output / "Fib_K14_KLF4_positive_target_ORA_FDR05.tsv", sep="\t", index=False
    )

    # Four strict targets are retained as a separate evidence tier.  The list
    # is too small for stable functional interpretation; run and retain its
    # table only as an audit, never merge it with the expanded analysis.
    strict_enr = gp.enrichr(
        gene_list=strict,
        gene_sets=libraries,
        organism="human",
        background=background,
        outdir=None,
        cutoff=1.0,
        verbose=True,
    )
    strict_result = strict_enr.results.copy()
    strict_result["query_genes_n"] = len(strict)
    strict_result["background_genes_n"] = len(background)
    strict_result.to_csv(args.output / "Fib_K14_KLF4_strict4_ORA_audit.tsv", sep="\t", index=False)

    checkpoint = {
        "status": "PASS",
        "engine": f"GSEApy {gp.__version__} / Enrichr ORA",
        "species": "human",
        "libraries": libraries,
        "query_positive_network_targets": len(genes),
        "strict_targets": len(strict),
        "background_genes": len(background),
        "expanded_FDR05_terms": int((result["Adjusted P-value"] < 0.05).sum()),
        "strict_FDR05_terms": int((strict_result["Adjusted P-value"] < 0.05).sum()),
        "statistical_label": "gene-set enrichment statistic; not donor-replicated clinical inference",
    }
    (args.output / "Fib_K14_KLF4_ORA_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
