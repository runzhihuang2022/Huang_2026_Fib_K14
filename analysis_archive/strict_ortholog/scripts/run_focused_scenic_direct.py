# -*- coding: utf-8 -*-
from __future__ import annotations

import json
import sys
from functools import partial
from pathlib import Path

import pandas as pd
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.aucell import aucell
from pyscenic.prune import df2regulons
from pyscenic.transform import module2features_auc1st_impl, modules2df
from pyscenic.utils import load_motif_annotations, modules_from_adjacencies


ROOT = Path(r"data/external/local_a82f375325f6")
RES = ROOT / "resources" / "scenic"
TAB = ROOT / "tables"


def main() -> None:
    prefix = sys.argv[1] if len(sys.argv) > 1 else "KLF4"
    if prefix == "KLF4_human_cervid_IGFL1":
        expr_name = "human_cervid_FibK14_IGFL1_expression.csv"
    else:
        expr_name = "focused_FibK14_human_harmonized_expression.csv"
    expr = pd.read_csv(RES / expr_name, index_col=0)
    adj = pd.read_csv(RES / f"{prefix}_adjacencies_top300.tsv", sep="\t")
    modules = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=True, keep_only_activating=False))
    db = FeatherRankingDatabase(
        str(RES / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
        name="hg38_10kb_gene_based",
    )
    annotations = load_motif_annotations(
        str(RES / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),
        motif_similarity_fdr=0.001,
        orthologous_identity_threshold=0.0,
    )
    transform = partial(
        module2features_auc1st_impl,
        rank_threshold=1500,
        auc_threshold=0.05,
        nes_threshold=3.0,
        filter_for_annotation=True,
    )
    motifs = modules2df(db, modules, annotations, module2features_func=transform)
    motifs.to_csv(RES / f"{prefix}_motif_enrichment.csv")
    regulons = df2regulons(motifs)
    if not regulons:
        raise RuntimeError("No motif-pruned KLF4 regulon passed NES >= 3.0")
    auc = aucell(expr, regulons, num_workers=1)
    auc.to_csv(RES / f"{prefix}_regulon_AUC.csv")

    target_rows = []
    for reg in regulons:
        for target, weight in reg.gene2weight.items():
            target_rows.append({
                "regulon": reg.name, "transcription_factor": reg.transcription_factor,
                "target_gene": target, "weight": weight,
                "is_IGFL1": target.upper() == "IGFL1",
            })
    targets = pd.DataFrame(target_rows)
    targets.to_csv(TAB / f"{prefix}_SCENIC_motif_pruned_targets.csv", index=False)
    summary = {
        "input_cells": int(expr.shape[0]), "input_genes": int(expr.shape[1]),
        "GRNBoost2_edges_total": 2193, "focused_edges_for_pruning": int(len(adj)),
        "coexpression_modules": len(modules), "motif_pruned_regulons": len(regulons),
        "motif_pruned_target_rows": int(len(targets)),
        "IGFL1_in_common_matrix": "IGFL1" in expr.columns,
        "IGFL1_in_KLF4_regulon": bool(targets["is_IGFL1"].any()) if len(targets) else False,
        "method_note": "human-harmonized Fib_K14-only focused SCENIC; hg38 v10 gene-based cisTarget; NES>=3",
    }
    (RES / f"{prefix}_focused_SCENIC_run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
