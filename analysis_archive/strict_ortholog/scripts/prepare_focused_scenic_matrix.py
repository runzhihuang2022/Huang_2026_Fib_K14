# -*- coding: utf-8 -*-
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(r"data/external/local_a82f375325f6")
RES = ROOT / "resources" / "scenic"
TAB = ROOT / "tables"
WORK = Path(r"data/external/local_7f2e175187a6")
sys.path.insert(0, str(WORK))

from build_fibk14_density_reference_style import marker_call, subset_condition  # noqa: E402
from build_fibk14_umap_plus_conservation_fullwidth import specs11  # noqa: E402
from build_panelA_annotation_human_rat_flexible_others_planarian_alias import fibk14_annotation_call  # noqa: E402


KEEP_GROUPS = [
    "Human normal", "Human wound", "Rat embryo", "Rat adult",
    "Reindeer velvet", "Reindeer back",
]


def display(label: str) -> str:
    return label.replace("hpa1", " h").replace("dpa1", " d")


def gene_map(a: ad.AnnData) -> dict[str, str]:
    columns = [c for c in ["human_gene_symbol", "gene_symbol", "original_gene_name"] if c in a.var]
    mapping: dict[str, str] = {}
    for idx in a.var_names:
        symbol = str(idx)
        for col in columns:
            value = str(a.var.loc[idx, col])
            if value and value.lower() not in {"nan", "none", ""}:
                symbol = value
                break
        symbol = symbol.split(";")[0].strip().upper()
        if symbol and symbol not in mapping:
            mapping[symbol] = str(idx)
    return mapping


def fib_mask(a: ad.AnnData, label: str) -> np.ndarray:
    if label.startswith(("Human", "Rat")):
        return fibk14_annotation_call(a)[0]
    if "Fib_K14" in a.obs:
        return a.obs["Fib_K14"].astype(str).to_numpy() == "Fib_K14"
    return marker_call(a, allow_krt5_fallback=False)[0]


def sample_values(a: ad.AnnData) -> np.ndarray:
    for column in ["sample_id", "samples", "sample", "h5ad_sample", "orig.ident", "orig.ident.x", "source_accession"]:
        if column in a.obs:
            return a.obs[column].astype(str).to_numpy()
    return np.asarray(["sample_unavailable"] * a.n_obs)


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    selected = [s for s in specs11() if s.label in KEEP_GROUPS]
    opened = []
    maps = []
    for spec in selected:
        a = ad.read_h5ad(spec.path, backed="r")
        opened.append((spec, a))
        maps.append(gene_map(a))

    common = set(maps[0])
    for m in maps[1:]:
        common &= set(m)
    # Keep human cisTarget-compatible symbols and prioritize the reindeer HVG order.
    reindeer = opened[4][1]
    reindeer_map = maps[4]
    rank = {}
    if "highly_variable_rank" in reindeer.var:
        reverse = {v: k for k, v in reindeer_map.items()}
        for idx, value in reindeer.var["highly_variable_rank"].items():
            symbol = reverse.get(str(idx))
            if symbol and np.isfinite(value):
                rank[symbol] = float(value)
    ordered = sorted(common, key=lambda g: (rank.get(g, 1e9), g))
    mandatory = [g for g in ["KLF4", "IGFL1", "KRT14", "KRT5", "KRT17", "EPCAM", "TACSTD2",
                                  "TP63", "GRHL3", "AREG", "COL1A1", "COL17A1", "CD34"] if g in common]
    genes = list(dict.fromkeys(mandatory + ordered))[:2200]

    rng = np.random.default_rng(20260713)
    frames = []
    meta = []
    audit = []
    for (spec, a), mapping in zip(opened, maps):
        keep = subset_condition(a, spec)
        fib = fib_mask(a, spec.label)
        idx = np.flatnonzero(keep & fib)
        if len(idx) > 600:
            idx = np.sort(rng.choice(idx, 600, replace=False))
        feature_ids = [mapping[g] for g in genes]
        x = a[idx, feature_ids].X
        x = x.toarray() if sp.issparse(x) else np.asarray(x)
        x = np.maximum(x, 0).astype(np.float32)
        names = [f"{spec.label.replace(' ', '_')}::{a.obs_names[i]}" for i in idx]
        frames.append(pd.DataFrame(x, index=names, columns=genes))
        samples = sample_values(a)
        meta.extend({"cell_id": name, "group": spec.label, "species": spec.cohort,
                     "sample_id": samples[i], "analysis_context": "Fib_K14 only"}
                    for name, i in zip(names, idx))
        audit.append({"group": spec.label, "input_path": str(spec.path), "fibk14_cells_available": int((keep & fib).sum()),
                      "cells_subsampled": len(idx), "common_genes": len(genes), "KLF4_in_matrix": "KLF4" in genes,
                      "IGFL1_in_matrix": "IGFL1" in genes})

    expr = pd.concat(frames, axis=0)
    expr.to_csv(RES / "focused_FibK14_human_harmonized_expression.csv")
    pd.DataFrame(meta).to_csv(RES / "focused_FibK14_cell_metadata.csv", index=False)
    pd.DataFrame(audit).to_csv(TAB / "focused_SCENIC_input_audit.csv", index=False)
    (RES / "KLF4_only_TFs.txt").write_text("KLF4\n", encoding="ascii")
    manifest = {
        "matrix_shape": list(expr.shape), "groups": KEEP_GROUPS, "genes": genes,
        "scope": "human-harmonized Fib_K14-only focused SCENIC",
        "TFs": ["KLF4"], "IGFL1_role": "candidate KLF4 target; not a TF",
    }
    (RES / "focused_SCENIC_input_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"matrix": str(RES / 'focused_FibK14_human_harmonized_expression.csv'),
                      "shape": list(expr.shape), "audit": audit}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
