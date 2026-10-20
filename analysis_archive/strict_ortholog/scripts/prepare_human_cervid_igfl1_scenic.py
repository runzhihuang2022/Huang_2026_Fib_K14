# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(r"data/external/local_a82f375325f6")
RES, TAB = ROOT / "resources" / "scenic", ROOT / "tables"
WORK = Path(r"data/external/local_7f2e175187a6")
sys.path.insert(0, str(WORK))
from build_fibk14_density_reference_style import subset_condition  # noqa: E402
from build_fibk14_umap_plus_conservation_fullwidth import specs11  # noqa: E402
from build_panelA_annotation_human_rat_flexible_others_planarian_alias import fibk14_annotation_call  # noqa: E402

GROUPS = ["Human normal", "Human wound", "Reindeer velvet", "Reindeer back"]


def mapping(a):
    out = {}
    for idx in a.var_names:
        values = [a.var.loc[idx, c] for c in ["human_gene_symbol", "gene_symbol", "original_gene_name"] if c in a.var]
        symbol = next((str(x) for x in values if str(x).lower() not in {"nan", "none", ""}), str(idx))
        symbol = symbol.split(";")[0].strip().upper()
        out.setdefault(symbol, str(idx))
    return out


def fib(a, label):
    if label.startswith("Human"):
        return fibk14_annotation_call(a)[0]
    return a.obs["Fib_K14"].astype(str).to_numpy() == "Fib_K14"


def sample(a):
    for c in ["sample_id", "samples", "sample", "h5ad_sample", "source_accession"]:
        if c in a.obs:
            return a.obs[c].astype(str).to_numpy()
    return np.asarray(["sample_unavailable"] * a.n_obs)


def main():
    specs = [s for s in specs11() if s.label in GROUPS]
    opened = [(s, ad.read_h5ad(s.path, backed="r")) for s in specs]
    maps = [mapping(a) for _, a in opened]
    common = set.intersection(*(set(m) for m in maps))
    mandatory = [g for g in ["KLF4", "IGFL1", "KRT14", "KRT5", "KRT17", "EPCAM", "TACSTD2", "TP63", "GRHL3"] if g in common]
    reindeer = opened[2][1]
    rev = {v: k for k, v in maps[2].items()}
    ranks = {rev[str(i)]: float(v) for i, v in reindeer.var.get("highly_variable_rank", pd.Series(dtype=float)).items()
             if str(i) in rev and np.isfinite(v)}
    genes = list(dict.fromkeys(mandatory + sorted(common, key=lambda g: (ranks.get(g, 1e9), g))))[:2200]
    if "IGFL1" not in genes:
        raise RuntimeError("IGFL1 is not jointly available in human and cervid Fib_K14 matrices")
    rng = np.random.default_rng(20260713)
    frames, metas = [], []
    for (spec, a), m in zip(opened, maps):
        idx = np.flatnonzero(subset_condition(a, spec) & fib(a, spec.label))
        idx = np.sort(rng.choice(idx, min(600, len(idx)), replace=False))
        x = a[idx, [m[g] for g in genes]].X
        x = x.toarray() if sp.issparse(x) else np.asarray(x)
        ids = [f"{spec.label.replace(' ', '_')}::{a.obs_names[i]}" for i in idx]
        frames.append(pd.DataFrame(np.maximum(x, 0).astype(np.float32), index=ids, columns=genes))
        sv = sample(a)
        metas.extend({"cell_id": cid, "group": spec.label, "sample_id": sv[i]} for cid, i in zip(ids, idx))
    expr = pd.concat(frames)
    expr.to_csv(RES / "human_cervid_FibK14_IGFL1_expression.csv")
    pd.DataFrame(metas).to_csv(RES / "human_cervid_FibK14_IGFL1_metadata.csv", index=False)
    pd.DataFrame({"gene": genes}).to_csv(TAB / "human_cervid_SCENIC_gene_universe.csv", index=False)
    print(expr.shape, "IGFL1", float((expr["IGFL1"] > 0).mean()), "KLF4", float((expr["KLF4"] > 0).mean()))


if __name__ == "__main__":
    main()
