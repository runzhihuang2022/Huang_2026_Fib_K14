from __future__ import annotations

import csv
import gc
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from figure_svg_utils import ROOT


ANNOTATED_H5AD = Path(os.environ.get(
    "FIG1_ANNOTATED_H5AD",
    r"data/external/local_f1a3b8b8d8e7",
))
CACHE = ROOT / "source_data/Figure1EF_registered11_scores"
SUMMARY = ROOT / "source_data/Figure1F_section_level_scores.csv"
METADATA = ROOT / "source_data/Figure1EF_score_provenance.json"

SAMPLES = {
    "5dpb": "BW32_A01597A3_SDSDB_5dpb",
    "12dpb_DPTDI1": "BW13_1_B1_DSDB_12dpb",
    "12dpb_SPTDI1": "BW14_1_C1_SSDB_12dpb",
    "12dpb_DPTDI2": "BW13_A3_DSDB_12dpb",
    "12dpb_SPTDI2": "BW14_B3_SSDB_12dpb",
    "19dpb": "BW15D_C6_SDSDB_19dpb",
    "19dpb_p1": "BW15D_1_D1_SDSDB_19dpb",
    "26dpb_p1": "BW81_C02846B6_SDSDB_26dpb_part1",
    "26dpb_p2": "BW81_C02846B6_SDSDB_26dpb_part2",
    "2mph": "2mph_A03699G6.SCT",
    "Normal": "NS_C02847B1",
}
TIMEPOINT = {
    "5dpb": "5 dpb", "12dpb_DPTDI1": "12 dpb", "12dpb_SPTDI1": "12 dpb",
    "12dpb_DPTDI2": "12 dpb", "12dpb_SPTDI2": "12 dpb", "19dpb": "19 dpb",
    "19dpb_p1": "19 dpb", "26dpb_p1": "26 dpb", "26dpb_p2": "26 dpb",
    "2mph": "2 mph", "Normal": "Normal",
}
GENE_SETS = {
    "Epithelialization": ["KRT14", "KRT5", "KRT10", "KRT1", "IVL", "SBSN", "LOR", "FLG", "AREG", "MMP3", "MMP9", "KRT17", "KRT6A"],
    "Skin development": ["TP63", "WNT10A", "LEF1", "SHH", "BMP2", "EDAR", "SOX9"],
    "Inflammatory response": ["IL1B", "IL6", "TNF", "CXCL8", "CCL2", "CXCL1", "CXCL2"],
    "Vascularization": ["VEGFA", "PECAM1", "CD34", "KDR", "FLT1", "ANGPT1", "VWF"],
    "Matrix remodeling": ["ACTA2", "COL1A1", "COL1A2", "COL3A1", "FN1", "TGFB1", "CTGF"],
}
REPRESENTATIVE = {
    "5dpb": 135,
    "19dpb": -225,
    "12dpb_SPTDI2": 135,
    "2mph": 45,
}


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    print(f"Loading annotated object: {ANNOTATED_H5AD}", flush=True)
    adata = sc.read_h5ad(ANNOTATED_H5AD)
    sample_col = "sample_batch_new"
    registered = set(SAMPLES.values())
    mask = adata.obs[sample_col].astype(str).isin(registered).to_numpy()
    adata = adata[mask].copy()
    print(f"Registered subset: {adata.n_obs:,} bins x {adata.n_vars:,} genes", flush=True)

    valid_by_set = {}
    score_cols = []
    for name, genes in GENE_SETS.items():
        valid = [gene for gene in genes if gene in adata.var_names]
        if not valid:
            raise RuntimeError(f"No genes found for {name}")
        score = "Score_" + name.replace(" ", "_")
        print(f"Scoring {name}: {len(valid)}/{len(genes)} genes", flush=True)
        sc.tl.score_genes(
            adata,
            gene_list=valid,
            score_name=score,
            ctrl_size=50,
            n_bins=25,
            random_state=0,
            use_raw=False,
        )
        valid_by_set[name] = valid
        score_cols.append(score)

    id_to_key = {v: k for k, v in SAMPLES.items()}
    rows = []
    for sample_id in SAMPLES.values():
        key = id_to_key[sample_id]
        sub = adata.obs.loc[adata.obs[sample_col].astype(str).eq(sample_id), score_cols]
        row = {"sample_key": key, "sample_id": sample_id, "timepoint": TIMEPOINT[key], "n_bins": len(sub)}
        for name, score in zip(GENE_SETS, score_cols):
            row[name] = float(sub[score].mean())
        rows.append(row)
    summary = pd.DataFrame(rows)
    for name in GENE_SETS:
        sd = summary[name].std(ddof=1)
        summary[f"{name} z"] = (summary[name] - summary[name].mean()) / sd if sd > 0 else 0.0
    summary.to_csv(SUMMARY, index=False, encoding="utf-8-sig")

    all_epi = adata.obs[score_cols[0]].to_numpy()
    vmin, vmax = np.nanpercentile(all_epi, [1, 99])
    for key, angle in REPRESENTATIVE.items():
        sid = SAMPLES[key]
        keep = adata.obs[sample_col].astype(str).eq(sid).to_numpy()
        coords = np.asarray(adata.obsm["spatial"])[keep]
        theta = np.deg2rad(angle)
        x = coords[:, 0] * np.cos(theta) - coords[:, 1] * np.sin(theta)
        y = coords[:, 0] * np.sin(theta) + coords[:, 1] * np.cos(theta)
        score = adata.obs.loc[keep, score_cols[0]].to_numpy()
        np.savez_compressed(CACHE / f"{key}.npz", x=x, y=y, score=score)

    METADATA.write_text(json.dumps({
        "source_h5ad": str(ANNOTATED_H5AD),
        "subset_rule": "exact whitelist of 11 registered human Stereo-seq sample_batch_new values",
        "n_bins": int(adata.n_obs),
        "sample_counts": summary.set_index("sample_key")["n_bins"].astype(int).to_dict(),
        "gene_sets_requested": GENE_SETS,
        "genes_found": valid_by_set,
        "score_method": "scanpy.tl.score_genes; ctrl_size=50; n_bins=25; random_state=0; use_raw=False",
        "epithelialization_display_limits_registered11_percentile_1_99": [float(vmin), float(vmax)],
        "temporal_unit": "section-level mean; z-standardized across 11 sections within each gene set for display",
        "excluded_timepoints": ["6 mph", "9 mph"],
    }, indent=2), encoding="utf-8")
    print(SUMMARY, flush=True)
    print(METADATA, flush=True)
    del adata
    gc.collect()


if __name__ == "__main__":
    main()
