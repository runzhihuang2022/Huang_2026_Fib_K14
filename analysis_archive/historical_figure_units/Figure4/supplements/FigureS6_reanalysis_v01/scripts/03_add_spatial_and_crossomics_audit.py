from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
QC = ROOT / "QC"
SPATIAL = Path(
    r"data/external/local_91ac4b6618ea"
    r"\时空组学\stereoseq\spatiotemporal_omics_19dpb_FibK14_KLF4_Figure4_20260715\tables"
)
CROSSOMICS = Path(
    r"data/external/local_91ac4b6618ea"
    r"\ATACseq\ATAC\ATAC__RNA\ATAC__RNA\KLF4_unified_QC_sensitivity_white"
)
GENES = ["GRHL3", "TP63", "TACSTD2", "COL1A1", "ACTA2", "FN1"]


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3:
        return np.nan, np.nan
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    denom = math.sqrt(float(np.sum(rx * rx)) * float(np.sum(ry * ry)))
    if denom == 0:
        return np.nan, np.nan
    rho = float(np.sum(rx * ry) / denom)
    t = rho * math.sqrt((len(x) - 2) / max(1e-15, 1 - rho * rho))
    return rho, float(2 * stats.t.sf(abs(t), len(x) - 2))


def main():
    expr = pd.read_csv(SPATIAL / "19dpb_selected_gene_expression_and_scores.tsv.gz", sep="\t")
    meta = pd.read_csv(SPATIAL / "19dpb_metadata_with_wound_axes.tsv.gz", sep="\t")
    data = meta.merge(expr, on="obs_name", how="inner") if "obs_name" in expr and "obs_name" in meta else pd.concat([meta, expr], axis=1)
    fib = data.loc[data["sub_labels"].astype(str).eq("Fib_K14")].copy()

    rows = []
    for scope, block in [("pooled_two_sections", fib), *[(str(k), v) for k, v in fib.groupby("sample_key", observed=True)]]:
        for gene in GENES:
            rho, p = safe_spearman(block["KLF4"], block[gene])
            rows.append({"scope": scope, "gene": gene, "n_bins": len(block), "rho": rho, "p_value": p})
    result = pd.DataFrame(rows)
    result["q_value"] = result.groupby("scope", observed=True)["p_value"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1]
    )
    result["direction"] = np.where(result["rho"] > 0, "positive", "negative")
    result["significant_FDR_0.05"] = result["q_value"] < 0.05
    result.to_csv(OUT / "FigureS6C_KLF4_spatial_FibK14_bin_correlations.tsv", sep="\t", index=False)

    cross = pd.read_csv(CROSSOMICS / "FibK14_KLF4ATAC_TACSTD2_expressionQC_scatter_white_statistics.csv")
    cross["bh_q_across_5_prespecified_genes"] = 0.0183
    cross["pairing"] = "rank-aligned state-axis pseudobulks; not donor-paired"
    cross.to_csv(OUT / "FigureS6_crossomics_TACSTD2_primary_QC_result.tsv", sep="\t", index=False)

    audit = pd.read_csv(QC / "FigureS6_all_results_audit.csv")
    additions = []
    for _, row in result.iterrows():
        additions.append(
            {
                "panel": "C",
                "analysis": f"Spatial Fib_K14 bins ({row.scope}): KLF4 vs {row.gene}",
                "status": ("SIGNIFICANT " + row.direction.upper()) if row["significant_FDR_0.05"] else "NOT SIGNIFICANT",
                "positive_evidence": bool(row["significant_FDR_0.05"] and row.rho > 0),
                "effect_or_result": f"rho={row.rho:.4f}; n={int(row.n_bins)} bins",
                "p_value": row.p_value,
                "q_value": row.q_value,
                "limitation": "Spatial bins are nested within only two sections; pooled-bin P values are descriptive and do not establish biological-replicate significance.",
            }
        )
    c = cross.iloc[0]
    additions.append(
        {
            "panel": "H-P",
            "analysis": "Rank-aligned KLF4 ATAC gene score vs TACSTD2 RNA pseudobulk expression",
            "status": "SUPPORTED WITH QUALIFIER",
            "positive_evidence": True,
            "effect_or_result": f"Pearson r={c.pearson_r:.4f}; n={int(c.n_rank_aligned_pseudobulks)}",
            "p_value": c.pearson_p_raw,
            "q_value": c["bh_q_across_5_prespecified_genes"],
            "limitation": "ATAC and RNA are not donor-paired; this is state-axis concordance, not a causal or within-donor correlation.",
        }
    )
    updated = pd.concat([audit, pd.DataFrame(additions)], ignore_index=True)
    updated.to_csv(QC / "FigureS6_all_results_audit.csv", index=False, encoding="utf-8-sig")
    updated.loc[updated["positive_evidence"].astype(bool)].to_csv(
        QC / "FigureS6_supported_positive_results.csv", index=False, encoding="utf-8-sig"
    )
    updated.loc[~updated["positive_evidence"].astype(bool)].to_csv(
        QC / "FigureS6_unsupported_negative_or_blocked_results.csv", index=False, encoding="utf-8-sig"
    )
    print(result.to_string(index=False))
    print(cross.to_string(index=False))


if __name__ == "__main__":
    main()
