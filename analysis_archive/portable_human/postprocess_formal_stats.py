#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-dir", required=True, type=Path)
    args = parser.parse_args()
    root = args.formal_dir.resolve()

    da_path = root / "M04_differential_abundance.tsv"
    da = pd.read_csv(da_path, sep="\t")
    da["approx_se_from_qlF"] = np.where(
        da["F"] > 0, np.abs(da["logFC"]) / np.sqrt(da["F"]), np.nan
    )
    da["ci95_low"] = da["logFC"] - 1.96 * da["approx_se_from_qlF"]
    da["ci95_high"] = da["logFC"] + 1.96 * da["approx_se_from_qlF"]
    da["ci_method"] = "approximate_Wald_from_edgeR_QLF"
    da.to_csv(da_path, sep="\t", index=False)

    index_path = root / "M05_result_index.tsv"
    index = pd.read_csv(index_path, sep="\t")
    rewritten = 0
    for path_text in index["result_file"]:
        path = Path(path_text)
        result = pd.read_csv(path, sep="\t")
        result["ci95_low"] = result["log2FoldChange"] - 1.96 * result["lfcSE"]
        result["ci95_high"] = result["log2FoldChange"] + 1.96 * result["lfcSE"]
        result["ci_method"] = "DESeq2_Wald_95CI"
        result.to_csv(path, sep="\t", index=False, compression="gzip")
        rewritten += 1

    model_registry = pd.DataFrame([
        {
            "module": "M04", "analysis": "all_patient_unpaired",
            "statistical_unit": "patient_group",
            "model_formula": "~ group",
            "offset": "log(total_cells_per_patient_group)",
            "effect": "edgeR log2 fold change in relative abundance",
            "test": "edgeR quasi-likelihood F test",
            "multiplicity": "Benjamini-Hochberg within track and contrast",
        },
        {
            "module": "M04", "analysis": "paired_patient_sensitivity",
            "statistical_unit": "paired_patient_group",
            "model_formula": "~ patient + group",
            "offset": "log(total_cells_per_patient_group)",
            "effect": "edgeR log2 fold change in relative abundance",
            "test": "edgeR quasi-likelihood F test",
            "multiplicity": "Benjamini-Hochberg within track and contrast",
        },
        {
            "module": "M05", "analysis": "all_patient_unpaired",
            "statistical_unit": "patient_group_celltype_pseudobulk",
            "model_formula": "~ group",
            "offset": "DESeq2 median-ratio size factor",
            "effect": "DESeq2 log2FoldChange",
            "test": "DESeq2 Wald test",
            "multiplicity": "Benjamini-Hochberg within label and contrast",
        },
        {
            "module": "M05", "analysis": "paired_patient_sensitivity",
            "statistical_unit": "paired_patient_group_celltype_pseudobulk",
            "model_formula": "~ patient + group",
            "offset": "DESeq2 median-ratio size factor",
            "effect": "DESeq2 log2FoldChange",
            "test": "DESeq2 Wald test",
            "multiplicity": "Benjamini-Hochberg within label and contrast",
        },
    ])
    model_registry.to_csv(root / "formal_model_registry.tsv", sep="\t", index=False)
    checkpoint = {
        "status": "PASS",
        "da_rows": int(len(da)),
        "de_result_files_rewritten_with_ci": int(rewritten),
        "cells_as_replicates": False,
        "main_replicate_unit": "patient_group",
        "paired_sensitivity_patients": ["patient45", "patient52", "patient57"],
        "confidence_intervals": {
            "M04": "approximate Wald CI reconstructed from one-df edgeR QL F",
            "M05": "DESeq2 Wald CI from lfcSE",
        },
    }
    (root / "formal_statistics_postprocess_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
