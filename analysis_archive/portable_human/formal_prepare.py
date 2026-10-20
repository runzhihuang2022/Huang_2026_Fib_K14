#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


GROUP_ORDER = ["Normal", "SPTDI", "DPTDI"]
TRACKS = ["bulk_labels", "sub_labels"]
LINEAGE_SUBTYPES = {
    "epithelial": [
        "KC_Basal", "KC_Basal_Mig", "KC_Spinous_Mig",
    ],
    "fibroblast": [
        "Fib_CD45", "Fib_EN1", "Fib_Fasci", "Fib_Inflama", "Fib_K14",
        "Fib_Myo", "Fib_Papi", "Fib_Prolif", "Fib_SFRP2",
    ],
    "immune": [
        "M1_Macrophage", "M2_Macrophage", "LAM",
    ],
    "endothelial": [
        "Endo_Arterial", "Endo_Capillary", "Endo_Lymphatic",
        "Endo_Prolif", "Endo_Venous",
    ],
}


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def categorical_summary(obs: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    for key in keys:
        counts = obs[key].astype(str).value_counts(dropna=False)
        for value, count in counts.items():
            rows.append({"field": key, "value": value, "cells": int(count)})
    return pd.DataFrame(rows)


def metric_summary(obs: pd.DataFrame, group_keys: list[str]) -> pd.DataFrame:
    candidates = [
        "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
        "total_counts", "n_genes_by_counts", "pct_counts_mt",
    ]
    metrics = [x for x in candidates if x in obs.columns]
    rows = []
    if not metrics:
        return pd.DataFrame(
            columns=group_keys + ["metric", "n", "median", "q1", "q3", "min", "max"]
        )
    for identity, frame in obs.groupby(group_keys, observed=True, sort=False):
        identity = identity if isinstance(identity, tuple) else (identity,)
        base = dict(zip(group_keys, map(str, identity)))
        for metric in metrics:
            values = pd.to_numeric(frame[metric], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append({
                **base, "metric": metric, "n": int(values.size),
                "median": float(values.median()),
                "q1": float(values.quantile(0.25)),
                "q3": float(values.quantile(0.75)),
                "min": float(values.min()), "max": float(values.max()),
            })
    return pd.DataFrame(rows)


def composition_tables(obs: pd.DataFrame, track: str, out: Path) -> None:
    for unit in ["patient", "sample_id"]:
        counts = (
            obs.groupby([unit, "group", track], observed=True)
            .size().rename("cells").reset_index()
        )
        totals = (
            obs.groupby([unit, "group"], observed=True)
            .size().rename("total_cells").reset_index()
        )
        table = counts.merge(totals, on=[unit, "group"], how="left")
        table["proportion"] = table["cells"] / table["total_cells"]
        write_tsv(table, out / f"{track}_{unit}_composition.tsv")


def pseudobulk_h5(adata: ad.AnnData, out: Path) -> dict:
    obs = adata.obs
    x = adata.X
    if not sparse.issparse(x):
        raise RuntimeError("Canonical X is not sparse; refusing unexpected dense conversion.")
    x = x.tocsr()
    genes = np.asarray(adata.var_names.astype(str))
    h5_path = out / "pseudobulk_counts.h5"
    inventory = {}
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("genes", data=genes.astype(object), dtype=string_dtype)
        for track in TRACKS:
            meta = obs[["patient", "group", track]].astype(str).copy()
            meta["key"] = (
                meta["patient"] + "||" + meta["group"] + "||" + meta[track]
            )
            keys = pd.unique(meta["key"])
            grp = handle.create_group(track)
            counts_ds = grp.create_dataset(
                "counts", shape=(len(keys), adata.n_vars), dtype="i8",
                chunks=(1, min(adata.n_vars, 8192)),
                compression="gzip", compression_opts=4, shuffle=True,
            )
            records = []
            key_values = meta["key"].to_numpy()
            for row_index, key in enumerate(keys):
                cell_index = np.flatnonzero(key_values == key)
                aggregate = np.asarray(x[cell_index].sum(axis=0)).ravel()
                if not np.allclose(aggregate, np.rint(aggregate)):
                    raise RuntimeError(f"Non-integer pseudobulk aggregate for {key}")
                counts_ds[row_index, :] = np.rint(aggregate).astype(np.int64)
                patient, group, label = key.split("||", 2)
                records.append({
                    "row_index": row_index, "patient": patient, "group": group,
                    "label": label, "cells": int(cell_index.size),
                })
            metadata = pd.DataFrame(records)
            write_tsv(metadata, out / f"{track}_pseudobulk_metadata.tsv")
            inventory[track] = {
                "pseudobulk_rows": int(len(metadata)),
                "labels": int(metadata["label"].nunique()),
                "patients": int(metadata["patient"].nunique()),
            }
    inventory["h5_path"] = str(h5_path)
    inventory["genes"] = int(adata.n_vars)
    return inventory


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.canonical)
    required = [
        "sample_id", "patient", "group", "orig.ident", "time_point",
        "tissue_group", "bulk_labels", "sub_labels",
    ]
    missing = [x for x in required if x not in adata.obs]
    if missing:
        raise RuntimeError(f"Missing canonical obs fields: {missing}")
    if adata.n_obs != 279305 or adata.n_vars != 52399:
        raise RuntimeError(f"Frozen shape violated: {adata.shape}")
    if "X_umap" not in adata.obsm:
        raise RuntimeError("Frozen X_umap missing.")
    observed_groups = set(adata.obs["group"].astype(str))
    if observed_groups != set(GROUP_ORDER):
        raise RuntimeError(f"Unexpected clinical groups: {sorted(observed_groups)}")

    obs = adata.obs.copy()
    obs["group"] = pd.Categorical(obs["group"].astype(str), GROUP_ORDER, ordered=True)
    write_tsv(
        categorical_summary(
            obs,
            ["group", "sample_id", "patient", "orig.ident", "time_point",
             "tissue_group", "bulk_labels", "sub_labels"],
        ),
        out / "M01_categorical_counts.tsv",
    )
    write_tsv(
        metric_summary(obs, ["group", "sample_id", "patient"]),
        out / "M01_stored_qc_metric_summary.tsv",
    )
    sample_summary = (
        obs.groupby(["sample_id", "patient", "group"], observed=True)
        .size().rename("cells").reset_index()
    )
    write_tsv(sample_summary, out / "M01_sample_summary.tsv")
    patient_summary = (
        obs.groupby(["patient", "group"], observed=True)
        .size().rename("cells").reset_index()
    )
    write_tsv(patient_summary, out / "M01_patient_group_summary.tsv")

    for track in TRACKS:
        counts = obs[track].astype(str).value_counts().rename_axis(track).reset_index(name="cells")
        write_tsv(counts, out / f"M02_{track}_counts.tsv")
        composition_tables(obs, track, out)
    cross = pd.crosstab(
        obs["bulk_labels"].astype(str), obs["sub_labels"].astype(str),
        dropna=False,
    ).rename_axis("bulk_labels").reset_index()
    write_tsv(cross, out / "M02_bulk_sub_cross_table.tsv")

    lineage_rows = []
    observed_sub = set(obs["sub_labels"].astype(str))
    for lineage, labels in LINEAGE_SUBTYPES.items():
        for label in labels:
            lineage_rows.append({
                "lineage": lineage, "sub_label": label,
                "exists": label in observed_sub,
                "cells": int((obs["sub_labels"].astype(str) == label).sum()),
            })
    write_tsv(pd.DataFrame(lineage_rows), out / "M07_lineage_subtype_inventory.tsv")

    inventory = pseudobulk_h5(adata, out)
    status = {
        "status": "PASS",
        "canonical": str(args.canonical),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "preserved_reductions": sorted(adata.obsm.keys()),
        "tracks": TRACKS,
        "pseudobulk": inventory,
        "prohibited_operations_run": [],
    }
    (out / "formal_prepare_checkpoint.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
