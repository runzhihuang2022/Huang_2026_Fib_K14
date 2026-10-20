from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


ROOT = Path(r"data/external/local_4a7bea8abb97")
MANIFEST = ROOT / "00_manifest" / "csc_portable_input.yaml"
ATAC_MAP = ROOT / "00_manifest" / "atac_library_map.csv"
ATAC_META = Path(r"data/external/local_07bb298389cf")
OUT = ROOT / "01_preflight"


def write_json(name: str, value: dict) -> None:
    (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def run_ssh(script: str, timeout: int = 180) -> str:
    cp = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "csc-analysis-server", "bash", "-s"],
        input=script,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    return cp.stdout


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    failures: list[str] = []
    warnings: list[str] = []
    if manifest["execution"]["mode"] != "derived_integration":
        failures.append("execution.mode is not derived_integration")

    rna_record, atac_record = manifest["inputs"]
    rna_path = Path(rna_record["path"])
    for path in (rna_path, Path(atac_record["path"]), ATAC_MAP, ATAC_META):
        if not path.is_file() or path.stat().st_size == 0:
            failures.append(f"missing or zero-byte local input: {path}")

    a = ad.read_h5ad(rna_path, backed="r")
    obs = a.obs.copy()
    required_obs = ["sample_id", "Grade", "bulk_labels.x", "sub_labels"]
    for col in required_obs:
        if col not in obs.columns:
            failures.append(f"RNA obs missing {col}")
    if not a.obs_names.is_unique:
        failures.append("RNA cell IDs are not unique")
    if not a.var_names.is_unique:
        failures.append("RNA gene IDs are not unique")
    names = pd.Index(a.var_names.astype(str))
    requested = ["Klf4", "Tp63", "Krt5", "Krt14", "Epcam", "Tacstd2", "Vim", "Pdgfra", "Col1a1"]
    presence = {g: bool((names.str.lower() == g.lower()).any()) for g in requested}
    absent = [g for g, ok in presence.items() if not ok]
    if absent:
        warnings.append(f"Requested marker genes absent and skipped: {absent}")

    present = [g for g in requested if presence[g]]
    idx = [int(np.flatnonzero(names.str.lower() == g.lower())[0]) for g in present]
    small = a[:, idx].to_memory().X
    if hasattr(small, "toarray"):
        small = small.toarray()
    small = np.asarray(small) > 0
    marker = {g: small[:, j] for j, g in enumerate(present)}
    fib = (obs["bulk_labels.x"].astype(str) == "Fibroblast").to_numpy()
    epi = np.column_stack([marker[g] for g in ["Tp63", "Krt5", "Epcam", "Tacstd2"] if g in marker]).any(axis=1)
    mes = np.column_stack([marker[g] for g in ["Vim", "Pdgfra", "Col1a1"] if g in marker]).any(axis=1)
    strict = fib & marker["Klf4"] & epi & mes
    strict_table = obs.loc[strict, ["sample_id", "Grade", "sub_labels"]].copy()
    strict_by_sample = strict_table.groupby(["sample_id", "Grade"], observed=True).size().rename("strict_Fib_K14_cells").reset_index()
    strict_by_source = strict_table["sub_labels"].astype(str).value_counts().rename_axis("source_sub_labels").reset_index(name="strict_Fib_K14_cells")
    strict_by_sample.to_csv(OUT / "preflight_strict_FibK14_by_sample.csv", index=False)
    strict_by_source.to_csv(OUT / "preflight_strict_FibK14_by_source_subtype.csv", index=False)

    rna_batch = obs.groupby(["sample_id", "Grade"], observed=True).size().rename("cells").reset_index()
    rna_batch.to_csv(OUT / "rna_batch_by_group.csv", index=False)
    rna_conf = {
        "batch_key": "sample_id",
        "biological_unit": "sample_id",
        "classification": "perfect_confounding_with_biological_unit_partial_overlap_by_group",
        "decision": "Harmony sensitivity only; retain uncorrected PCA/UMAP.",
    }
    write_json("rna_harmony_confounding_audit.json", rna_conf)

    amap = pd.read_csv(ATAC_MAP)
    ameta = pd.read_csv(ATAC_META, index_col=0)
    libs = sorted(ameta["library_id"].astype(str).unique())
    mapped = sorted(amap["library_id"].astype(str).unique())
    if libs != mapped:
        failures.append(f"ATAC library map mismatch: metadata-only={sorted(set(libs)-set(mapped))}; map-only={sorted(set(mapped)-set(libs))}")
    atac_batch = ameta.groupby("library_id", observed=True).size().rename("postqc_nuclei").reset_index().merge(amap, on="library_id", how="left")
    atac_batch.to_csv(OUT / "atac_library_by_group_donor.csv", index=False)
    atac_conf = {
        "batch_key": "library_id",
        "biological_unit": "donor_id",
        "classification": "perfect_confounding_with_library_and_partial_pseudoreplication",
        "decision": "Harmony sensitivity only; use uncorrected IterativeLSI as primary ATAC geometry unless sensitivity preserves biology.",
        "paper_note": "Two libraries per rat/time point except control and day 15 with two rats/four libraries.",
    }
    write_json("atac_harmony_confounding_audit.json", atac_conf)

    remote_paths = [rna_record["remote_path"], atac_record["remote_path"]]
    remote_paths += [
        manifest["resources"]["scenic_tf_list"],
        manifest["resources"]["scenic_ranking_database"],
        manifest["resources"]["scenic_motif_annotations"],
    ]
    q = "\n".join(
        ["set -e", "test -x /root/miniconda3/envs/csc-portable-v1/bin/python", "test -x /root/miniconda3/envs/csc-pyscenic-0112/bin/pyscenic", "test -w external_asset_c91875faf0b2"]
        + [f"test -s '{p}'" for p in remote_paths]
        + ["echo PASS"]
    )
    try:
        remote_status = run_ssh(q, timeout=30).strip()
    except Exception as exc:
        remote_status = f"FAIL: {exc}"
        failures.append("remote runtime/input/resource gate failed")

    capability = {
        "scrna_qc": "applicable",
        "scatac_qc": "applicable",
        "harmony": "applicable_as_sensitivity_due_to_confounding",
        "multimodal_anchor_mapping": "applicable",
        "three_track_annotation": "applicable",
        "fibroblast_subclustering": "applicable",
        "formal_pyscenic": "applicable_after_rat_to_mouse_ortholog_generation",
        "spatial": "not_applicable",
        "cellchat": "not_requested",
        "monocle3": "not_requested",
    }
    write_json("capability_matrix.json", capability)
    audit = {
        "status": "FAIL" if failures else "PASS",
        "route": "server",
        "remote_gate": remote_status,
        "rna_shape": [int(a.n_obs), int(a.n_vars)],
        "rna_samples": int(obs["sample_id"].nunique()),
        "rna_groups": obs["Grade"].astype(str).value_counts().to_dict(),
        "atac_postqc_nuclei": int(len(ameta)),
        "atac_libraries": int(len(libs)),
        "requested_marker_presence": presence,
        "strict_Fib_K14_preflight_cells": int(strict.sum()),
        "failures": failures,
        "warnings": warnings,
        "governance_decisions": [rna_conf["decision"], atac_conf["decision"]],
    }
    write_json("derived_integration_preflight.json", audit)
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
