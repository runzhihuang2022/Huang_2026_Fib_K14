#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd


BASELINE_HASHES = {
    "render_top5_tf_network_go_reference_v26.py": "7f7410d1815c1de832d482208cc0bbfdc8ef5e4a7169b34b4007b1075b0b3133",
    "render_cellchat_default_chord_v26.R": "65488c8b88e9963ba2a183b2cdb02038c647daefe2610bd504fc9af4cc4d307a",
    "compose_cellchat_default_chord_v26.py": "54943bca8b75760b39f3966fe803949b87689ce5717c02eaa2dfe991d6c06f40",
    "render_targeted_cellchat_v25.R": "cf4fd4436fc3af65d7d4ea0393172d91d2c82918762d0f43dd0f29121adf9f31",
    "render_cellchat_lr_enrichment_v25.py": "439726d3b11cad11ab470aa56ca3f5f4e1ccf9f734098e87c406262c2a52df5d",
    "build_complete_html_v26.py": "72cc645bd6c4743581d826830817c870cab3933ca599e178d1afff3af6db338a",
    "final_qa_v26.py": "28e91a77f1c4af640bb06abf73219f8776273943f2d0e38a2362cde48c6fd831",
}
GROUPS = ["Normal", "SPTDI", "DPTDI", "Mode_1", "Mode_2", "Mode_3", "Mode_4"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    pristine = root / "00_manifest/locked_renderer"
    checks = []
    for name, expected in BASELINE_HASHES.items():
        path = pristine / name
        observed = sha256(path) if path.exists() else "MISSING"
        checks.append({"file": name, "expected_sha256": expected, "observed_sha256": observed, "status": "PASS" if observed == expected else "FAIL"})
    checks_frame = pd.DataFrame(checks)
    checks_frame.to_csv(root / "00_manifest/LOCKED_RENDERER_BASELINE_SHA256.tsv", sep="\t", index=False)
    if not checks_frame.status.eq("PASS").all():
        raise RuntimeError("locked renderer baseline hash mismatch")

    marker = pd.read_csv(root / "00_manifest/marker_palette_allowlist.tsv", sep="\t")
    palette = marker[["sub_labels", "color"]].drop_duplicates("sub_labels").rename(columns={"sub_labels": "subtype"})
    palette.to_csv(root / "00_manifest/subtype_palette.tsv", sep="\t", index=False)

    rds_stage = root / "03_analysis/cellchat_locked_v27_inputs"
    rds_stage.mkdir(parents=True, exist_ok=True)
    rds_records = []
    for group in GROUPS:
        source = root / f"03_analysis/cellchat/{group}/CellChat_{group}.rds"
        target = rds_stage / source.name
        if target.exists() or target.is_symlink():
            target.unlink()
        os.symlink(source, target)
        rds_records.append({"group": group, "source": str(source), "staged": str(target), "exists": source.exists()})
    pd.DataFrame(rds_records).to_csv(root / "00_manifest/CELLCHAT_RDS_INPUT_REGISTRY.tsv", sep="\t", index=False)
    if not all(row["exists"] for row in rds_records):
        raise RuntimeError("missing CellChat RDS input")

    allowed = pd.DataFrame([
        {"working_file": "08_code/render_cellchat_default_chord_v26.R", "baseline_file": "00_manifest/locked_renderer/render_cellchat_default_chord_v26.R", "allowed_change": "real group registry and RDS filename only"},
        {"working_file": "08_code/render_targeted_cellchat_v25.R", "baseline_file": "00_manifest/locked_renderer/render_targeted_cellchat_v25.R", "allowed_change": "real group registry, RDS filename, DPTDI reference context, audit text only"},
        {"working_file": "08_code/run_locked_v27_networks.py", "baseline_file": "00_manifest/locked_renderer/render_top5_tf_network_go_reference_v26.py", "allowed_change": "external schema adapter and subtype-isolated invocation; locked renderer itself unchanged"},
        {"working_file": "08_code/compose_network_go_pairs.py", "baseline_file": "not_applicable", "allowed_change": "composition only; inserts unchanged locked PDF pages without redrawing"},
    ])
    allowed.to_csv(root / "00_manifest/LOCKED_RENDERER_ALLOWED_DIFF.tsv", sep="\t", index=False)
    contract = {
        "status": "PASS", "project": "human_epithelialization",
        "version": "v03_locked_v27_20260903", "groups": GROUPS,
        "subtype_palette_rows": int(len(palette)), "source_v02_preserved": True,
        "plotting_boundary": "v27 locked renderers; external schema adapters only",
    }
    (root / "00_manifest/LOCKED_RENDERER_CONTRACT.json").write_text(json.dumps(contract, indent=2), encoding="utf-8")
    print(json.dumps(contract, indent=2))


if __name__ == "__main__":
    main()
