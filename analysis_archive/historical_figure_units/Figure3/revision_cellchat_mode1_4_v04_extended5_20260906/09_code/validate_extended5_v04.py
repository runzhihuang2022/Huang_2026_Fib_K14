#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED = {
    "Fib_SAC": ["Fib_Papi", "Fib_K14", "SAC_SG_Progenitor", "SAC_SG_Ductal", "SAC_SG_Dark"],
    "Fib_KC": ["Fib_Papi", "Fib_K14", "KC_Basal_Mig", "KC_Basal", "KC_Spinous_Mig"],
    "Fib_Endo": ["Fib_Prolif", "Fib_Myo", "Endo_Capillary", "Endo_Venous", "Endo_Arterial"],
    "Fib_Immune": ["Fib_Inflama", "Fib_CD45", "M1_Macrophage", "M2_Macrophage", "Mast"],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    table_dir = root / "02_analysis" / "required_chord" / "tables_extended5"
    audit = pd.read_csv(table_dir / "REQUIRED_CHORD_INSTANCE_AUDIT.tsv", sep="\t")
    displayed = pd.read_csv(table_dir / "ALL_DISPLAYED_REQUIRED_CHORD_LR.tsv", sep="\t")
    rows = []
    for _, item in audit.iterrows():
        block = displayed[(displayed.group == item.group) & (displayed.interface == item.interface)]
        observed = sorted(set(block.source.astype(str)) | set(block.target.astype(str))) if len(block) else []
        required = REQUIRED[item.interface]
        missing = sorted(set(required) - set(observed))
        rows.append({
            "group": item.group, "interface": item.interface, "status": item.status,
            "display_lr_n": len(block), "required_subtype_n": len(required),
            "observed_required_subtype_n": len(set(required) & set(observed)),
            "missing_required_subtypes": ";".join(missing),
            "passes_8lr_5subtype_gate": bool(item.status == "PASS" and len(block) == 8 and not missing),
        })
    detail = pd.DataFrame(rows)
    qa_dir = root / "10_qa"
    qa_dir.mkdir(parents=True, exist_ok=True)
    detail.to_csv(qa_dir / "EXTENDED5_CHORD_COMPLETENESS.tsv", sep="\t", index=False)
    ai = pd.read_csv(root / "07_figure_s3" / "Editable_AI_v08_extended5_positiveGO_final" / "QC" / "ai_reaudit.csv")
    formal = detail[detail.status == "PASS"]
    report = {
        "status": "PASS" if (len(formal) == 14 and formal.passes_8lr_5subtype_gate.all() and
                               len(detail[detail.status == "NOT_ESTIMABLE"]) == 2 and
                               (ai.Opened.astype(str).str.lower() == "true").all() and
                               (pd.to_numeric(ai.PlacedItems) == 0).all() and
                               (pd.to_numeric(ai.RasterItems) == 0).all()) else "FAIL",
        "formal_chord_panels": int(len(formal)),
        "not_estimable_positions": int(len(detail[detail.status == "NOT_ESTIMABLE"])),
        "formal_panels_with_exactly_8_lr_and_all_5_required_subtypes": int(formal.passes_8lr_5subtype_gate.sum()),
        "displayed_lr_total": int(len(displayed)),
        "ai_files": int(len(ai)),
        "ai_all_opened": bool((ai.Opened.astype(str).str.lower() == "true").all()),
        "ai_total_placed_items": int(pd.to_numeric(ai.PlacedItems).sum()),
        "ai_total_raster_items": int(pd.to_numeric(ai.RasterItems).sum()),
        "statistical_boundary": "CellChat pooled-mode predictions; GO BP one-sided hypergeometric tests with within-interface BH correction; no donor-level between-mode test and no causal claim.",
    }
    (qa_dir / "FINAL_QA_EXTENDED5_V04.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
