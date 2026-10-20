#!/usr/bin/env python
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd


def table_html(path: Path, max_rows: int = 100) -> str:
    try:
        frame = pd.read_csv(path, sep="\t")
    except Exception as error:
        return f"<p class='error'>{html.escape(str(error))}</p>"
    note = ""
    if len(frame) > max_rows:
        note = f"<p>Showing {max_rows} of {len(frame):,} rows.</p>"
        frame = frame.head(max_rows)
    return note + frame.to_html(index=False, border=0, classes="data")


def json_html(path: Path) -> str:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return f"<pre>{html.escape(json.dumps(payload, indent=2, ensure_ascii=False))}</pre>"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.analysis_root.resolve()
    sections = [
        ("Canonical preparation", root / "00_prepare" / "formal_prepare_checkpoint.json", "json"),
        ("Differential abundance", root / "01_formal_inference" / "M04_checkpoint.tsv", "tsv"),
        ("Pseudobulk DE", root / "01_formal_inference" / "M05_checkpoint.tsv", "tsv"),
        ("Formal inference modules", root / "01_formal_inference" / "formal_inference_checkpoint.tsv", "tsv"),
        ("GSEA", root / "01_formal_inference" / "M06_GSEA" / "M06_GSEA_checkpoint.tsv", "tsv"),
        ("Advanced module inputs", root / "02_advanced_prepare" / "advanced_prepare_checkpoint.json", "json"),
        ("CellChat", root / "03_cellchat" / "M08_cellchat_checkpoint.tsv", "tsv"),
        ("Monocle3", root / "04_monocle3" / "M09_monocle3_checkpoint.tsv", "tsv"),
        ("Stereo-seq integration", root / "05_spatial_projection" / "M12_spatial_projection_checkpoint.json", "json"),
    ]
    missing = [str(path) for _, path, _ in sections if not path.exists()]
    if missing:
        raise RuntimeError("Missing formal checkpoints:\n" + "\n".join(missing))
    body = []
    for title, path, kind in sections:
        rendered = json_html(path) if kind == "json" else table_html(path)
        body.append(
            f"<section><h2>{html.escape(title)}</h2>"
            f"<p><code>{html.escape(str(path.relative_to(root)))}</code></p>"
            f"{rendered}</section>"
        )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>human_epithelialization formal statistics checkpoint</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 1120px; margin: 24px auto;
       color: #202124; line-height: 1.4; }}
h1 {{ font-size: 22px; }} h2 {{ font-size: 16px; margin-top: 28px; }}
.notice {{ border-left: 5px solid #1b9e77; background: #edf8f3; padding: 12px; }}
table.data {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
table.data th, table.data td {{ border: 1px solid #ddd; padding: 4px 6px; }}
table.data th {{ background: #f5f5f5; position: sticky; top: 0; }}
pre {{ white-space: pre-wrap; background: #f7f7f7; padding: 12px; }}
code {{ font-family: Consolas, monospace; }}
</style>
</head>
<body>
<h1>human_epithelialization — formal statistics checkpoint</h1>
<div class="notice">
Statistics use patient or patient×group pseudobulk units. Cell counts are
descriptive only. No QC filtering, PCA, Harmony, neighbors, clustering, UMAP,
automatic annotation, or approved lock creation was run.
</div>
{''.join(body)}
</body></html>"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
