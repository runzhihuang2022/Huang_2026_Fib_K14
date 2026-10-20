#!/usr/bin/env python3
"""Batch KFBF to lossless pyramidal OME-TIFF with per-file audits."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from zarr_to_ome_tiff import (
    load_store,
    sha256_file,
    validate_pixel_equality,
    write_ome_tiff,
)


FIB_CHANNEL_MAP = {
    "DAPI": "DAPI",
    "520": "PDGFRα",
    "570": "P63",
    "620": "KLF4",
    "690": "KRT14",
}

KC_CHANNEL_MAP = {
    "DAPI": "DAPI",
    "520": "MMP3",
    "570": "SOX9",
    "690": "KRT6B",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_remove_store(store: Path, work_dir: Path) -> None:
    resolved = store.resolve()
    if not resolved.is_relative_to(work_dir.resolve()):
        raise RuntimeError(f"Refusing to remove temporary store outside work directory: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


def audit_existing(store: Path, output: Path) -> dict:
    _, _, arrays, channels = load_store(store)
    report = {
        "output_path": str(output),
        "output_bytes": output.stat().st_size,
        "shape_cyx": list(arrays[0].shape),
        "dtype": str(arrays[0].dtype),
        "channel_names": channels,
        "pyramid_levels": len(arrays),
        "compression": "ADOBE_DEFLATE with horizontal predictor",
    }
    report["validation"] = validate_pixel_equality(store, output)
    report["output_sha256"] = sha256_file(output)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--kfb2zarr", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--name-contains", default="")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    work_dir = output_dir / ".work"
    audit_dir = output_dir / "00_conversion_audit"
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    inputs = sorted(
        (
            item
            for item in input_dir.glob("*.kfbf")
            if not args.name_contains or args.name_contains.casefold() in item.name.casefold()
        ),
        key=lambda item: item.name.casefold(),
    )
    if not inputs:
        raise FileNotFoundError(f"No .kfbf files found in {input_dir}")

    rows: list[dict] = []
    for index, source in enumerate(inputs, start=1):
        started = time.time()
        channel_map = KC_CHANNEL_MAP if source.name.casefold().startswith("kc_mig") else FIB_CHANNEL_MAP
        output = output_dir / f"{source.stem}.ome.tiff"
        store = work_dir / f"slide_{index:03d}.ome.zarr"
        audit_path = audit_dir / f"slide_{index:03d}.audit.json"
        row = {
            "index": index,
            "source_path": str(source),
            "source_bytes": source.stat().st_size,
            "source_mtime_utc": datetime.fromtimestamp(source.stat().st_mtime, timezone.utc).isoformat(),
            "source_sha256": "",
            "output_path": str(output),
            "output_bytes": "",
            "output_sha256": "",
            "status": "running",
            "error": "",
        }
        print(f"[{index}/{len(inputs)}] HASH {source.name}", flush=True)
        try:
            row["source_sha256"] = sha256_file(source)
            if not store.exists():
                print(f"[{index}/{len(inputs)}] DECODE {source.name}", flush=True)
                subprocess.run(
                    [str(args.kfb2zarr), "--threads", str(args.threads), str(source), str(store)],
                    check=True,
                )
            if output.exists():
                print(f"[{index}/{len(inputs)}] REVALIDATE {output.name}", flush=True)
                report = audit_existing(store, output)
            else:
                print(f"[{index}/{len(inputs)}] WRITE {output.name}", flush=True)
                report = write_ome_tiff(store, output)
                print(f"[{index}/{len(inputs)}] VALIDATE {output.name}", flush=True)
                report["validation"] = validate_pixel_equality(store, output)
                report["output_sha256"] = sha256_file(output)

            source_after = source.stat()
            if source_after.st_size != row["source_bytes"]:
                raise RuntimeError("Source size changed during conversion")
            if datetime.fromtimestamp(source_after.st_mtime, timezone.utc).isoformat() != row["source_mtime_utc"]:
                raise RuntimeError("Source modification time changed during conversion")

            row.update(
                {
                    "output_bytes": report["output_bytes"],
                    "output_sha256": report["output_sha256"],
                    "shape_cyx": "x".join(map(str, report["shape_cyx"])),
                    "dtype": report["dtype"],
                    "source_channel_names": "|".join(report["channel_names"]),
                    "marker_mapping": "|".join(
                        f"{name}={channel_map.get(name, 'unmapped')}" for name in report["channel_names"]
                    ),
                    "pyramid_levels": report["pyramid_levels"],
                    "physical_size_x_um": report.get("physical_size_x_um", ""),
                    "physical_size_y_um": report.get("physical_size_y_um", ""),
                    "compression": report["compression"],
                    "pixel_blocks_compared": report["validation"]["pixel_blocks_compared"],
                    "validation_status": report["validation"]["status"],
                    "status": "pass",
                }
            )
            full_report = {
                "created_utc": utc_now(),
                "source": {
                    "path": str(source),
                    "bytes": row["source_bytes"],
                    "mtime_utc": row["source_mtime_utc"],
                    "sha256": row["source_sha256"],
                },
                "output": report,
                "channel_marker_mapping": channel_map,
                "source_preservation_check": "size_and_mtime_unchanged",
            }
            audit_path.write_text(json.dumps(full_report, ensure_ascii=False, indent=2), encoding="utf-8")
            safe_remove_store(store, work_dir)
            print(
                f"[{index}/{len(inputs)}] PASS {output.name} "
                f"({report['output_bytes'] / 1024**2:.1f} MiB, {time.time() - started:.1f}s)",
                flush=True,
            )
        except Exception as exc:
            row["status"] = "failed"
            row["error"] = repr(exc)
            print(f"[{index}/{len(inputs)}] FAILED {source.name}: {exc}", file=sys.stderr, flush=True)
        row["elapsed_seconds"] = round(time.time() - started, 3)
        rows.append(row)

        # Checkpoint after every slide.
        (audit_dir / "batch_manifest.json").write_text(
            json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    fieldnames = sorted({key for row in rows for key in row})
    with (audit_dir / "batch_manifest.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "created_utc": utc_now(),
        "input_directory": str(input_dir),
        "output_directory": str(output_dir),
        "source_count": len(inputs),
        "passed": sum(row["status"] == "pass" for row in rows),
        "failed": sum(row["status"] == "failed" for row in rows),
        "lossless_definition": "No new lossy encoding; TIFF uses Deflate compression. Source KFBF JPEG encoding is preserved only after decoding to uint8 pixels.",
        "imagej_opening": "Use Fiji/ImageJ > Plugins > Bio-Formats > Bio-Formats Importer and enable virtual stack for large slides.",
    }
    (audit_dir / "conversion_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return 0 if summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
