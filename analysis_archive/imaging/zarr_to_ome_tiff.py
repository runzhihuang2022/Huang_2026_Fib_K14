#!/usr/bin/env python3
"""Stream a multichannel OME-Zarr pyramid to lossless pyramidal OME-TIFF."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Iterator

import numpy as np
import tifffile
import zarr


TILE_SIZE = 256


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_store(store: Path):
    attrs = json.loads((store / ".zattrs").read_text(encoding="utf-8"))
    group = zarr.open_group(str(store), mode="r")
    multiscale = attrs["multiscales"][0]
    datasets = multiscale["datasets"]
    arrays = [group[str(item["path"])] for item in datasets]
    channels = [str(item.get("label", f"Channel {i + 1}")) for i, item in enumerate(attrs["omero"]["channels"])]
    if not arrays:
        raise ValueError("OME-Zarr contains no pyramid datasets")
    if any(array.ndim != 3 for array in arrays):
        raise ValueError("Expected every pyramid level to have C,Y,X dimensions")
    if any(array.shape[0] != len(channels) for array in arrays):
        raise ValueError("Channel count differs between OME metadata and pixel arrays")
    if any(array.dtype != np.dtype("uint8") for array in arrays):
        raise ValueError("This workflow is locked to the source KFBF uint8 pixel type")
    return attrs, datasets, arrays, channels


def tile_iterator(array, tile_size: int = TILE_SIZE) -> Iterator[np.ndarray]:
    channels, height, width = array.shape
    for channel in range(channels):
        for y0 in range(0, height, tile_size):
            for x0 in range(0, width, tile_size):
                yield np.ascontiguousarray(
                    array[channel, y0 : y0 + tile_size, x0 : x0 + tile_size]
                )


def level_scale(dataset: dict) -> tuple[float, float]:
    transforms = dataset.get("coordinateTransformations", [])
    scale = next((item["scale"] for item in transforms if item.get("type") == "scale"), None)
    if scale is None or len(scale) != 3:
        raise ValueError("Missing C,Y,X scale transform in OME-Zarr metadata")
    return float(scale[1]), float(scale[2])


def write_ome_tiff(store: Path, output: Path, compression_level: int = 6) -> dict:
    _, datasets, arrays, channel_names = load_store(store)
    base_y_um, base_x_um = level_scale(datasets[0])
    partial = output.with_name(output.name + ".partial")
    if partial.exists():
        partial.unlink()
    output.parent.mkdir(parents=True, exist_ok=True)

    metadata = {
        "axes": "CYX",
        "SignificantBits": 8,
        "PhysicalSizeX": base_x_um,
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": base_y_um,
        "PhysicalSizeYUnit": "µm",
        "Channel": {"Name": channel_names},
        "Description": "Lossless pyramidal OME-TIFF converted from KFBF via OME-Zarr",
        "MapAnnotation": {
            "Namespace": "openmicroscopy.org/PyramidResolution",
            **{
                str(index): f"{array.shape[2]} {array.shape[1]}"
                for index, array in enumerate(arrays[1:], start=1)
            },
        },
    }
    options = {
        "photometric": "minisblack",
        "tile": (TILE_SIZE, TILE_SIZE),
        "compression": "zlib",
        "compressionargs": {"level": compression_level},
        "predictor": True,
        "resolutionunit": "CENTIMETER",
        "maxworkers": min(4, os.cpu_count() or 1),
        "software": "kfbf-to-ome-tiff (kfb2zarr + tifffile)",
    }

    started = time.time()
    with tifffile.TiffWriter(partial, bigtiff=True, ome=True) as writer:
        base = arrays[0]
        writer.write(
            tile_iterator(base),
            shape=base.shape,
            dtype=base.dtype,
            subifds=len(arrays) - 1,
            resolution=(1e4 / base_x_um, 1e4 / base_y_um),
            metadata=metadata,
            **options,
        )
        for dataset, array in zip(datasets[1:], arrays[1:], strict=True):
            y_um, x_um = level_scale(dataset)
            writer.write(
                tile_iterator(array),
                shape=array.shape,
                dtype=array.dtype,
                subfiletype=1,
                resolution=(1e4 / x_um, 1e4 / y_um),
                metadata=None,
                **options,
            )
    partial.replace(output)
    return {
        "output_path": str(output),
        "output_bytes": output.stat().st_size,
        "elapsed_seconds": round(time.time() - started, 3),
        "shape_cyx": list(arrays[0].shape),
        "dtype": str(arrays[0].dtype),
        "channel_names": channel_names,
        "pyramid_levels": len(arrays),
        "physical_size_x_um": base_x_um,
        "physical_size_y_um": base_y_um,
        "compression": f"ADOBE_DEFLATE level {compression_level} with horizontal predictor",
    }


def validation_points(height: int, width: int) -> list[tuple[int, int, int, int]]:
    size = min(64, height, width)
    return [
        (0, 0, size, size),
        (max(0, height // 2 - size // 2), max(0, width // 2 - size // 2), size, size),
        (max(0, height - size), max(0, width - size), size, size),
    ]


def validate_pixel_equality(store: Path, output: Path) -> dict:
    _, _, source_arrays, channel_names = load_store(store)
    with tifffile.TiffFile(output) as tif:
        if not tif.is_ome:
            raise ValueError("Output is not recognized as OME-TIFF")
        series = tif.series[0]
        if series.axes != "CYX":
            raise ValueError(f"Unexpected TIFF axes: {series.axes}")
        if tuple(series.shape) != tuple(source_arrays[0].shape):
            raise ValueError(f"Base shape mismatch: {series.shape} != {source_arrays[0].shape}")
        if len(series.levels) != len(source_arrays):
            raise ValueError(f"Pyramid level mismatch: {len(series.levels)} != {len(source_arrays)}")
        ome_xml = tif.ome_metadata or ""
        missing = [name for name in channel_names if f'Name="{name}"' not in ome_xml]
        if missing:
            raise ValueError(f"Channel names missing from OME-XML: {missing}")

    comparisons = 0
    for level_index in sorted(set([0, len(source_arrays) // 2, len(source_arrays) - 1])):
        source = source_arrays[level_index]
        store_view = tifffile.imread(output, series=0, level=level_index, aszarr=True)
        try:
            target = zarr.open(store_view, mode="r")
            for channel in range(source.shape[0]):
                for y0, x0, height, width in validation_points(source.shape[1], source.shape[2]):
                    expected = np.asarray(source[channel, y0 : y0 + height, x0 : x0 + width])
                    observed = np.asarray(target[channel, y0 : y0 + height, x0 : x0 + width])
                    if not np.array_equal(expected, observed):
                        raise ValueError(
                            f"Pixel mismatch at level={level_index}, channel={channel}, x={x0}, y={y0}"
                        )
                    comparisons += 1
        finally:
            close = getattr(store_view, "close", None)
            if close is not None:
                close()
    return {
        "status": "pass",
        "pixel_blocks_compared": comparisons,
        "levels_sampled": sorted(set([0, len(source_arrays) // 2, len(source_arrays) - 1])),
        "channels_compared": len(channel_names),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("store", type=Path, help="Input .ome.zarr directory")
    parser.add_argument("output", type=Path, help="Output .ome.tiff file")
    parser.add_argument("--compression-level", type=int, default=6, choices=range(0, 10))
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")
    report = write_ome_tiff(args.store.resolve(), args.output.resolve(), args.compression_level)
    report["validation"] = validate_pixel_equality(args.store.resolve(), args.output.resolve())
    report["output_sha256"] = sha256_file(args.output.resolve())
    report_path = args.report or args.output.with_suffix(args.output.suffix + ".audit.json")
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
