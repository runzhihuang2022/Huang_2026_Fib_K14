from __future__ import annotations

from pathlib import Path
import json
import h5py
import numpy as np


PATHS = {
    "scrna": Path(r"data/external/local_decd995edda5"),
    "stereo": Path(r"data/external/local_b3c2bd3ee952"),
    "visium": Path(r"data/external/local_a3b4e565693f"),
}


def categorical(group: h5py.Group) -> list[str]:
    categories = [
        item.decode("utf-8") if isinstance(item, bytes) else str(item)
        for item in group["categories"][:]
    ]
    codes = np.asarray(group["codes"], dtype=int)
    return sorted({categories[code] for code in codes if code >= 0})


report = {}
for name, file_path in PATHS.items():
    with h5py.File(file_path, "r") as handle:
        obs = handle["obs"]
        entry = {
            "path": str(file_path),
            "shape": [int(handle["X"].attrs.get("shape", [len(obs["_index"]), 0])[0]), int(handle["X"].attrs.get("shape", [0, 0])[1])],
            "obs_fields": sorted(obs.keys()),
            "obsm_fields": sorted(handle.get("obsm", {}).keys()),
            "categorical_examples": {},
        }
        for field in ("sample_id", "samples", "orig.ident", "batch", "library_id"):
            if field in obs and isinstance(obs[field], h5py.Group):
                entry["categorical_examples"][field] = categorical(obs[field])
        report[name] = entry

print(json.dumps(report, indent=2, ensure_ascii=False))
