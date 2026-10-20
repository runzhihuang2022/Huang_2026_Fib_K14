from __future__ import annotations

import json
from pathlib import Path

import anndata as ad


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)

H5AD = Path(r"data/external/local_decd995edda5")


def main() -> None:
    adata = ad.read_h5ad(H5AD, backed="r")
    report = {
        "source": str(H5AD),
        "shape": list(adata.shape),
        "obs_columns": list(map(str, adata.obs.columns)),
        "var_columns": list(map(str, adata.var.columns)),
        "layers": list(map(str, adata.layers.keys())),
        "obsm": list(map(str, adata.obsm.keys())),
        "uns": list(map(str, adata.uns.keys())),
        "raw_shape": list(adata.raw.shape) if adata.raw is not None else None,
        "candidate_obs_values": {},
    }
    for column in (
        "sub_labels",
        "bulk_labels",
        "samples",
        "sample",
        "orig.ident",
        "patient",
        "patient_id",
        "Time",
        "Time_category",
        "Grade",
    ):
        if column in adata.obs:
            values = adata.obs[column].astype(str)
            report["candidate_obs_values"][column] = {
                "n_unique": int(values.nunique(dropna=False)),
                "head": values.value_counts(dropna=False).head(60).to_dict(),
            }
    (OUT / "human_scrna_inventory.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
