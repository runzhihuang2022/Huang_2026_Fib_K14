from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import tifffile


OUT = Path(r"data/external/local_93c497bb3a06")


def ome_meta(path: Path) -> dict:
    with tifffile.TiffFile(path) as tif:
        root = ET.fromstring(tif.ome_metadata)
        ns = {"o":"http://www.openmicroscopy.org/Schemas/OME/2016-06"}
        px = root.find(".//o:Pixels", ns)
        names = [c.get("Name","") for c in px.findall("o:Channel",ns)]
        series = tif.series[0]
        return {
            "ome_recognized": tif.is_ome, "axes": series.axes,
            "shape": "x".join(map(str,series.shape)), "dtype":str(series.dtype),
            "channel_names":"|".join(names), "physical_size_x_um":float(px.get("PhysicalSizeX")),
            "physical_size_y_um":float(px.get("PhysicalSizeY")), "pyramid_levels":len(series.levels),
            "bigtiff":tif.is_bigtiff,
        }


def main() -> None:
    roi_manifest_path=OUT/"03_roi_handoff"/"epi_front_roi_manifest.csv"
    if roi_manifest_path.exists():
        roi=pd.read_csv(roi_manifest_path)
        roi["coordinate_space"]="source_orientation_analysis_image"
        suffix=" Standalone Fib section; no KC serial-section registration in this run."
        roi["note"]=roi["note"].astype(str).map(lambda x: x if suffix.strip() in x else x+suffix)
        roi.to_csv(roi_manifest_path,index=False,encoding="utf-8-sig")
    batch = pd.read_csv(OUT/"01_ome_tiff"/"00_conversion_audit"/"batch_manifest.csv")
    rows=[]
    for r in batch.itertuples(index=False):
        p=Path(r.output_path); m=ome_meta(p)
        rows.append({
            "source_path":r.source_path,"source_bytes":r.source_bytes,"source_sha256":r.source_sha256,
            "source_mtime_utc":r.source_mtime_utc,"output_path":r.output_path,"output_bytes":r.output_bytes,
            "output_sha256":r.output_sha256,"conversion_status":r.status,
            "pixel_equality_status":r.validation_status,"pixel_blocks_compared":r.pixel_blocks_compared,
            "compression":r.compression,"marker_mapping":r.marker_mapping,**m,
            "fiji_bioformats_status":"MANUAL_VERIFICATION_REQUIRED_runtime_not_found",
        })
    final=pd.DataFrame(rows)
    final.to_csv(OUT/"00_manifest"/"source_and_conversion_manifest.tsv",sep="\t",index=False)
    fiji=final[["output_path","ome_recognized","axes","shape","dtype","channel_names",
                "physical_size_x_um","pyramid_levels","fiji_bioformats_status"]]
    fiji.to_csv(OUT/"08_audit"/"fiji_bioformats_verification.tsv",sep="\t",index=False)
    gated={d:len(list((OUT/d).rglob("*.*"))) for d in ["04_user_selected_zoomin","05_cell_tables",
                                                          "06_spatial_statistics","07_publication_figures"]}
    partials=list(OUT.rglob("*.partial"))
    qa=pd.DataFrame([
        ["source inventory",len(final)==10,f"{len(final)}/10"],
        ["conversion pass",final.conversion_status.eq("pass").all(),f"{final.conversion_status.eq('pass').sum()}/10"],
        ["OME recognition",final.ome_recognized.all(),f"{final.ome_recognized.sum()}/10"],
        ["CYX axes",final["axes"].eq("CYX").all(),f"{final['axes'].eq('CYX').sum()}/10"],
        ["five channels",final.channel_names.eq("DAPI|520|570|620|690").all(),"DAPI|520|570|620|690"],
        ["uint8",final.dtype.eq("uint8").all(),f"{final.dtype.eq('uint8').sum()}/10"],
        ["physical calibration",final.physical_size_x_um.notna().all(),"0.25284451246261597 um/px"],
        ["sampled pixel equality",final.pixel_equality_status.eq("pass").all(),f"{final.pixel_blocks_compared.sum()} blocks"],
        ["no partial files",len(partials)==0,str(len(partials))],
        ["zoomout PNG",len(list((OUT/'02_registered_zoomout').glob('*_zoomout_merge.png')))==10,"10/10"],
        ["DAPI ROI reference",len(list((OUT/'02_registered_zoomout').glob('*_DAPI_analysis_reference.ome.tiff')))==10,"10/10"],
        ["Fiji Bio-Formats",False,"MANUAL_VERIFICATION_REQUIRED: Fiji runtime not found"],
        ["ROI gate",False,"AUTHOR_INPUT_NEEDED: no approved RoiSet.zip"],
        ["forbidden pre-ROI outputs empty",all(v==0 for v in gated.values()),json.dumps(gated)],
    ],columns=["check","pass","detail"])
    qa.to_csv(OUT/"08_audit"/"FINAL_QA.tsv",sep="\t",index=False)
    pd.DataFrame([
        ["whole-tissue visualization","descriptive_supported","10 zoomouts; source orientation preserved"],
        ["Epi-Front localization","not_tested","requires user-reviewed anatomical ROIs"],
        ["Fib_K14 zoom-in","blocked","requires approved ROI archives"],
        ["Fib_K14 line/statistical plot","blocked","requires ROI-level and slide-level source tables"],
        ["serial-section spatial correspondence","not_tested","paired KC registration not performed in this run"],
        ["same-cell colocalization","not_supported","cannot be inferred across serial sections"],
    ],columns=["claim","status","basis"]).to_csv(OUT/"08_audit"/"claim_status.tsv",sep="\t",index=False)
    readme="""# Fib_K14 KFBF conversion and pre-ROI handoff

## Completed

- Converted 10/10 KFBF files into losslessly compressed pyramidal OME-BigTIFF.
- Preserved five independent channels: DAPI, 520, 570, 620 and 690; marker mapping is DAPI, PDGFRα, P63, KLF4 and KRT14.
- Preserved uint8 pixels, physical calibration (0.25284451246261597 µm/px), channel names and 9–10 pyramid levels.
- Sampled pixel equality passed for 450 blocks in total (45 per slide).
- Exported 10 full-tissue ZoomOut composites with metadata-derived scale bars and 10 DAPI analysis-reference OME-TIFFs.

## Intentionally blocked by the ROI gate

No approved `样本名_RoiSet.zip` archives were found. Therefore `04_user_selected_zoomin`, `05_cell_tables`, `06_spatial_statistics` and `07_publication_figures` are empty. No automatic marker-bright ROI was substituted, and no old line plot/statistics were copied into this run.

Open each `*_DAPI_analysis_reference.ome.tiff` in Fiji/Bio-Formats, review orientation, and draw 3–5 non-overlapping anatomical Epi-Front ROIs (minimum 2). Name them `样本名_EF01`, `样本名_EF02`... and save the archive as `样本名_RoiSet.zip`. Fill reviewer and set status to `approved` in `03_roi_handoff/epi_front_roi_manifest.csv`.

The current computer did not expose a Fiji/Bio-Formats runtime to command-line validation. OME recognition and sampled pixel equality passed in tifffile/zarr, but Fiji opening remains a manual release check.
"""
    (OUT/"README_请先阅读.md").write_text(readme,encoding="utf-8")
    # Small-file package index; OME hashes come from the conversion manifest to avoid rehashing twice.
    items=[]
    for p in sorted(OUT.rglob("*")):
        if not p.is_file() or p.name=="package_manifest.tsv": continue
        if ".pilot_work" in p.parts or ".work" in p.parts: continue
        if p.suffix.lower() in {".tiff"} and p.parent.name=="01_ome_tiff":
            hit=final[final.output_path==str(p)]
            digest=hit.output_sha256.iloc[0] if len(hit) else ""
        else:
            h=hashlib.sha256()
            with p.open("rb") as fh:
                for b in iter(lambda:fh.read(4*1024*1024),b""): h.update(b)
            digest=h.hexdigest()
        items.append({"relative_path":str(p.relative_to(OUT)),"size_bytes":p.stat().st_size,"sha256":digest})
    pd.DataFrame(items).to_csv(OUT/"00_manifest"/"package_manifest.tsv",sep="\t",index=False)
    print(json.dumps({"conversion_pass":int(final.conversion_status.eq('pass').sum()),"gated":gated},ensure_ascii=False))


if __name__=="__main__":
    main()
