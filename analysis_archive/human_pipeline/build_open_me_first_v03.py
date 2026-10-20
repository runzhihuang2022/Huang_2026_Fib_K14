#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path

import pandas as pd


FAMILIES = ["KC", "SAC", "Endo", "Fib", "Myeloid", "Lymphoid"]


def safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def url(path: Path, html_dir: Path) -> str:
    return "../" + path.relative_to(html_dir.parent).as_posix().replace("%", "%25").replace("#", "%23")


def card(pdf: Path, html_dir: Path, title: str, note: str = "") -> str:
    png = pdf.with_suffix(".png")
    if not pdf.exists():
        return f'<figure class="blocked"><b>{html.escape(title)}</b><p>NOT_ESTIMABLE or missing registered panel</p></figure>'
    preview = f'<a href="{url(pdf, html_dir)}"><img loading="lazy" src="{url(png, html_dir)}"></a>' if png.exists() else ""
    links = f'<a class="button" href="{url(pdf, html_dir)}">PDF</a>'
    if png.exists(): links += f'<a class="button" href="{url(png, html_dir)}">PNG</a>'
    tiff = pdf.with_suffix(".tiff")
    if tiff.exists(): links += f'<a class="button" href="{url(tiff, html_dir)}">TIFF</a>'
    return f'<figure>{preview}<figcaption><b>{html.escape(title)}</b><br><span class="small">{html.escape(note)}</span><br>{links}</figcaption></figure>'


def file_link(path: Path, html_dir: Path, label: str) -> str:
    return f'<a class="button" href="{url(path, html_dir)}">{html.escape(label)}</a>' if path.exists() else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve(); html_dir = root / "06_html"; html_dir.mkdir(exist_ok=True)
    output = html_dir / "OPEN_ME_FIRST_human_epithelialization_complete_v03.html"
    top = pd.read_csv(root / "03_analysis/top5_TFs_by_family_state.tsv", sep="\t")
    top = top[top.annotation_track.eq("sub_labels")].sort_values(["family", "display_level", "rank"])
    network = pd.read_csv(root / "05_tables/tf_network_enrichment/Top5_TF_network_GO_instance_audit.tsv", sep="\t")
    oracle = pd.read_csv(root / "05_tables/CellOracle_v02/subtype_visualization_registry.tsv", sep="\t")
    monocle = pd.read_csv(root / "05_tables/Monocle3_v02/subtype_visualization_registry.tsv", sep="\t")
    companions = pd.read_csv(root / "05_tables/CellChat_v02/reclustered_umap_visualization_registry.tsv", sep="\t")
    chords = pd.read_csv(root / "05_tables/cellchat_default_chord/CELLCHAT_DEFAULT_CHORD_INSTANCE_AUDIT.tsv", sep="\t")

    css = """body{font-family:Arial,sans-serif;margin:0;color:#20242a;background:#f5f7fa}header,main{max-width:1500px;margin:auto;padding:18px}header{background:#14213d;color:white;max-width:none}h1{margin:0;font-size:24px}h2{border-bottom:2px solid #d9e1ea;padding-bottom:6px}h3{margin-top:22px}details{background:white;border:1px solid #d9e1ea;border-radius:8px;margin:12px 0;padding:10px}summary{font-weight:bold;cursor:pointer}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(245px,1fr));gap:12px}figure{margin:0;background:white;border:1px solid #e0e4e8;padding:8px;border-radius:6px}img{width:100%;height:auto;display:block}figcaption{font-size:12px;margin-top:5px}.button{display:inline-block;padding:5px 9px;margin:3px;background:#0b6e99;color:white;text-decoration:none;border-radius:4px;font-size:12px}.notice{padding:10px;border-left:4px solid #c69214;background:#fff8df;margin:10px 0}.small{font-size:12px;color:#56606b}.blocked{background:#fff1f1}nav{position:sticky;top:0;background:white;padding:8px;z-index:4;border-bottom:1px solid #ccd}nav a{margin-right:8px}"""
    parts = [f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Human epithelialization complete analysis v03</title><style>{css}</style></head><body><header><h1>人创面上皮化单细胞完整亚型分析 v03</h1><p>Figure 2 + Monocle3 + pySCENIC + 锁定版TF network/GO + CellOracle + 锁定版CellChat Top8显著LR</p></header><main><div class="notice"><b>代码边界：</b>network/GO和CellChat采用用户指定v27锁定渲染器；只适配人上皮化输入、sub_labels、Normal/SPTDI/DPTDI及Mode 1–4分组。正负边表示表达相关方向；Monocle3、SCENIC、CellOracle和CellChat均不是因果验证。</div><nav><a href="#lineages">完整亚型</a><a href="#networks">TF network/GO</a><a href="#cellchat">CellChat</a><a href="#qa">代码与QA</a></nav><h2 id="lineages">一、Figure 2、Monocle3、pySCENIC与CellOracle</h2>''']
    for family in FAMILIES:
        parts.append(f'<h3>{family}</h3><div class="grid">')
        for path in [root / f"04_figures/Shared_geometry_v02/{family}_sub_labels_and_reclustered_leiden.pdf", root / f"04_figures/Figure2_v02_reclustered/Figure2_{family}_subtypes_Mode1-4_A-E.pdf"]:
            if path.exists(): parts.append(card(path, html_dir, path.stem, "frozen lineage-reclustered UMAP"))
        parts.append('</div>')
        subtypes = top[top.family.eq(family)].display_level.astype(str).drop_duplicates().tolist()
        for subtype in subtypes:
            parts.append(f'<details><summary>{html.escape(subtype)} — Monocle3、pySCENIC Top5、CellOracle</summary><div class="grid">')
            mono = root / f"04_figures/{family}/Monocle3_v02/{safe(subtype)}/{safe(subtype)}_Monocle3_reclustered_UMAP.pdf"
            parts.append(card(mono, html_dir, f"{subtype} Monocle3", "same frozen lineage geometry"))
            for row in top[(top.family.eq(family)) & (top.display_level.astype(str).eq(subtype))].itertuples(index=False):
                scenic = root / f"04_figures/{family}/SCENIC_top5/{safe(subtype)}__rank{int(row.rank)}__{safe(row.TF)}.pdf"
                parts.append(card(scenic, html_dir, f"{subtype} · {row.TF} pySCENIC", "RSS + subtype UMAP + AUCell"))
            oblock = oracle[(oracle.family.eq(family)) & (oracle.subtype.astype(str).eq(subtype))]
            if len(oblock) and oblock.iloc[0].status == "PASS":
                for condition in ["KO", "OE"]:
                    parts.append(card(root / str(oblock.iloc[0][f"{condition}_figure"]), html_dir, f"{subtype} · {oblock.iloc[0].TF} {condition}", "CellOracle computational prediction"))
            else:
                reason = str(oblock.iloc[0].reason) if len(oblock) else "missing registry"
                parts.append(f'<figure class="blocked"><b>{html.escape(subtype)} CellOracle</b><p>NOT_ESTIMABLE: {html.escape(reason)}</p></figure>')
            parts.append('</div></details>')

    parts.append('<h2 id="networks">二、锁定v27代码重绘的Top5 TF network与GO</h2><p>45个亚型×Top5 TF，共225个network、225个GO和225个network+GO配对页。下列结果由SHA锁定渲染器生成，不使用v02通用替代绘图器。</p>')
    fig_root = root / "04_figures/tf_network_enrichment_reference"
    for (family, subtype), block in network.groupby(["family", "subtype"], sort=False):
        folder = fig_root / str(family) / safe(str(subtype))
        parts.append(f'<details><summary>{html.escape(str(family))} · {html.escape(str(subtype))}</summary>')
        parts.append(file_link(folder / f"{safe(str(subtype))}_Top5_TF_signed_networks_reference_v26.pdf", html_dir, "Top5 network multipage PDF"))
        parts.append(file_link(folder / f"{safe(str(subtype))}_Top5_TF_GO_BP_barplots_reference_v26.pdf", html_dir, "Top5 GO multipage PDF"))
        parts.append('<div class="grid">')
        for row in block.sort_values("rank").itertuples(index=False):
            prefix = f"{safe(str(subtype))}_{safe(str(row.tf_symbol))}"
            parts.append(card(folder / f"{prefix}_signed_network_reference_v26.pdf", html_dir, f"{subtype} · {row.tf_symbol} network"))
            parts.append(card(folder / f"{prefix}_GO_BP_barplot_reference_v26.pdf", html_dir, f"{subtype} · {row.tf_symbol} GO BP"))
            parts.append(card(folder / f"{prefix}_network_GO_paired_150x70mm.pdf", html_dir, f"{subtype} · {row.tf_symbol} paired network+GO"))
        parts.append('</div></details>')

    parts.append('<h2 id="cellchat">三、锁定v27代码重绘的CellChat</h2><p>CellChat主体图由7个正式RDS重新渲染；每个可估计界面只展示确定性排序的前8对显著LR，不补入不显著结果。UMAP companion沿用对应谱系冻结坐标。</p><h3>Top8显著LR默认弦图</h3><div class="grid">')
    chord_root = root / "04_figures/cellchat_default_chord"
    for row in chords.itertuples(index=False):
        stem = chord_root / f"{row.group}_{row.interface}_CellChat_default_LR_chord.pdf"
        if row.status == "PASS": parts.append(card(stem, html_dir, f"{row.group} · {row.interface}", f"{row.n_display_lr} significant LR pairs"))
        else: parts.append(f'<figure class="blocked"><b>{html.escape(row.group)} · {html.escape(row.interface)}</b><p>NOT_ESTIMABLE: {html.escape(str(row.reason))}</p></figure>')
    parts.append('</div><h3>CellChat network、中心性、LR chord、violin和GO</h3><div class="grid">')
    targeted = root / "04_figures/cellchat_targeted_locked_v27"
    for pdf in sorted(targeted.glob("*.pdf")):
        parts.append(card(pdf, html_dir, pdf.stem, "locked v27 CellChat renderer"))
    for pdf in sorted((root / "04_figures/cellchat_lr_enrichment_locked_v27").glob("*.pdf")):
        parts.append(card(pdf, html_dir, pdf.stem, "GO BP over CellChat signaling background"))
    parts.append('</div><h3>谱系重聚类UMAP companion</h3><div class="grid">')
    for row in companions.itertuples(index=False):
        if row.status != "PASS": continue
        for col, direction in [("outgoing_figure", "outgoing"), ("incoming_figure", "incoming")]:
            parts.append(card(root / str(getattr(row, col)), html_dir, f"{row.context} · {row.interface} · {row.family} {direction}", "subtype-level CellChat summary painted onto frozen UMAP"))
    parts.append('</div><h2 id="qa">四、数据表、代码与QA</h2>')
    links = [
        (root / "00_manifest/LOCKED_RENDERER_BASELINE_SHA256.tsv", "锁定代码SHA-256"),
        (root / "00_manifest/LOCKED_RENDERER_ALLOWED_DIFF.tsv", "允许改动清单"),
        (root / "05_tables/tf_network_enrichment/Top5_TF_all_signed_target_genelist.tsv.gz", "全部TF-target表"),
        (root / "05_tables/tf_network_enrichment/Top5_TF_GO_BP_enrichment_expression_background_all.tsv.gz", "全部GO表"),
        (root / "05_tables/cellchat_default_chord/ALL_DISPLAYED_DEFAULT_CHORD_LR.tsv", "CellChat展示LR表"),
        (root / "05_tables/cellchat_default_chord/CELLCHAT_DEFAULT_CHORD_INSTANCE_AUDIT.tsv", "CellChat实例审计"),
        (root / "10_qa/FINAL_V03_LOCKED_RENDER_QA.json", "最终QA"),
        (root / "08_code/run_locked_v27_networks.py", "network输入适配与锁定调用"),
        (root / "08_code/render_cellchat_default_chord_v26.R", "CellChat锁定弦图代码"),
    ]
    parts.extend(file_link(path, html_dir, label) for path, label in links)
    parts.append('</main></body></html>')
    output.write_text("".join(parts), encoding="utf-8")
    print(json.dumps({"html": str(output), "bytes": output.stat().st_size}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
