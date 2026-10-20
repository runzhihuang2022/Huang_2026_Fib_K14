#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd


LOCKED_SHA256 = "7f7410d1815c1de832d482208cc0bbfdc8ef5e4a7169b34b4007b1075b0b3133"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def copy_tree(source: Path, target: Path) -> None:
    for path in source.rglob("*"):
        if not path.is_file():
            continue
        out = target / path.relative_to(source)
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, out)


def display_gene(symbol: str, width: int = 9) -> str:
    symbol = str(symbol)
    if len(symbol) <= width:
        return symbol
    split = (len(symbol) + 1) // 2
    return symbol[:split] + "\n" + symbol[split:]


def display_term(term: str, width: int = 34) -> str:
    clean = re.sub(r"\s*\(GO:\d+\)\s*$", "", str(term))
    return textwrap.shorten(clean, width=width, placeholder="…")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--only", action="append", default=[], help="Optional FAMILY/SUBTYPE selector; repeatable")
    args = parser.parse_args()
    root = args.project.resolve()
    renderer = root / "00_manifest/locked_renderer/render_top5_tf_network_go_reference_v26.py"
    if sha256(renderer) != LOCKED_SHA256:
        raise RuntimeError("locked TF network/GO renderer SHA-256 mismatch")

    top = pd.read_csv(root / "03_analysis/top5_TFs_by_family_state.tsv", sep="\t")
    top = top[top.annotation_track.eq("sub_labels")].copy()
    all_top = top.copy()
    if args.only:
        selectors = {tuple(item.split("/", 1)) for item in args.only}
        top = top[[((str(f), str(s)) in selectors) for f, s in zip(top.family, top.display_level)]].copy()
    edges = pd.read_csv(root / "03_analysis/top5_TF_target_edges_signed.tsv", sep="\t")
    edges = edges[edges.annotation_track.eq("sub_labels")].copy()
    go = pd.read_csv(root / "03_analysis/functional_enrichment_v02/TF_instance_GO_BP_all.tsv", sep="\t")

    edge_out = edges.rename(columns={
        "TF": "tf_symbol", "target_gene": "target_symbol",
        "motif_edge_weight": "edge_weight",
    })
    edge_out["edge_weight"] = edge_out.edge_weight.astype(float).abs()
    wanted = ["family", "display_level", "rank", "tf_symbol", "RSS", "target_symbol",
              "direction", "edge_weight", "expression_correlation"]
    edge_out = edge_out[wanted].drop_duplicates()

    go_out = go.rename(columns={
        "TF": "tf_symbol", "Description": "Term", "p.adjust": "Adjusted P-value",
        "pvalue": "P-value",
    }).copy()
    go_out["Overlap"] = go_out.Count.astype(int).astype(str) + "/" + go_out.mapped_target_count.astype(int).astype(str)
    go_wanted = ["family", "subtype", "rank", "tf_symbol", "Term", "Adjusted P-value", "P-value", "Overlap"]
    go_out = go_out[go_wanted].drop_duplicates()

    table_dir = root / "05_tables/tf_network_enrichment"
    table_dir.mkdir(parents=True, exist_ok=True)
    edge_out.to_csv(table_dir / "Top5_TF_all_signed_target_genelist.tsv.gz", sep="\t", index=False, compression="gzip")
    go_out.to_csv(table_dir / "Top5_TF_GO_BP_enrichment_expression_background_all.tsv.gz", sep="\t", index=False, compression="gzip")

    final_figures = root / "04_figures/tf_network_enrichment_reference"
    if final_figures.exists() and not args.only:
        shutil.rmtree(final_figures)
    final_figures.mkdir(parents=True, exist_ok=True)
    jobs = root / "03_analysis/locked_v27_network_jobs"
    if jobs.exists():
        shutil.rmtree(jobs)
    records = []
    gene_label_ledger = []
    term_label_ledger = []
    for (family, subtype), block in top.groupby(["family", "display_level"], sort=False):
        job = jobs / f"{safe(family)}__{safe(subtype)}"
        job_table = job / "05_tables/tf_network_enrichment"
        (job / "10_qa").mkdir(parents=True, exist_ok=True)
        job_table.mkdir(parents=True, exist_ok=True)
        tfs = set(block.TF.astype(str))
        one_edges = edge_out[(edge_out.family.eq(family)) & (edge_out.display_level.astype(str).eq(str(subtype))) & edge_out.tf_symbol.isin(tfs)]
        one_go = go_out[(go_out.family.eq(family)) & (go_out.subtype.astype(str).eq(str(subtype))) & go_out.tf_symbol.isin(tfs)]
        job_edges = one_edges.copy()
        job_edges["target_symbol_original"] = job_edges.target_symbol.astype(str)
        job_edges["target_symbol"] = job_edges.target_symbol.astype(str).map(display_gene)
        for row in job_edges[job_edges.target_symbol.ne(job_edges.target_symbol_original)].itertuples(index=False):
            gene_label_ledger.append({"family": family, "subtype": subtype, "tf_symbol": row.tf_symbol, "original": row.target_symbol_original, "display": row.target_symbol, "rule": "line wrap only; exact symbol retained in master table"})
        job_edges.drop(columns=["target_symbol_original"]).to_csv(job_table / "Top5_TF_all_signed_target_genelist.tsv.gz", sep="\t", index=False, compression="gzip")
        job_go = one_go[["tf_symbol", "Term", "Adjusted P-value", "P-value", "Overlap"]].copy()
        job_go["Term_original"] = job_go.Term.astype(str)
        job_go["Term"] = job_go.Term.astype(str).map(display_term)
        for row in job_go[job_go.Term.ne(job_go.Term_original)].itertuples(index=False):
            term_label_ledger.append({"family": family, "subtype": subtype, "tf_symbol": row.tf_symbol, "original": row.Term_original, "display": row.Term, "rule": "display-only deterministic shortening to at most 34 characters; exact term retained in master table"})
        job_go.drop(columns=["Term_original"]).to_csv(
            job_table / "Top5_TF_GO_BP_enrichment_expression_background_all.tsv.gz",
            sep="\t", index=False, compression="gzip"
        )
        proc = subprocess.run([sys.executable, str(renderer), "--project", str(job)], text=True, capture_output=True)
        if proc.returncode:
            raise RuntimeError(f"locked renderer failed for {family}/{subtype}: {proc.stderr}")
        generated = job / "04_figures/tf_network_enrichment_reference"
        copy_tree(generated, final_figures)
        for row in block.itertuples(index=False):
            folder = final_figures / str(family) / safe(str(subtype))
            stem = f"{safe(str(subtype))}_{safe(str(row.TF))}"
            net = folder / f"{stem}_signed_network_reference_v26.pdf"
            gofig = folder / f"{stem}_GO_BP_barplot_reference_v26.pdf"
            records.append({
                "family": family, "subtype": subtype, "rank": int(row.rank),
                "tf_symbol": row.TF, "n_targets": int(len(one_edges[one_edges.tf_symbol.eq(row.TF)])),
                "n_go_terms": int(len(one_go[one_go.tf_symbol.eq(row.TF)])),
                "network_pdf": str(net.relative_to(root)), "go_pdf": str(gofig.relative_to(root)),
                "status": "PASS" if net.exists() and gofig.exists() else "FAIL",
            })

    new_audit = pd.DataFrame(records).sort_values(["family", "subtype", "rank"])
    prior_path = table_dir / "Top5_TF_network_GO_instance_audit.tsv"
    if args.only and prior_path.exists():
        prior = pd.read_csv(prior_path, sep="\t")
        replaced = {(str(f), str(s)) for f, s in zip(new_audit.family, new_audit.subtype)}
        prior = prior[[((str(f), str(s)) not in replaced) for f, s in zip(prior.family, prior.subtype)]]
        audit = pd.concat([prior, new_audit], ignore_index=True).sort_values(["family", "subtype", "rank"])
    else:
        audit = new_audit
    audit.to_csv(table_dir / "Top5_TF_network_GO_instance_audit.tsv", sep="\t", index=False)
    for rows, name in [(gene_label_ledger, "NETWORK_DISPLAY_GENE_LABEL_LEDGER.tsv"), (term_label_ledger, "GO_DISPLAY_TERM_LABEL_LEDGER.tsv")]:
        fresh = pd.DataFrame(rows)
        ledger_path = table_dir / name
        if args.only and ledger_path.exists():
            prior = pd.read_csv(ledger_path, sep="\t")
            replaced = {(str(f), str(s)) for f, s in zip(top.family, top.display_level)}
            prior = prior[[((str(f), str(s)) not in replaced) for f, s in zip(prior.family, prior.subtype)]]
            fresh = pd.concat([prior, fresh], ignore_index=True)
        fresh.to_csv(ledger_path, sep="\t", index=False)
    report = {
        "status": "PASS" if len(audit) == len(all_top) and audit.status.eq("PASS").all() else "FAIL",
        "renderer": str(renderer.relative_to(root)), "renderer_sha256": sha256(renderer),
        "execution_mode": "unchanged locked renderer executed in subtype-isolated input adapters",
        "instances": int(len(audit)), "expected_instances": int(len(all_top)),
        "network_pdf": len(list(final_figures.rglob("*_signed_network_reference_v26.pdf"))),
        "go_pdf": len(list(final_figures.rglob("*_GO_BP_barplot_reference_v26.pdf"))),
    }
    (root / "10_qa/TF_NETWORK_GO_LOCKED_V27_EXECUTION_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
