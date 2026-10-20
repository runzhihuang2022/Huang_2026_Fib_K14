# -*- coding: utf-8 -*-
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfamsa import Aligner, Sequence
from scipy import stats


ROOT = Path(r"data/external/local_a82f375325f6")
FIG = ROOT / "figures"
TAB = ROOT / "tables"
SEQ = ROOT / "resources" / "sequences"
REP = ROOT / "reports"
ATAC = Path(r"data/external/local_44ea9b8fef30")

mpl.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
})

TARGETS = {
    "KLF4": [
        "homo_sapiens", "rattus_norvegicus", "mus_musculus",
        "cervus_hanglu_yarkandensis", "cavia_porcellus",
        "heterocephalus_glaber_female", "moschus_moschiferus",
    ],
    "IGFL1": [
        "homo_sapiens", "cervus_hanglu_yarkandensis", "cavia_porcellus",
        "heterocephalus_glaber_female", "moschus_moschiferus",
        "bos_taurus", "capra_hircus",
    ],
}

DISPLAY = {
    "homo_sapiens": "Human", "rattus_norvegicus": "Rat",
    "mus_musculus": "Mouse", "cervus_hanglu_yarkandensis": "Cervid",
    "cavia_porcellus": "Guinea pig proxy", "heterocephalus_glaber_female": "Naked mole-rat proxy",
    "moschus_moschiferus": "Musk deer (cervid proxy)", "bos_taurus": "Cattle", "capra_hircus": "Goat",
}


def ensure_dirs() -> None:
    for p in [FIG, TAB, SEQ, REP]:
        p.mkdir(parents=True, exist_ok=True)


def ensembl_homologies(gene: str) -> tuple[list[dict], dict]:
    url = (
        f"https://rest.ensembl.org/homology/symbol/homo_sapiens/{gene}"
        "?content-type=application/json;type=orthologues;sequence=protein"
    )
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    payload = response.json()
    return payload["data"][0]["homologies"], payload


def select_records(gene: str, homologies: list[dict]) -> tuple[list[SeqRecord], pd.DataFrame]:
    by_species: dict[str, list[dict]] = {}
    for h in homologies:
        by_species.setdefault(h["target"]["species"], []).append(h)

    source = homologies[0]["source"]
    source_seq = source["align_seq"].replace("-", "")
    records = [SeqRecord(Seq(source_seq), id="Human", description=f"{gene}|{source['protein_id']}")]
    rows = [{
        "gene": gene, "requested_species": "homo_sapiens", "display_species": "Human",
        "strict_ortholog_present": True, "orthology_type": "reference",
        "gene_id": source["id"], "protein_id": source["protein_id"],
        "percent_identity_to_human": 100.0, "sequence_length_aa": len(source_seq),
        "tree_included": True, "interpretation": "human reference sequence",
    }]

    for species in TARGETS[gene][1:]:
        candidates = by_species.get(species, [])
        if not candidates:
            rows.append({
                "gene": gene, "requested_species": species,
                "display_species": DISPLAY.get(species, species), "strict_ortholog_present": False,
                "orthology_type": "not returned by Ensembl Compara", "gene_id": "", "protein_id": "",
                "percent_identity_to_human": np.nan, "sequence_length_aa": np.nan,
                "tree_included": False, "interpretation": "strict sequence evidence unavailable; do not place on tree",
            })
            continue
        # Prefer one-to-one orthology, then the highest protein identity.
        candidates.sort(key=lambda h: (h["type"] == "ortholog_one2one",
                                       len(h["target"].get("align_seq", "").replace("-", "")),
                                       h["target"].get("perc_id", 0)), reverse=True)
        h = candidates[0]
        target = h["target"]
        seq = target["align_seq"].replace("-", "")
        label = DISPLAY.get(species, species)
        tree_included = not (gene == "KLF4" and len(seq) < 0.65 * len(source_seq))
        if tree_included:
            records.append(SeqRecord(Seq(seq), id=label.replace(" ", "_"), description=f"{gene}|{target['protein_id']}"))
        rows.append({
            "gene": gene, "requested_species": species, "display_species": label,
            "strict_ortholog_present": True, "orthology_type": h["type"],
            "gene_id": target["id"], "protein_id": target["protein_id"],
            "percent_identity_to_human": target.get("perc_id", np.nan), "sequence_length_aa": len(seq),
            "tree_included": tree_included,
            "interpretation": ("strict Ensembl Compara protein ortholog" if tree_included else
                               "strict ortholog record is partial; retained in audit but excluded from full-length tree"),
        })
    return records, pd.DataFrame(rows)


def famsa_align(records: list[SeqRecord]) -> MultipleSeqAlignment:
    aligner = Aligner(guide_tree="upgma")
    seqs = [Sequence(r.id.encode(), bytes(str(r.seq), "ascii")) for r in records]
    aligned = aligner.align(seqs)
    out = []
    for item in aligned:
        out.append(SeqRecord(Seq(item.sequence.decode()), id=item.id.decode(), description=""))
    return MultipleSeqAlignment(out)


def bootstrap_support(aln: MultipleSeqAlignment, tree, n_boot: int = 300, seed: int = 20260713) -> None:
    rng = np.random.default_rng(seed)
    calc = DistanceCalculator("blosum62")
    constructor = DistanceTreeConstructor()
    labels = {x.name for x in tree.get_terminals()}
    counts: Counter[frozenset[str]] = Counter()
    width = aln.get_alignment_length()
    for _ in range(n_boot):
        cols = rng.integers(0, width, width)
        boot = MultipleSeqAlignment([
            SeqRecord(Seq("".join(str(rec.seq)[i] for i in cols)), id=rec.id) for rec in aln
        ])
        bt = constructor.nj(calc.get_distance(boot))
        for clade in bt.get_nonterminals():
            members = frozenset(x.name for x in clade.get_terminals())
            if 1 < len(members) < len(labels):
                counts[members] += 1
    for clade in tree.get_nonterminals():
        members = frozenset(x.name for x in clade.get_terminals())
        if 1 < len(members) < len(labels):
            clade.confidence = round(100 * counts[members] / n_boot)


def phylogeny(gene: str, records: list[SeqRecord]) -> tuple[pd.DataFrame, Path]:
    aln = famsa_align(records)
    fasta = SEQ / f"{gene}_strict_proteins.fasta"
    aligned = SEQ / f"{gene}_strict_proteins_aligned.fasta"
    with fasta.open("w", encoding="ascii") as handle:
        for r in records:
            handle.write(f">{r.id} {r.description}\n{str(r.seq)}\n")
    with aligned.open("w", encoding="ascii") as handle:
        AlignIO.write(aln, handle, "fasta")

    calculator = DistanceCalculator("blosum62")
    dm = calculator.get_distance(aln)
    tree = DistanceTreeConstructor().nj(dm)
    tree.root_at_midpoint()
    bootstrap_support(aln, tree)
    newick = SEQ / f"{gene}_strict_NJ_bootstrap.newick"
    Phylo.write(tree, str(newick), "newick")

    names = list(dm.names)
    ident = pd.DataFrame(index=names, columns=names, dtype=float)
    seqs = {r.id: str(r.seq) for r in aln}
    for a in names:
        for b in names:
            pairs = [(x, y) for x, y in zip(seqs[a], seqs[b]) if x != "-" and y != "-"]
            ident.loc[a, b] = 100 * sum(x == y for x, y in pairs) / max(len(pairs), 1)
    ident.to_csv(TAB / f"{gene}_strict_protein_percent_identity.csv")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), gridspec_kw={"width_ratios": [1.05, 1]})
    Phylo.draw(tree, axes=axes[0], do_show=False, show_confidence=False,
               label_func=lambda c: c.name.replace("_", " ") if c.is_terminal() else None,
               branch_labels=lambda c: int(c.confidence) if c.confidence is not None else None)
    axes[0].set_title(f"{gene} strict protein orthology")
    axes[0].set_xlabel("Amino-acid distance")
    axes[0].set_ylabel("")
    im = axes[1].imshow(ident.to_numpy(float), cmap="Blues", vmin=0, vmax=100, aspect="equal")
    axes[1].set_xticks(range(len(names)), [x.replace("_", " ") for x in names], rotation=45, ha="right")
    axes[1].set_yticks(range(len(names)), [x.replace("_", " ") for x in names])
    for i in range(len(names)):
        for j in range(len(names)):
            axes[1].text(j, i, f"{ident.iloc[i, j]:.0f}", ha="center", va="center", fontsize=7,
                         color="white" if ident.iloc[i, j] > 70 else "black")
    cb = fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    cb.set_label("Protein identity (%)")
    axes[1].set_title("Pairwise identity")
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"01_{gene}_strict_sequence_tree.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return ident, newick


def atac_evidence() -> dict:
    sample = pd.read_csv(ATAC / "KLF4_chromVAR_FibK14_by_sample.csv", skipinitialspace=True)
    peaks = pd.read_csv(ATAC / "fibk14_klf4_target_peak_partial_spearman.csv")
    stage = pd.read_csv(ATAC / "stage4_statistical_summary.csv")
    sample.to_csv(TAB / "KLF4_chromVAR_FibK14_by_sample.csv", index=False)
    peaks.to_csv(TAB / "KLF4_target_peak_partial_spearman.csv", index=False)
    stage.to_csv(TAB / "KLF4_ATAC_stage_summary.csv", index=False)

    peaks["minus_log10_FDR"] = -np.log10(peaks["FDR"].clip(lower=np.finfo(float).tiny))
    peaks["significant"] = peaks["FDR"] < 0.05
    key = {"Krt14", "Krt5", "Krt17", "Epcam", "Tacstd2", "Grhl3", "Areg", "Col1a1", "Col17a1", "Cd34"}
    highlights = peaks.loc[peaks["target_gene"].isin(key)].copy()
    highlights.to_csv(TAB / "KLF4_ATAC_highlighted_target_peaks.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.0), gridspec_kw={"width_ratios": [0.92, 1.35]})
    order = sample.sort_values("KLF4_chromVAR_z.mean")
    colors = np.where(order["KLF4_chromVAR_z.mean"] >= 0, "#D55E00", "#4C78A8")
    axes[0].barh(np.arange(len(order)), order["KLF4_chromVAR_z.mean"], color=colors, edgecolor="#444444", linewidth=.5)
    axes[0].axvline(0, color="#333333", lw=.8)
    axes[0].set_yticks(np.arange(len(order)), [x.replace("230907YYLib_V2_", "").replace("231008YYLib_V2_", "") for x in order["sample_id"]], fontsize=6.5)
    axes[0].set_xlabel("KLF4 motif chromVAR z-score")
    axes[0].set_title("Fib_K14 sample-level accessibility")

    base = peaks.loc[~peaks["significant"]]
    sig = peaks.loc[peaks["significant"]]
    axes[1].scatter(base["partial_spearman"], base["minus_log10_FDR"], s=12, color="#BDBDBD", alpha=.65, linewidths=0)
    axes[1].scatter(sig["partial_spearman"], sig["minus_log10_FDR"], s=18,
                    c=np.where(sig["partial_spearman"] > 0, "#D55E00", "#0072B2"), alpha=.8, linewidths=0)
    axes[1].axhline(-math.log10(.05), color="#555555", ls="--", lw=.8)
    axes[1].axvline(0, color="#555555", lw=.7)
    label_rows = highlights.sort_values("FDR").groupby("target_gene", as_index=False).first()
    for _, row in label_rows.iterrows():
        axes[1].annotate(row["target_gene"], (row["partial_spearman"], row["minus_log10_FDR"]),
                         xytext=(3, 3), textcoords="offset points", fontsize=7)
    axes[1].set_xlabel("Partial Spearman rho\n(KLF4 motif accessibility vs target peak)")
    axes[1].set_ylabel("-log10(FDR)")
    axes[1].set_title("KLF4 motif-linked target peaks")
    fig.tight_layout(w_pad=1.2)
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"02_KLF4_ATAC_motif_validation.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {
        "fibk14_samples": int(len(sample)),
        "tested_target_peaks": int(len(peaks)),
        "fdr_significant_peaks": int((peaks["FDR"] < .05).sum()),
        "positive_fdr_significant_peaks": int(((peaks["FDR"] < .05) & (peaks["partial_spearman"] > 0)).sum()),
        "negative_fdr_significant_peaks": int(((peaks["FDR"] < .05) & (peaks["partial_spearman"] < 0)).sum()),
    }


def evidence_matrix(audit: pd.DataFrame, atac: dict) -> None:
    groups = ["Human", "Rat", "Cervid", "Acomys", "Planarian"]
    rows = []
    for gene in ["KLF4", "IGFL1"]:
        for group in groups:
            if group == "Human": species = "homo_sapiens"
            elif group == "Rat": species = "rattus_norvegicus"
            elif group == "Cervid": species = "cervus_hanglu_yarkandensis"
            elif group == "Acomys": species = "cavia_porcellus"
            else: species = "planarian_unavailable"
            hit = audit[(audit.gene == gene) & (audit.requested_species == species)]
            strict = bool(hit["strict_ortholog_present"].iloc[0]) if len(hit) else False
            rows.append({
                "gene": gene, "study_group": group, "strict_sequence_ortholog": strict,
                "sequence_evidence_note": (hit["interpretation"].iloc[0] if len(hit) else "no strict sequence entered"),
                "motif_ATAC_evidence": "rat Fib_K14 KLF4 chromVAR/footprint/peak links" if gene == "KLF4" and group == "Rat" else "not tested",
                "regulon_role": "transcription factor" if gene == "KLF4" else "candidate downstream effector, not a TF",
            })
    pd.DataFrame(rows).to_csv(TAB / "strict_sequence_ATAC_SCENIC_evidence_matrix.csv", index=False)


def main() -> None:
    ensure_dirs()
    audits = []
    raw = {}
    for gene in ["KLF4", "IGFL1"]:
        homologies, payload = ensembl_homologies(gene)
        raw[gene] = payload
        records, audit = select_records(gene, homologies)
        audits.append(audit)
        phylogeny(gene, records)
    audit = pd.concat(audits, ignore_index=True)
    audit.to_csv(TAB / "KLF4_IGFL1_strict_orthology_audit.csv", index=False)
    (SEQ / "Ensembl_Compara_homology_raw.json").write_text(json.dumps(raw, ensure_ascii=False), encoding="utf-8")
    atac = atac_evidence()
    evidence_matrix(audit, atac)
    report = {
        "strict_tree_method": "Ensembl Compara protein orthologs; FAMSA MSA; BLOSUM62 NJ; 300 column bootstraps",
        "important_absence": "IGFL1 was not returned for rat or mouse and is therefore excluded from the IGFL1 tree",
        "proxy_rule": "Guinea pig and naked mole-rat are labeled rodent proxies; they are not claimed as Acomys sequences",
        "planarian_rule": "No validated planarian protein ortholog is placed on either strict tree",
        "ATAC": atac,
    }
    (REP / "strict_orthology_ATAC_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"root": str(ROOT), "audit_rows": len(audit), **atac}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
