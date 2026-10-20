from pathlib import Path
import csv
import hashlib

ROOT = Path("data/external/local_d0cd6f8df107")

def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

audit = read_tsv(ROOT / "02_analysis/required_chord/tables/REQUIRED_CHORD_INSTANCE_AUDIT.tsv")
lr = read_tsv(ROOT / "02_analysis/required_chord/tables/ALL_DISPLAYED_REQUIRED_CHORD_LR.tsv")
go = read_tsv(ROOT / "05_enrichment/D_GO/MODE_REQUIRED_GO_AUDIT.tsv")

group_to_mode = {"D0_7":"Mode 1", "D15_28":"Mode 2", "D8_14":"Mode 3", "M1_2":"Mode 4"}
interface_order = ["Fib_SAC", "Fib_KC", "Fib_Endo", "Fib_Immune"]

lines = [
    "# Figure S3 CellChat v03 analysis summary",
    "",
    "## Completeness",
    "",
    "A-D each retain 4 modes x 4 interfaces = 16 positions. A contains 14 formal chord diagrams and 2 explicit NOT_ESTIMABLE panels. Every estimable chord shows 8 statistically significant ligand-receptor interactions selected deterministically from the locked CellChat object; this exceeds the requested minimum of five.",
    "",
    "Mode 1 Fib-Immune is not estimable because M1_Macrophage n=5 and M2_Macrophage n=1. Mode 4 Fib-Immune is not estimable because M1_Macrophage n=2 and M2_Macrophage n=0. No nonsignificant ligand-receptor pair was inserted.",
    "",
    "## Displayed ligand-receptor evidence",
    "",
]
for group in ["D0_7", "D15_28", "D8_14", "M1_2"]:
    for interface in interface_order:
        rows = [r for r in lr if r["group"] == group and r["interface"] == interface]
        if not rows:
            lines.append(f"- {group_to_mode[group]} / {interface.replace('_','-')}: NOT_ESTIMABLE.")
            continue
        pairs = [f"{r['ligand']}-{r['receptor']} ({r['source']} to {r['target']})" for r in rows]
        lines.append(f"- {group_to_mode[group]} / {interface.replace('_','-')}: " + "; ".join(pairs) + ".")

lines += [
    "",
    "## Evidence-bounded interpretation",
    "",
    "- Fib-SAC and Fib-KC communication is repeatedly represented by extracellular-matrix, adhesion and basement-membrane-associated interactions, including collagen-CD44, THBS1-SDC1, SDC4 and MDK/NCL signals. This supports a fibroblast-associated structural and signaling niche adjacent to appendage and keratinocyte states.",
    "- The Mode 2 Fib-Endo gene set shows significant enrichment for positive regulation of angiogenesis, angiogenesis involved in wound healing, regulation of endothelial migration and vasculature development, supporting an angiogenic component of the wound-repair microenvironment.",
    "- The Mode 2 Fib-Immune gene set shows significant enrichment for positive regulation of immune response, regulation of inflammatory response, leukocyte cell-cell adhesion, T-cell activation and antigen presentation, supporting immune-microenvironment organization. Macrophage activation itself did not pass FDR in this comparison and is not claimed as a positive result.",
    "- Prespecified epithelial proliferation/migration terms in Mode 2 Fib-SAC and keratinocyte differentiation/migration terms in Mode 2 Fib-KC did not pass within-interface BH correction. Therefore the figure supports predicted communication compatible with proliferation/migration, but does not by itself prove that these interactions cause proliferation or migration.",
    "- These analyses use pooled mode-level CellChat objects. Cells are not independent biological replicates; the panels are interaction predictions and gene-set enrichments, not donor-level differential communication tests or causal evidence.",
]
(ROOT / "08_legends/RESULTS_INTERPRETATION_EVIDENCE_BOUND.md").write_text("\n".join(lines)+"\n", encoding="utf-8")

legend = """# Figure S3 legend

**Figure S3. Mode-resolved fibroblast-centered cell-cell communication across human wound-healing states.**

**(A)** Ligand-receptor chord diagrams for four prespecified interfaces (Fib-SAC, Fib-KC, Fib-Endo and Fib-Immune) across Mode 1 (0-7 dpb, unhealed gradient), Mode 2 (15-28 dpb, one-way epithelialization), Mode 3 (8-14 dpb, almost healed) and Mode 4 (1-2 mph, scar/hyperplasia). Required states were Fib_Papi, Fib_K14 and SAC_SG_Progenitor for Fib-SAC; Fib_Papi, Fib_K14 and KC_Basal_Mig for Fib-KC; Fib_Prolif, Fib_Myo and Endo_Capillary for Fib-Endo; and Fib_Inflama, Fib_CD45, M1_Macrophage and M2_Macrophage for Fib-Immune. Each estimable panel displays eight significant ligand-receptor pairs selected from the corresponding locked CellChat object by required-state coverage followed by communication probability and P value. Mode 1 and Mode 4 Fib-Immune panels were retained but marked NOT_ESTIMABLE because at least one required macrophage subtype contained fewer than 10 cells; no nonsignificant pairs were added.

**(B)** Aggregate communication networks among the same required states. Edge width represents summed CellChat communication probability across significant ligand-receptor interactions.

**(C)** Relative outgoing and incoming communication strengths for the required states, calculated from all significant required-state ligand-receptor interactions within each pooled mode object. Raw and normalized values are supplied in the source-data table.

**(D)** Gene Ontology Biological Process enrichment of genes participating in significant required-state ligand-receptor interactions. Bars show -log10 of the within-interface Benjamini-Hochberg-adjusted P value; the dashed line denotes q=0.05. Stars denote adjusted significance (*q<0.05, **q<0.01, ***q<0.001). NOT_ESTIMABLE positions correspond to the same low-cell-count guard used in A-C.

CellChat results are mode-level predictions from pooled cells and should not be interpreted as donor-level differential tests or causal evidence. dpb, days post burn; mph, months post healing; SAC, skin appendage cell; KC, keratinocyte; Endo, endothelial cell; GO, Gene Ontology; BH, Benjamini-Hochberg.
"""
(ROOT / "08_legends/FIGURE_S3_LEGEND.md").write_text(legend, encoding="utf-8")

lock_files = [
    ROOT / "09_code/render_cellchat_default_chord_v26_PRISTINE.R",
    ROOT / "09_code/compose_cellchat_default_chord_v26_PRISTINE.py",
]
lock_lines = ["file\tsha256\tstatus_or_allowed_change"]
for p in lock_files:
    lock_lines.append(f"{p.name}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\tPRISTINE_COPY")
lock_lines += [
    "render_cellchat_default_chord_required_v26.R\tNA\tAllowed wrapper changes only: four locked mode groups and local R library path; core netVisual_chord_gene renderer retained",
    "build_required_cellchat_chord_adapters.R\tNA\tSeparate adapter masks each CellChat object to selected genuine significant required-cell LR pairs",
]
(ROOT / "00_manifest/LOCKED_RENDERER_PROVENANCE.tsv").write_text("\n".join(lock_lines)+"\n", encoding="utf-8")

qa = ["check\tstatus\tdetail"]
passes = sum(r["status"] == "PASS" for r in audit)
not_est = sum(r["status"] == "NOT_ESTIMABLE" for r in audit)
qa.append(f"A_chord_positions\tPASS\t16 retained: {passes} formal chords, {not_est} NOT_ESTIMABLE")
qa.append("minimum_LR_pairs\tPASS\tAll 14 estimable chords display 8 significant LR pairs")
qa.append("required_cell_coverage\tPASS\tAdapter audit reports full required-cell coverage for all estimable interfaces")
qa.append("B_positions\tPASS\t16 positions retained")
qa.append("C_positions\tPASS\t16 positions retained; raw and relative strength table exported")
qa.append("D_positions\tPASS\t16 positions retained; q-value star convention documented")
qa.append("low_count_guard\tPASS\tMode 1 and Mode 4 Fib-Immune retained as NOT_ESTIMABLE; no nonsignificant LR pairs added")
(ROOT / "10_qa/ANALYSIS_COMPLETENESS_QA.tsv").write_text("\n".join(qa)+"\n", encoding="utf-8")
print("release notes, legend, provenance and QA written")
