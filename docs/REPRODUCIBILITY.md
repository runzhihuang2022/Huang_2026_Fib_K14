# Human wound epithelialization: analysis and figure reproduction

This repository accompanies **A spatiotemporal human wound atlas reveals injury-induced fibroblast-to-keratinocyte transdifferentiation during epithelialization**. The scientific baseline is Manuscript V12 (19 October 2026); this consolidated code/data package was prepared on 20 October 2026.

## Start here

The executable modules rebuild processed-data figures and tables and refit the two bulk count models. The analysis archive preserves additional atlas, spatial, trajectory, CellOracle, chromatin, communication, cross-species and imaging code. **An archived script or a supplied image is not an executed end-to-end reproduction.** Coverage and unresolved inputs are reported per analysis and panel.

```text
python -m venv .venv
# Activate .venv using the command for your operating system.
python -m pip install -r requirements.txt
python scripts/validate_release.py
python run_reproduce.py
```

The default command runs six core tasks, the V12 Figure S7 build, patient-level mechanism analyses and their six plots, and retained CellOracle/orthology checks. It creates a separate `outputs/rerun/` tree and refuses to replace a previous execution report. Use `--output-dir outputs/rerun_02` for a later run. See [execution instructions](../docs/RUNNING.md) for optional R models and GSEA.

## Contents

| Directory | Contents |
|---|---|
| `modules/core/` | Public Table S1 derivative; atlas QC/retained embeddings; KC/SAC summaries; H&E summaries; primary PyDESeq2 count models; spatial-protein figure reconstruction |
| `modules/figure_s7/` | V12 S7A–F source tables, unchanged experimental H&E fields, and assembly code |
| `modules/mechanism/` | Two count matrices, complete DESeq2/edgeR and enrichment result tables, patient-level inputs, and original R/Python analysis code |
| `analysis_archive/` | Additional scientific analysis code and processed source tables; historical figure numbering retained |
| `paper_figures/` | Current manuscript image references; these images do not replace numerical source data |
| `expected_outputs/` | Retained comparison tables and S7 vector references |
| `data/` | External dataset index and schema documentation |
| `provenance/` | Panel inventory, analysis registry, source transformations, checksums and reproduction manifest |
| `validation/` | Actual rerun scope, comparisons, release checks and logs |
| `environment/` | Pinned Python requirements and recorded R/Python environments |

## Scope and sample units

The registered human data comprise 38 scRNA-seq libraries (279,305 cells), 11 Stereo-seq sections (387,061 bins) and two Visium sections (4,449 spots). These are different assay units, not 51 independent participants. The frozen adult-rat packages comprise 53 scRNA-seq libraries (450,282 cells), six adult spatial sections (324,382 bins), and 19 scATAC libraries (28,618 post-filter cells). Historical 30-section spatial outputs are retained only where explicitly labeled as history; they do not extend the current adult scope. Spatial proteomics contains one section and 2,029 pixels.

Mode 1–4 morphology and sampling-time groups are distinct variables. Several older KC/SAC summaries use temporal categories; their rerun is not proof of a morphology-stratified analysis. Likewise, the legacy nearest-neighbor density contrast must not be relabeled as an executed Milo model.

## Scientific interpretation

The manuscript reports mouse dual-recombinase evidence for directional fibroblast-to-keratinocyte transdifferentiation. Human transcriptional profiles, spatial proximity and CellOracle predictions provide different evidence types. The current regulatory interpretation is partial convergence of IGFL1-responsive and KLF4-associated repair programs. The package does not establish a linear KLF4→IGFL1→IGFLR1 pathway, receptor dependence, or MAPK/PI3K–AKT activation.

The primary PyDESeq2 figure analysis and the complementary R DESeq2/edgeR analysis use different documented filters. Their estimates must not be silently substituted for one another. All-library GSEA adjusted P values are retained; the displayed subset is never used as a new multiple-testing family.

## Data access and identifiers

Large H5AD/RDS matrices, raw reads/fragments and full-slide images are external. [The external index](../data/external_assets.csv) records intended destinations, registered hashes and access status; blank accession fields mean that a verified public accession was not supplied. Public deposition and access permissions remain author responsibilities. No pending accession is presented as downloadable data.

Release tables use module-specific pseudonyms; original linkage keys and workstation locations are excluded. See [identifier policy](../docs/IDENTIFIERS.md). `Table_S1_source.xlsx` in this repository is a public derivative, with its own hash, and is not the unchanged submission workbook.

## Reproduction status and known gaps

Read [the analysis index](../docs/ANALYSIS_INDEX.md), [panel coverage](../provenance/figure_panel_coverage.csv), and [reproduction report](../validation/reproduction_report.md). The exact Figure 3M and Figure 6E/F/L artwork/quantitative source tables remain absent from the current baseline. Full raw ATAC reconstruction also lacks fragment/Arrow inputs. These gaps remain explicit even when another module completes successfully.

## Citation, rights and upload

`CITATION.cff` records the manuscript author list; no journal acceptance, DOI or repository URL is invented. `LICENSE.proposed.txt` contains an unadopted code-license suggestion. Third-party data/resources retain their own terms; no blanket license is assigned to all research data.

Upload this repository directory only. The parent working directory contains private source locations and linkage maps and is not part of the release. See [GitHub upload instructions](../docs/GITHUB_UPLOAD.md).
