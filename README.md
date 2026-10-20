# Huang_2026_Fib_K14

Code and processed data for the human wound epithelialization and Fib_K14 study.

**Manuscript:** A spatiotemporal human wound atlas reveals injury-induced fibroblast-to-keratinocyte transdifferentiation during epithelialization  
**Authors:** Huang et al. — full author metadata in [CITATION.cff](CITATION.cff)  
**Repository:** https://github.com/runzhihuang2022/Huang_2026_Fib_K14  
**Scientific baseline:** Manuscript V12, 19 October 2026. Repository organized on 20 October 2026. Publication DOI and final journal details will be added when assigned.

The study combines a human spatiotemporal wound atlas, protein imaging, cross-species profiling, lineage tracing and fibroblast perturbations to investigate a stromal contribution to epithelialization. This repository provides the analysis code, plot-ready source tables and figure references.

## Repository organization

Following a figure-based organization, each `Figure` directory contains the figure reference, panel descriptions, source-data index, available execution entry point, and an HTML overview. Shared code and canonical inputs remain in their original Python/R modules to preserve one traceable analysis source.

| Folder | Content | Source data and report |
|---|---|---|
| [Figure1](Figure1/) | Human wound epithelialization atlas | [Data](Figure1/source_data/README.md) · [Report](Figure1/Figure1.html) |
| [Figure2](Figure2/) | Identification of Fib_K14 | [Data](Figure2/source_data/README.md) · [Report](Figure2/Figure2.html) |
| [Figure3](Figure3/) | Fibroblast-to-keratinocyte transdifferentiation | [Data](Figure3/source_data/README.md) · [Report](Figure3/Figure3.html) |
| [Figure4](Figure4/) | KLF4-associated epithelial regulation | [Data](Figure4/source_data/README.md) · [Report](Figure4/Figure4.html) |
| [Figure5](Figure5/) | The epithelialization-front niche | [Data](Figure5/source_data/README.md) · [Report](Figure5/Figure5.html) |
| [Figure6](Figure6/) | IGFL1-responsive epithelial features | [Data](Figure6/source_data/README.md) · [Report](Figure6/Figure6.html) |
| [FigureS1](FigureS1/) | Atlas sampling and quality control | [Data](FigureS1/source_data/README.md) · [Report](FigureS1/FigureS1.html) |
| [FigureS2](FigureS2/) | Keratinocyte state analysis | [Data](FigureS2/source_data/README.md) · [Report](FigureS2/FigureS2.html) |
| [FigureS3](FigureS3/) | Skin appendage cell analysis | [Data](FigureS3/source_data/README.md) · [Report](FigureS3/FigureS3.html) |
| [FigureS4](FigureS4/) | Fibroblast landscape and supporting comparisons | [Data](FigureS4/source_data/README.md) · [Report](FigureS4/FigureS4.html) |
| [FigureS5](FigureS5/) | Supporting spatial analyses | [Data](FigureS5/source_data/README.md) · [Report](FigureS5/FigureS5.html) |
| [FigureS6](FigureS6/) | Adult-rat multi-omics validation | [Data](FigureS6/source_data/README.md) · [Report](FigureS6/FigureS6.html) |
| [FigureS7](FigureS7/) | KLF4 and IGFL1 response comparison | [Data](FigureS7/source_data/README.md) · [Report](FigureS7/FigureS7.html) |
| [TableS1](TableS1/) | Complete public workbook and CSV exports | [Workbook](TableS1/Table_S1.xlsx) |

- [`data/`](data/processed/README.md): processed-data entry points, data dictionary and external large-data index.
- [`script/`](script/README.md): common analysis and figure-run entry points.
- [`modules/`](modules/): executable analyses and canonical inputs.
- [`analysis_archive/`](analysis_archive/README.md): original atlas, spatial, trajectory, CellOracle, chromatin, communication, cross-species and imaging code.
- [`environment/`](environment/): pinned Python dependencies, resolved environment and recorded R sessions.
- [`validation/`](validation/reproduction_report.md): execution logs and numerical comparisons.
- [`provenance/`](provenance/analysis_registry.csv): 24 analysis scopes, 115 panel mappings, checksums and source transformations.

## Requirements

The tested Python environment is **Python 3.10.20**; dependencies are pinned in [requirements.txt](requirements.txt), with a resolved Windows environment in [environment/requirements-lock-windows-py310.txt](environment/requirements-lock-windows-py310.txt). Complementary R analyses used **R 4.4.1**, DESeq2 1.46.0, edgeR 4.4.2, limma 3.62.2 and fgsea 1.32.4. See [detailed execution instructions](docs/RUNNING.md).

## How to run

```bash
git clone https://github.com/runzhihuang2022/Huang_2026_Fib_K14.git
cd Huang_2026_Fib_K14
python -m venv .venv
# Activate the environment using the command appropriate for your system.
python -m pip install -r requirements.txt
python scripts/validate_release.py
python run_reproduce.py
```

For a figure-specific entry, for example:

```bash
python FigureS7/run_FigureS7.py
```

The default reproduction rebuilds the registered processed-data/count-model tasks in a separate output directory. The local verification executed **12 tasks across six groups**, with **44 result tables matching after stable-key alignment**. All six S7 panels also matched their retained render at the recorded comparison resolution. See [the reproduction report](validation/reproduction_report.md) for the exact scope.

The archive also contains historical stages that were not re-executed. Missing raw/large inputs and unresolved experimental panel sources remain explicit in the [panel coverage table](provenance/figure_panel_coverage.csv) and [external data index](data/external_assets.csv). Figure-specific commands rebuild registered supporting analyses; they do not manufacture missing experimental measurements.

## Data and code availability

Analysis code and the processed data bundled here are available at this repository. Large H5AD/RDS matrices, raw reads/fragments, full-slide images and some external reference libraries are indexed separately with their available hashes and access status. The repository does not substitute for a raw-sequencing archive. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

## Citation

Please cite the associated manuscript:

Huang et al. *A spatiotemporal human wound atlas reveals injury-induced fibroblast-to-keratinocyte transdifferentiation during epithelialization.* Manuscript; publication details pending.

Machine-readable author and repository metadata are provided in [CITATION.cff](CITATION.cff). Rights information is recorded in [LICENSE.proposed.txt](LICENSE.proposed.txt); a final open-source code license has not been assigned. Third-party resources retain their original terms.
