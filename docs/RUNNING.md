# Running the analyses

## Python environment

Tested with Python 3.10.20 and the versions in `requirements.txt`. The desktop application's bundled Python was a different, incomplete environment; validation used a separate environment without modifying the user's scientific environments. Arial is used when installed; a fallback may change layout and PDF hashes on another OS. Numerical comparisons are distinct from pixel equality.

```text
python -m pip install -r requirements.txt
python scripts/validate_release.py
python run_reproduce.py --modules core,figure_s7,mechanism,derived --output-dir outputs/rerun
```

Module `core` runs Table S1 export, H&E summaries, PyDESeq2 fits, atlas checks, spatial-protein rendering, and KC/SAC summaries. `mechanism` recomputes patient aggregation/correlations and redraws six mechanism plots from the complete stored result tables. `derived` checks per-cell CellOracle score aggregation against its retained subtype summaries; it does not refit the CellOracle model.

## Optional R count models

Use R 4.4.1 with DESeq2 1.46.0, edgeR 4.4.2 and limma 3.62.2, as recorded in `environment/historical/mechanism_R_sessionInfo.txt`. R package versions are recorded evidence, not an automatically installed lockfile. Required package installation is explicit in `environment/install_R_packages.R`.

```text
python run_reproduce.py --modules mechanism_r --rscript Rscript --output-dir outputs/rerun_R
```

This refits both experiments from integer counts: gene total >=10, design ~condition, treated/control contrast, DESeq2 independentFiltering=FALSE, and edgeR robust QL sensitivity. It is complementary to the primary PyDESeq2 figures, whose filter is >=10 counts in >=2 libraries with independent filtering enabled.

## Optional full-library GSEA

The exact frozen GO+KEGG GMT is indexed as an external resource because redistribution terms were not established in the supplied material. Supply an authorized copy matching the indexed SHA-256. fgsea 1.32.4 is the recorded analysis version.

```text
python run_reproduce.py --modules gsea --rscript Rscript --external-gmt data/external/GO_KEGG_version_locked.gmt --output-dir outputs/rerun_GSEA
```

The GSEA module starts from the bundled complementary DESeq2 statistics. Use newly generated tables explicitly if replacing these inputs after a new R fit; do not mix filters. Seed 20261017; minSize=10, maxSize=500, eps=1e-6, nproc=1. GO+KEGG sets are adjusted together before selecting display pathways.

## Upstream atlas and historical code

Additional Python/R/IPython/Groovy scripts are under `analysis_archive/`. They preserve scientific procedures but include several historical stages. Workstation paths have been replaced by external-asset aliases; archival scripts require input binding, dependency matching and version adjudication before execution. Some original notebooks contain historical exploratory cells; they are not run by `run_reproduce.py`.

Examples of retained entry points include human `formal_prepare.py`, `formal_inference.R`, `run_cellchat.R`, `run_monocle3.R`, and `run_spatial_projection.py`, the rat RNA/ATAC integration stages, CellOracle preparation/simulation scripts, and original spatial neighborhood analysis. Consult each script's command-line arguments. Large author objects, trained models, reference genomes and motif libraries are not automatically fetched.

Raw-read alignment, full atlas integration/clustering, raw-fragment ATAC, new microscopy segmentation and biological experiments are outside this consolidation run. Availability of code is reported separately from execution.
