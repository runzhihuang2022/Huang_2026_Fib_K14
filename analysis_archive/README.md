# Scientific analysis archive

This directory preserves author-project scientific code and processed result tables. It excludes installed packages, runtime binaries, Office documents, reviewer correspondence and raw large objects. Original source hashes and release transformations are recorded in `../provenance/source_transformations.csv`; full workstation paths are retained privately outside the release.

The archive spans earlier figure numbering and analysis versions. Historical scripts are **not_run in this consolidation**, except for the explicitly named derived checks. Their results are not automatically reassigned to present manuscript panels. Portable aliases require configuration and external assets; read-only source preservation does not turn an archival script into a complete executable pipeline.

The `original_atlas` R script and RCTD notebook provide the older broad workflow. `portable_human` contains parameterized human pseudobulk, enrichment, CellChat, trajectory and spatial-projection code. `human_celloracle` retains model inputs, scored vectors, GRN tables and simulation code. `rat_pipeline` and `rat_atac` retain processed-stage RNA/ATAC analyses. `historical_figure_units` retains spatial and communication analyses. `imaging` contains image-registration, segmentation and quantification code; raw registered microscopy data are external.

The legacy 20260402 Stereo-seq export contains an incomplete dictionary block in its source and is retained as .py.fragment.txt. It is not executable code; the parameterized spatial pipeline is indexed separately.
