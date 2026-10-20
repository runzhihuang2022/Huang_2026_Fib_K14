# Figure6. IGFL1-responsive epithelial features

Local signaling candidates, fibroblast perturbation, bulk expression responses and wound-repair observations.

## Files

- [Figure6.png](Figure6.png): current manuscript figure reference.
- [source_data/](source_data/README.md): processed inputs and checksums.
- [panel_manifest.csv](panel_manifest.csv): current panel descriptions and evidence/reproduction status.
- [Figure6.html](Figure6.html): local HTML overview with the figure, panel status and source links.

## Run

```text
python Figure6/run_Figure6.py
```

Executed supporting analyses and/or archived sources; see panel-specific status. The available command rebuilds the supporting analyses listed in `analysis_config.json`; it does not reconstruct an unresolved experimental panel from a raster. Output goes to a new directory under `outputs/`.

The final Figure 6E/F/L source tables/artwork remain pending. Primary PyDESeq2 figure statistics and complementary R DESeq2/edgeR results retain their separate parameter sets.

The archived full analysis code and external data requirements are indexed in [the analysis registry](../provenance/analysis_registry.csv). See [the executed reproduction report](../validation/reproduction_report.md) for the exact verified scope.
