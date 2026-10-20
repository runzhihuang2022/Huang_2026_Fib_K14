# Figure3. Fibroblast-to-keratinocyte transdifferentiation

Expression relationships, cross-species support, dermal isolation and mouse dual-recombinase lineage tracing.

## Files

- [Figure3.png](Figure3.png): current manuscript figure reference.
- [source_data/](source_data/README.md): processed inputs and checksums.
- [panel_manifest.csv](panel_manifest.csv): current panel descriptions and evidence/reproduction status.
- [Figure3.html](Figure3.html): local HTML overview with the figure, panel status and source links.

## Run

```text
python Figure3/run_Figure3.py
```

Executed supporting analyses and/or archived sources; see panel-specific status. The available command rebuilds the supporting analyses listed in `analysis_config.json`; it does not reconstruct an unresolved experimental panel from a raster. Output goes to a new directory under `outputs/`.

The final Figure 3M quantitative source/artwork remains pending in the baseline. Genetic lineage tracing, transcriptional trajectories and marker colocalization are recorded as distinct evidence types.

The archived full analysis code and external data requirements are indexed in [the analysis registry](../provenance/analysis_registry.csv). See [the executed reproduction report](../validation/reproduction_report.md) for the exact verified scope.
