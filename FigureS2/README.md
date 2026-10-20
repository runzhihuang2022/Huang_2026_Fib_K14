# FigureS2. Keratinocyte state analysis

Keratinocyte states, temporal composition, marker resources and selected functional programs.

## Files

- [FigureS2.png](FigureS2.png): current manuscript figure reference.
- [source_data/](source_data/README.md): processed inputs and checksums.
- [panel_manifest.csv](panel_manifest.csv): current panel descriptions and evidence/reproduction status.
- [FigureS2.html](FigureS2.html): local HTML overview with the figure, panel status and source links.

## Run

```text
python FigureS2/run_FigureS2.py
```

Executed supporting analyses and/or archived sources; see panel-specific status. The available command rebuilds the supporting analyses listed in `analysis_config.json`; it does not reconstruct an unresolved experimental panel from a raster. Output goes to a new directory under `outputs/`.

Temporal sampling groups and morphological Mode 1–4 are separate variables; historical N labels and the method/figure version must be reconciled as recorded in the technical audit.

The archived full analysis code and external data requirements are indexed in [the analysis registry](../provenance/analysis_registry.csv). See [the executed reproduction report](../validation/reproduction_report.md) for the exact verified scope.
