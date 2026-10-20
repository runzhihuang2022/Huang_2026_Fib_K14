# FigureS1. Atlas sampling and quality control

Human assay sample aliases, retained QC observations, and supporting atlas resources.

## Files

- [FigureS1.png](FigureS1.png): current manuscript figure reference.
- [source_data/](source_data/README.md): processed inputs and checksums.
- [panel_manifest.csv](panel_manifest.csv): current panel descriptions and evidence/reproduction status.
- [FigureS1.html](FigureS1.html): local HTML overview with the figure, panel status and source links.

## Run

```text
python FigureS1/run_FigureS1.py
```

Executed supporting analyses and/or archived sources; see panel-specific status. The available command rebuilds the supporting analyses listed in `analysis_config.json`; it does not reconstruct an unresolved experimental panel from a raster. Output goes to a new directory under `outputs/`.



The archived full analysis code and external data requirements are indexed in [the analysis registry](../provenance/analysis_registry.csv). See [the executed reproduction report](../validation/reproduction_report.md) for the exact verified scope.
