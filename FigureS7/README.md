# FigureS7. KLF4 and IGFL1 response comparison

Retained H&E fields, shared gene responses, full-library GSEA, patient-level associations and receptor-expression distribution.

## Files

- [FigureS7.png](FigureS7.png): current manuscript figure reference.
- [source_data/](source_data/README.md): processed inputs and checksums.
- [panel_manifest.csv](panel_manifest.csv): current panel descriptions and evidence/reproduction status.
- [FigureS7.html](FigureS7.html): local HTML overview with the figure, panel status and source links.

## Run

```text
python FigureS7/run_FigureS7.py
```

Rebuilt computational panels. The available command rebuilds the supporting analyses listed in `analysis_config.json`; it does not reconstruct an unresolved experimental panel from a raster. Output goes to a new directory under `outputs/`.

S7A/B preserve supplied experimental rasters. S7C–F are data-derived plots. The regulatory interpretation is partial convergence of IGFL1-responsive and KLF4-associated repair programs.

The archived full analysis code and external data requirements are indexed in [the analysis registry](../provenance/analysis_registry.csv). See [the executed reproduction report](../validation/reproduction_report.md) for the exact verified scope.
