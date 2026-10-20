# Reproduction report — 20 October 2026

Indexed 24 scientific analysis scopes and 115 current panels. Executed 12 commands in 6 groups; every command exited 0.

Compared 44 numerical tables: 44 matched after stable-key alignment; 0 require interpretation. Numeric tolerance is rtol=1e-7/atol=1e-9; probability/FDR fields use rtol=1e-6/atol=1e-300. Public Table S1 is compared to its public derivative, with original identity translations kept privately. See numerical_comparison.csv.

## Actually executed

- core / audit_table_s1.py (0.695 s)
- core / rebuild_current_histology.py (2.204 s)
- core / reanalyze_bulk.py (48.039 s)
- core / rebuild_atlas_audits.py (6.792 s)
- core / rebuild.py (11.881 s)
- core / rebuild_current_supplements.py (6.96 s)
- figure_s7 / build_figure_s7.py (10.094 s)
- mechanism / analyze_human.py (21.942 s)
- mechanism / make_figures.py (3.355 s)
- mechanism_r / run_bulk.R (19.579 s)
- gsea / run_gsea.R (10.277 s)
- derived / rebuild_derived.py (1.367 s)

## Limits

This is a consolidated processed-data and count-model release, not a complete raw-read-to-paper rerun. Raw alignment, full atlas integration/clustering, fragment-level ATAC, fresh CellOracle model fitting, and new imaging/biological experiments were not performed. Archive scripts remain not_run. Numeric agreement does not resolve the underlying experimental metadata gaps.

R count-model and enrichment runs used the frozen external GO+KEGG library locally. Its hash is indexed; the resource is not redistributed in this repository. Primary PyDESeq2 and complementary R model filters remain distinct.

Exact Figure 3M and Figure 6E/F/L final quantitative inputs remain missing. The public Table S1 is a pseudonymized derivative. Original source paths/linkage maps and failed development attempts remain outside the release.
