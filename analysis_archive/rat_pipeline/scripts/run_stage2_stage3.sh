#!/usr/bin/env bash
set -euo pipefail
ROOT=external_asset_6f6697f200d0
while [[ ! -f "$ROOT/working/stage1_rna/STAGE1_DONE" ]]; do
  if ! pgrep -f 'stage1_rna_qc_embedding.py' >/dev/null; then echo 'STAGE1 stopped without completion'; exit 1; fi
  sleep 30
done
if [[ ! -f "$ROOT/working/stage1_rna/STAGE1B_DONE" ]]; then
  env PYTHONPATH=/root/miniconda3/envs/codex-scrna/lib/python3.10/site-packages /root/miniconda3/envs/csc-portable-v1/bin/python "$ROOT/runtime/stage1b_harmony_sensitivity.py" > "$ROOT/logs/stage1b_harmony.log" 2>&1
fi
if [[ ! -f "$ROOT/working/stage2_joint/STAGE2_DONE" ]]; then
  /usr/bin/Rscript "$ROOT/runtime/stage2_atac_joint_integration.R" > "$ROOT/logs/stage2_joint.log" 2>&1
  env PYTHONPATH=/root/miniconda3/envs/codex-scrna/lib/python3.10/site-packages /root/miniconda3/envs/csc-portable-v1/bin/python "$ROOT/runtime/stage2b_joint_cluster.py" > "$ROOT/logs/stage2b_joint_cluster.log" 2>&1
fi
if [[ ! -f "$ROOT/working/stage3_scenic/fibroblast_mouse_ortholog_counts.loom" ]]; then
  /root/miniconda3/envs/csc-pyscenic-0112/bin/python "$ROOT/runtime/stage3_prepare_scenic.py" > "$ROOT/logs/stage3_prepare.log" 2>&1
fi
SC=/root/miniconda3/envs/csc-pyscenic-0112/bin/pyscenic
SOUT="$ROOT/working/stage3_scenic"
if [[ ! -f "$SOUT/adjacencies.tsv" ]]; then
  $SC grn "$SOUT/fibroblast_mouse_ortholog_counts_training_filtered.loom" external_asset_7e9373e27468 \
    -o "$SOUT/adjacencies.tsv" --num_workers 24 --seed 20260811 > "$ROOT/logs/stage3_grn.log" 2>&1
fi
if [[ ! -f "$SOUT/regulons.csv" ]]; then
  $SC ctx "$SOUT/adjacencies.tsv" external_asset_a41ba1021230 \
  --annotations_fname external_asset_c272118497ad \
  --expression_mtx_fname "$SOUT/fibroblast_mouse_ortholog_counts.loom" --output "$SOUT/regulons.csv" \
  --mask_dropouts --num_workers 24 > "$ROOT/logs/stage3_ctx.log" 2>&1
fi
if [[ ! -f "$SOUT/scenic_output.loom" ]]; then
  $SC aucell "$SOUT/fibroblast_mouse_ortholog_counts.loom" "$SOUT/regulons.csv" \
    --output "$SOUT/scenic_output.loom" --num_workers 24 > "$ROOT/logs/stage3_aucell.log" 2>&1
fi
if [[ ! -f "$SOUT/STAGE3_DONE" ]]; then
  /root/miniconda3/envs/csc-pyscenic-0112/bin/python "$ROOT/runtime/stage3_postprocess_scenic.py" > "$ROOT/logs/stage3_postprocess.log" 2>&1
fi
if [[ ! -f "$ROOT/outputs/STAGE4_DONE" ]]; then
  /root/miniconda3/envs/csc-portable-v1/bin/python "$ROOT/runtime/stage4_render_figures.py" > "$ROOT/logs/stage4_render.log" 2>&1
fi
