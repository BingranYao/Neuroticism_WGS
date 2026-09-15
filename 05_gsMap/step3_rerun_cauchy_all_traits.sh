#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/Ashoka/Goni"
WORKDIR="/Ashoka/Goni/ldscore"
PYTHON_BIN="/envs/Goni/bin/python"
GSMAP_BIN="/envs/Goni/bin/gsmap"
ANNOT="region"

SAMPLES=(WT12M01 WT17M01 WT2M02)
TRAITS=(26653-2.0 26663-2.0 26748-2.0 27211-2.0 27513-2.0)

for sample in "${SAMPLES[@]}"; do
  h5="${WORKDIR}/${sample}/find_latent_representations/${sample}_add_latent.h5ad"
  [[ -f "${h5}" ]] || { echo "[ERROR] missing ${h5}" >&2; exit 1; }
  "${PYTHON_BIN}" - "${h5}" "${ANNOT}" <<'PY'
import sys
import scanpy as sc
h5, annot = sys.argv[1], sys.argv[2]
ad = sc.read_h5ad(h5)
if annot not in ad.obs.columns:
    raise SystemExit(f"[ERROR] {h5} missing obs[{annot!r}]")
print(f"[OK] {annot} exists in {h5}")
PY
done

for trait in "${TRAITS[@]}"; do
  echo "[RUN] trait=${trait}"
  "${GSMAP_BIN}" run_cauchy_combination     --workdir "${WORKDIR}"     --trait_name "${trait}"     --annotation "${ANNOT}"     --sample_name_list "${SAMPLES[@]}"     --output_file "${PROJECT_DIR}/cauchy_combined/${trait}.combined.cauchy.csv.gz"
done

echo "[DONE] all cauchy finished"
