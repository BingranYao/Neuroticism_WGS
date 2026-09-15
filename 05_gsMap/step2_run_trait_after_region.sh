#!/usr/bin/env bash
set -euo pipefail

TRAIT="${1:?Usage: bash run_trait_after_region.sh <trait_id> [raw_tsv]}"
RAW_TSV="${2:-/gsmap_${TRAIT}.tsv}"

PROJECT_DIR="/Ashoka/Goni"
GSMAP_WORKDIR="/Ashoka/Goni/ldscore"
PYTHON_BIN="/envs/Goni/bin/python"
GSMAP_BIN="/envs/Goni/bin/gsmap"
WFILE="/Ashoka/gsmap/gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights."
ANNOT="region"
NUM_PROCESSES="${NUM_PROCESSES:-200}"
FORCE_SPATIAL_LDSC="${FORCE_SPATIAL_LDSC:-0}"

SAMPLES=("WT12M01" "WT17M01" "WT2M02")

mkdir -p "${PROJECT_DIR}/logs" "${PROJECT_DIR}/gwas_prepared" "${PROJECT_DIR}/gwas_sumstats" "${PROJECT_DIR}/cauchy_combined"

[[ -x "${PYTHON_BIN}" ]] || { echo "[ERROR] missing PYTHON_BIN=${PYTHON_BIN}" >&2; exit 1; }
[[ -x "${GSMAP_BIN}" ]] || { echo "[ERROR] missing GSMAP_BIN=${GSMAP_BIN}" >&2; exit 1; }
[[ -f "${RAW_TSV}" ]] || { echo "[ERROR] missing RAW_TSV=${RAW_TSV}" >&2; exit 1; }

for sample in "${SAMPLES[@]}"; do
  h5="${GSMAP_WORKDIR}/${sample}/ST/${sample}.gsmap.h5ad"
  [[ -f "${h5}" ]] || { echo "[ERROR] missing ${h5}" >&2; exit 1; }
  "${PYTHON_BIN}" - "${h5}" <<'PY'
import sys, scanpy as sc
h5 = sys.argv[1]
ad = sc.read_h5ad(h5)
if "region" not in ad.obs.columns:
    raise SystemExit(f"[ERROR] {h5} missing obs['region']; run notebook Cell 2-3 first")
print(f"[OK] region exists in {h5}")
PY
done

PREP_TSV="${PROJECT_DIR}/gwas_prepared/${TRAIT}.for_gsmap.tsv"
SUMSTATS_PREFIX="${PROJECT_DIR}/gwas_sumstats/${TRAIT}"
SUMSTATS_GZ="${SUMSTATS_PREFIX}.sumstats.gz"
PREP_PY="${PROJECT_DIR}/scripts/01_prepare_gwas_for_gsmap.py"

echo "[STEP] prepare GWAS"
"${PYTHON_BIN}" "${PREP_PY}" --in_tsv "${RAW_TSV}" --out_tsv "${PREP_TSV}"

echo "[STEP] format_sumstats"
"${GSMAP_BIN}" format_sumstats   --sumstats "${PREP_TSV}"   --out "${SUMSTATS_PREFIX}"   --snp SNP   --a1 A1   --a2 A2   --beta BETA   --se SE   --p P   --n N   --maf_min 0

echo "[STEP] run_spatial_ldsc"
for sample in "${SAMPLES[@]}"; do
  out_csv="${GSMAP_WORKDIR}/${sample}/spatial_ldsc/${sample}_${TRAIT}.csv.gz"
  if [[ "${FORCE_SPATIAL_LDSC}" != "1" && -f "${out_csv}" ]]; then
    echo "  -> skip existing ${out_csv}"
    continue
  fi
  echo "  -> run ${sample}"
  "${GSMAP_BIN}" run_spatial_ldsc     --workdir "${GSMAP_WORKDIR}"     --sample_name "${sample}"     --sumstats_file "${SUMSTATS_GZ}"     --w_file "${WFILE}"     --trait_name "${TRAIT}"     --num_processes "${NUM_PROCESSES}"
done

echo "[STEP] official multi-sample cauchy"
"${GSMAP_BIN}" run_cauchy_combination   --workdir "${GSMAP_WORKDIR}"   --trait_name "${TRAIT}"   --annotation "region"   --sample_name_list "${SAMPLES[@]}"   --output_file "${PROJECT_DIR}/cauchy_combined/${TRAIT}.combined.cauchy.csv.gz"

echo "[DONE] TRAIT=${TRAIT}"
