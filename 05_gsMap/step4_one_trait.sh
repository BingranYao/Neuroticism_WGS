#!/usr/bin/env bash
set -euo pipefail

TRAIT="${1:?Usage: bash one_trait.sh <trait_id> [raw_tsv] }"
RAW_TSV="${2:-gsmap_${TRAIT}.tsv}"

PROJECT_DIR="/Ashoka/Goni"
GSMAP_WORKDIR="${PROJECT_DIR}/ldscore"
RES="/Ashoka/gsmap/gsMap_resource"
WFILE="${RES}/LDSC_resource/weights_hm3_no_hla/weights."
PREP_PY="${PROJECT_DIR}/scripts/01_prepare_gwas_for_gsmap.py"

SAMPLES=("WT12M01" "WT17M01" "WT2M02")
ANNOT="region"
NUM_PROCESSES="${NUM_PROCESSES:-24}"

mkdir -p \
  "${PROJECT_DIR}/logs" \
  "${PROJECT_DIR}/gwas_prepared" \
  "${PROJECT_DIR}/gwas_sumstats" \
  "${PROJECT_DIR}/cauchy_combined"

# source "$(mamba info --base)/etc/profile.d/conda.sh"
# mamba activate gsMap
cd "${PROJECT_DIR}"

echo "[INFO] TRAIT=${TRAIT}"
echo "[INFO] RAW_TSV=${RAW_TSV}"
echo "[INFO] PROJECT_DIR=${PROJECT_DIR}"
echo "[INFO] GSMAP_WORKDIR=${GSMAP_WORKDIR}"
echo "[INFO] SAMPLES=${SAMPLES[*]}"

[[ -f "${RAW_TSV}" ]] || { echo "[ERROR] missing ${RAW_TSV}" >&2; exit 1; }
[[ -f "${PREP_PY}" ]] || { echo "[ERROR] missing ${PREP_PY}" >&2; exit 1; }

python3 - <<'PY'
import os
import scanpy as sc

project_dir = "/public/home/chenzhh/Ashoka/Goni"
workdir = f"{project_dir}/ldscore"
samples = ["WT12M01", "WT17M01", "WT2M02"]

for s in samples:
    h5 = f"{workdir}/{s}/ST/{s}.gsmap.h5ad"
    if not os.path.exists(h5):
        raise FileNotFoundError(h5)
    ad = sc.read_h5ad(h5)
    if "region" not in ad.obs.columns:
        raise ValueError(f"{h5} does not contain obs['region']; run notebook part ① first")
print("[OK] all selected gsmap h5ad files already contain obs['region']")
PY

PREP_TSV="${PROJECT_DIR}/gwas_prepared/${TRAIT}.for_gsmap.tsv"
SUMSTATS_PREFIX="${PROJECT_DIR}/gwas_sumstats/${TRAIT}"
SUMSTATS_GZ="${SUMSTATS_PREFIX}.sumstats.gz"

echo "[STEP] prepare GWAS"
python3 "${PREP_PY}" \
  --in_tsv "${RAW_TSV}" \
  --out_tsv "${PREP_TSV}"

echo "[STEP] format_sumstats"
gsmap format_sumstats \
  --sumstats "${PREP_TSV}" \
  --out "${SUMSTATS_PREFIX}" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --beta BETA \
  --se SE \
  --p P \
  --n N \
  --maf_min 0

echo "[STEP] run_spatial_ldsc on 3 samples"
for sample in "${SAMPLES[@]}"; do
  echo "  -> ${sample}"
  gsmap run_spatial_ldsc \
    --workdir "${GSMAP_WORKDIR}" \
    --sample_name "${sample}" \
    --sumstats_file "${SUMSTATS_GZ}" \
    --w_file "${WFILE}" \
    --trait_name "${TRAIT}" \
    --num_processes "${NUM_PROCESSES}"
done

echo "[STEP] try official multi-sample cauchy"
if gsmap run_cauchy_combination \
  --workdir "${GSMAP_WORKDIR}" \
  --trait_name "${TRAIT}" \
  --annotation "${ANNOT}" \
  --sample_name_list "${SAMPLES[@]}" \
  --output_file "${PROJECT_DIR}/cauchy_combined/${TRAIT}.combined.cauchy.csv.gz"
then
  echo "[DONE] official multi-sample cauchy finished"
else
  echo "[WARN] multi-sample cauchy failed; fallback to single-sample gsMap cauchy + python combine"

  for sample in "${SAMPLES[@]}"; do
    echo "  -> single-sample cauchy: ${sample}"
    gsmap run_cauchy_combination \
      --workdir "${GSMAP_WORKDIR}" \
      --sample_name "${sample}" \
      --trait_name "${TRAIT}" \
      --annotation "${ANNOT}"
  done

  python3 - "${TRAIT}" "${PROJECT_DIR}" "${GSMAP_WORKDIR}" "${SAMPLES[@]}" <<'PY'
import sys
import os
import numpy as np
import pandas as pd

trait = sys.argv[1]
project_dir = sys.argv[2]
workdir = sys.argv[3]
samples = sys.argv[4:]

def acat(pvals):
    p = pd.to_numeric(pd.Series(pvals), errors="coerce").dropna().astype(float).values
    p = p[(p > 0) & (p <= 1)]
    if len(p) == 0:
        return np.nan
    p = np.clip(p, 1e-300, 1 - 1e-16)
    t = np.tan((0.5 - p) * np.pi)
    T = np.mean(t)
    if T > 1e15:
        out = 1.0 / (np.pi * T)
    elif T < -1e15:
        out = 1.0
    else:
        out = 0.5 - np.arctan(T) / np.pi
    return float(np.clip(out, 1e-300, 1.0))

dfs = []
for s in samples:
    f = os.path.join(workdir, s, "cauchy_combina
    if not os.path.ex
        raise FileNotFoundError(f)
    df = pd.read_csv(f, compression="infer")
    
    df["p_cauchy"] 
    df["sample"] 


x = pd.concat(dfs, axis=0, ignore_index=True)

out = (
    x.groupby("annot
     .agg(
         p_cauchy=lambda v: acat(v),
         p_median=lambda v: float(

     )
)

out_file = os.path.join(project_dir,
out.to_csv(out_file, index=False, compression="gzip")
print("[DONE]
print(out.head())
PY
fi

echo "[ALL DONE] ${TRAI
