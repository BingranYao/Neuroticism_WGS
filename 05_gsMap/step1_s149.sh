#!/usr/bin/env bash
set -euo pipefail

CODEDIR="/Ashoka/Goni"
INDIR="/Sankara/New"
OUTROOT="${CODEDIR}/ldscore"

ANNOT="allen_acronym"
LAYER="count"
MAX_JOBS=5

RES="/Ashoka/gsmap/gsMap_resource"
HOMO="${RES}/homologs/mouse_human_homologs.txt"
BFILE_ROOT="${RES}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
KEEP_SNP_ROOT="${RES}/LDSC_resource/hapmap3_snps/hm"
GTF="${RES}/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf"
ENHANCER="${RES}/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed"
GENE_WINDOW=50000

command -v gsmap >/dev/null 2>&1 || { echo "[ERROR] gsmap not found in PATH" >&2; exit 1; }
[[ -f "${CODEDIR}/prepare_st_for_gsmap.py" ]] || { echo "[ERROR] missing ${CODEDIR}/prepare_st_for_gsmap.py" >&2; exit 1; }
[[ -d "${INDIR}" ]] || { echo "[ERROR] missing ${INDIR}" >&2; exit 1; }

mkdir -p "${OUTROOT}/logs/s149" "${OUTROOT}/nohup_logs"

mapfile -t H5ADS < <(find "${INDIR}" -maxdepth 1 -type f -name "*_s149_rot.h5ad" | sort)
(( ${#H5ADS[@]} > 0 )) || { echo "[ERROR] no s149 h5ad files found under ${INDIR}" >&2; exit 1; }

sample_from_h5ad() {
  local base
  base=$(basename "$1")
  printf '%s\n' "${base%%_cellbin_transfer*}"
}

job_count() {
  jobs -rp | wc -l | awk '{print $1}'
}

run_one() {
  local in_h5ad="$1"
  local sample sample_root st_h5ad chr

  sample=$(sample_from_h5ad "${in_h5ad}")
  sample_root="${OUTROOT}/${sample}"
  st_h5ad="${sample_root}/ST/${sample}.gsmap.h5ad"

  mkdir -p "${sample_root}/ST"

  echo "[$(date '+%F %T')] [START] ${sample}"
  echo "INPUT=${in_h5ad}"
  echo "SAMPLE_ROOT=${sample_root}"

  python3 "${CODEDIR}/prepare_st_for_gsmap.py" \
    --in_h5ad "${in_h5ad}" \
    --out_h5ad "${st_h5ad}" \
    --annotation "${ANNOT}" \
    --layer "${LAYER}"

  gsmap run_find_latent_representations \
    --workdir "${OUTROOT}" \
    --sample_name "${sample}" \
    --input_hdf5_path "${st_h5ad}" \
    --annotation "${ANNOT}" \
    --data_layer "${LAYER}"

  gsmap run_latent_to_gene \
    --workdir "${OUTROOT}" \
    --sample_name "${sample}" \
    --annotation "${ANNOT}" \
    --latent_representation "latent_GVAE" \
    --num_neighbour 51 \
    --num_neighbour_spatial 201 \
    --homolog_file "${HOMO}"

  for chr in $(seq 1 22); do
    echo "[$(date '+%F %T')] [${sample}] chr${chr}"
    gsmap run_generate_ldscore \
      --workdir "${OUTROOT}" \
      --sample_name "${sample}" \
      --chrom "${chr}" \
      --bfile_root "${BFILE_ROOT}" \
      --keep_snp_root "${KEEP_SNP_ROOT}" \
      --gtf_annotation_file "${GTF}" \
      --gene_window_size "${GENE_WINDOW}" \
      --enhancer_annotation_file "${ENHANCER}" \
      --snp_multiple_enhancer_strategy "max_mkscore" \
      --gene_window_enhancer_priority "gene_window_first"
  done

  echo "[$(date '+%F %T')] [DONE] ${sample}"
}

fail=0
declare -A PID_TO_SAMPLE=()

for h5ad in "${H5ADS[@]}"; do
  while [ "$(job_count)" -ge "${MAX_JOBS}" ]; do
    if ! wait -n; then
      fail=1
    fi
  done

  sample=$(sample_from_h5ad "${h5ad}")
  log_file="${OUTROOT}/logs/s149/${sample}.step1_3.log"

  (
    run_one "${h5ad}"
  ) > "${log_file}" 2>&1 &

  PID_TO_SAMPLE[$!]="${sample}"
  echo "[$(date '+%F %T')] [SUBMIT] ${sample} -> ${log_file}"
done

for pid in "${!PID_TO_SAMPLE[@]}"; do
  if ! wait "${pid}"; then
    echo "[ERROR] ${PID_TO_SAMPLE[$pid]} failed" >&2
    fail=1
  fi
done

(( fail == 0 )) || exit 1
echo "[DONE] all s149 samples finished"
