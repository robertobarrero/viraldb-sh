#!/bin/bash -l
#SBATCH --job-name=fetch_viralDB
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=fetch_viralDB.%j.out
#SBATCH --error=fetch_viralDB.%j.err

# Author: Roberto A. Barrero, eResearch, Research Infrastructure, Academic Division, QUT
# Date: 2026-03-10
# This version runs in a HPC cluster using PBS Pro scheduler
# See README.md for details on how to run this script

set -euo pipefail

# -----------------------------
# Pipeline installation dir (STATIC)
# -----------------------------
PIPELINE_DIR="${HOME}/git/viraldb-sh"

# -----------------------------
# Determine run directory (PBS / SLURM / local)
# -----------------------------
if [[ -n "${PBS_O_WORKDIR:-}" ]]; then
    RUN_DIR="$PBS_O_WORKDIR"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    RUN_DIR="$SLURM_SUBMIT_DIR"
else
    RUN_DIR="$(pwd)"
fi

cd "$RUN_DIR" || exit 1
LOGDIR="${RUN_DIR}/logs"

# -----------------------------
# Load config + utils
# -----------------------------
CONFIG_FILE="${PIPELINE_DIR}/config/config_viralDB.txt"
source "${CONFIG_FILE}"
source "${PIPELINE_DIR}/lib/pipeline_utils.sh"

parse_modes "$@"

# -----------------------------
# Helper: require directory exists
# -----------------------------
require_dir() {
  local d="$1"
  [[ -d "$d" ]] || die "Required directory missing: $d"
}

# -----------------------------
# Setup log dir (in run folder)
# -----------------------------
mkdir -p "${LOGDIR}"
LOGFILE="${LOGDIR}/viraldb_${DATE}.log"

# If dry-run, don't tee (keep clean). Otherwise log everything.
if [[ "${DRY_RUN}" -eq 0 ]]; then
  exec > >(tee -i "${LOGFILE}")
  exec 2>&1
fi

log "PIPELINE: ${PIPELINE_NAME} v${PIPELINE_VERSION}"
log "OUTDIR   : ${OUTDIR}"
log "DRY_RUN  : ${DRY_RUN}"
log "VALIDATE : ${VALIDATE_ONLY}"

# -----------------------------
# Activate conda environment with necessary tools
# -----------------------------
source ~/.bashrc
conda activate "$CONDAENV"

# -----------------------------
# Validation (fast, fail early)
# -----------------------------
log "Validating environment and inputs..."

require_exec python3
require_exec datasets
require_exec cd-hit-est

require_file "${VIRAL_FAMILIES_TSV}"

require_dir_writable "${OUTDIR}"
require_dir_writable "${LOGDIR}"

# Validate required scripts exist
require_file "${BIN_DIR}/download_viral_sequences_v0.9.py"
require_file "${BIN_DIR}/enrich_headers_with_lineage_from_datasets_jsonl_families_v0.8.py"
require_file "${BIN_DIR}/enrich_headers_for_unclassified_viroids_v0.1.py"
require_file "${BIN_DIR}/run_enrich_taxonomy_per_family.sh"
require_file "${BIN_DIR}/sort_fasta_by_species_v1.5.py"
require_file "${BIN_DIR}/cluster_and_select_reps_v0.2.py"
require_file "${BIN_DIR}/summarise_cdhit_clstr_with_fasta.py"
require_file "${BIN_DIR}/phase2_select_reps_from_cdhit_v0.3.py"

log "Validation OK."

if [[ "${VALIDATE_ONLY}" -eq 1 ]]; then
  log "Validation-only mode requested; exiting."
  exit 0
fi

# -----------------------------
# Manifest
# -----------------------------
if [[ "${DRY_RUN}" -eq 0 ]]; then
  init_manifest
  record_versions
fi

## -----------------------------
## STEP 1: Download
## -----------------------------
#log "🔹🔹🔹 Step 1: download NCBI viral + ViroidDB"
#
#run_cmd "python3 ${BIN_DIR}/download_viral_sequences_v0.9.py \
#  --ncbi_viral \
#  --viroid_db \
#  --viral_families ${VIRAL_FAMILIES_TSV} \
#  --merge_reports \
#  --datasets_retries 3 \
#  --sleep_between_downloads 10 \
#  -o ${OUTDIR} \
#  --verbose"

# -----------------------------
# STEP 2: Enrich NCBI FASTA per family
# -----------------------------
log "🔹🔹🔹 Step 2: enrich NCBI viral headers per family"

require_dir "${OUTDIR}/families"

log "Detected family directories:"
run_cmd "ls -1 ${OUTDIR}/families"

run_cmd "bash ${BIN_DIR}/run_enrich_taxonomy_per_family.sh \
  ${OUTDIR}/families \
  ${BIN_DIR}/enrich_headers_with_lineage_from_datasets_jsonl_families_v0.8.py"

# -----------------------------
# STEP 2b: Combine enriched family FASTA files
# -----------------------------
log "🔹🔹🔹 Step 2b: combine enriched family FASTA files"

run_cmd "find ${OUTDIR}/families -type f -name '*_ncbi_viral_taxonomy.fasta' -print | sort > ${OUTDIR}/family_fasta_list.txt"

run_cmd "if [[ -s ${OUTDIR}/family_fasta_list.txt ]]; then \
  xargs cat < ${OUTDIR}/family_fasta_list.txt > ${OUTDIR}/ncbi_viral__ALL_FAMILIES__taxonomy.fasta; \
else \
  echo '[ERROR] No enriched family FASTA files found'; exit 1; \
fi"

# -----------------------------
# STEP 3: Enrich viroids
# -----------------------------
log "🔹🔹🔹 Step 3: enrich unclassified viroids headers"

run_cmd "python3 ${BIN_DIR}/enrich_headers_for_unclassified_viroids_v0.1.py \
  --fasta ${OUTDIR}/unclassified.fasta \
  --datasets_report ${OUTDIR}/ncbi_viral__ALL_FAMILIES.data_report.jsonl \
  --source ViroidDB \
  --out_fasta ${OUTDIR}/unclassified_viroids_taxonomy.fasta \
  --out_metadata ${OUTDIR}/unclassified_viroids_taxonomy_metadata.tsv"

# -----------------------------
# STEP 4: Merge + filter + sort
# -----------------------------
log "🔹🔹🔹 Step 4: merge + filter + sort"

run_cmd "python3 ${BIN_DIR}/sort_fasta_by_species_v1.5.py \
  --fasta ${OUTDIR}/ncbi_viral__ALL_FAMILIES__taxonomy.fasta \
  --fasta2 ${OUTDIR}/unclassified_viroids_taxonomy.fasta \
  --min_len ${MIN_LEN} \
  --max_Ns_fraction ${MAX_NS_FRACTION} \
  --out_merged_fasta ${OUTDIR}/ncbi_viral_unclassifiedViroids__ALL_FAMILIES__taxonomy.fasta \
  --out_sorted_fasta ${OUTDIR}/ncbi_viral_unclassifiedViroids__ALL_FAMILIES__taxonomy_filtered.fasta"

# -----------------------------
# STEP 5: Cluster filtered sequences
# -----------------------------
log "🔹🔹🔹 Step 5: cluster filtered sequences"

# cd-hit-est -M expects MB; config provides MB
: "${CDHIT_MEM_MB:=0}"   # 0 = unlimited (matches cd-hit-est behavior)

run_cmd "python3 ${BIN_DIR}/cluster_and_select_reps_v0.2.py \
  --fasta ${OUTDIR}/ncbi_viral_unclassifiedViroids__ALL_FAMILIES__taxonomy_filtered.fasta \
  --outdir ${OUTDIR}/representatives \
  --group all \
  --identity 1.0 0.995 0.99 0.98 0.95 \
  --threads ${CPUS} \
  --memory ${CDHIT_MEM_MB} \
  --verbose"

# -----------------------------
# STEP 6: Summarise clusters
# -----------------------------
log "🔹🔹🔹 Step 6: summarise clusters"

IDENTITIES=("1.0" "0.995" "0.99" "0.98" "0.95")

for ID in "${IDENTITIES[@]}"; do

  CDIR=$(printf "c%.6f" "${ID}")

  if [[ "${ID}" == "1.0" ]]; then
    LABEL="c1000"
  elif [[ "${ID}" == "0.995" ]]; then
    LABEL="c995"
  elif [[ "${ID}" == "0.99" ]]; then
    LABEL="c990"
  elif [[ "${ID}" == "0.98" ]]; then
    LABEL="c980"	  
  elif [[ "${ID}" == "0.95" ]]; then
    LABEL="c950"
  else
    LABEL="${CDIR}"
  fi

  CLSTR_FILE="${OUTDIR}/representatives/all_groups/ALL/${CDIR}/reps.fasta.clstr"

  if [[ ! -s "${CLSTR_FILE}" ]]; then
    log "[WARN] Missing .clstr for ${CDIR}: ${CLSTR_FILE} (skipping)"
    continue
  fi

  log "🔹🔹🔹 Summarising ${LABEL} using ${CLSTR_FILE}"

  run_cmd "python3 ${BIN_DIR}/summarise_cdhit_clstr_with_fasta.py \
    --clstr ${CLSTR_FILE} \
    --fasta ${OUTDIR}/ncbi_viral_unclassifiedViroids__ALL_FAMILIES__taxonomy_filtered.fasta \
    --identity_label ${LABEL} \
    --group ALL \
    --out_summary ${OUTDIR}/clusters_summary_${LABEL}.tsv \
    --out_mixed ${OUTDIR}/clusters_mixed_species_members_${LABEL}.tsv \
    --out_flat ${OUTDIR}/clusters_${LABEL}.tsv"
done

# -----------------------------
# STEP 7: Selection of representatives
# -----------------------------
log "🔹🔹🔹 Step 7: Selection of representatives"

run_cmd "python3 ${BIN_DIR}/phase2_select_reps_from_cdhit_v0.3.py \
  --reps_all_dir ${OUTDIR}/representatives \
  --policy ${PHASE2_POLICY} \
  --refseq_min_len_frac ${REFSEQ_MIN_LEN_FRAC} \
  --prefix ${DATE}_viralDB \
  --verbose"

# -----------------------------
# Checksums of key outputs + finish
# -----------------------------
if [[ "${DRY_RUN}" -eq 0 ]]; then
  checksum_files \
    "${OUTDIR}/ncbi_viral__ALL_FAMILIES.fasta" \
    "${OUTDIR}/ncbi_viral__ALL_FAMILIES.data_report.jsonl" \
    "${OUTDIR}/unclassified.fasta" \
    "${OUTDIR}/ncbi_viral__ALL_FAMILIES__taxonomy.fasta" \
    "${OUTDIR}/unclassified_viroids_taxonomy.fasta" \
    "${OUTDIR}/ncbi_viral_unclassifiedViroids__ALL_FAMILIES__taxonomy_filtered.fasta" \
    "${OUTDIR}/clusters_summary_c1000.tsv" \
    "${OUTDIR}/clusters_mixed_species_members_c1000.tsv" \
    "${OUTDIR}/${DATE}_viralDB_representatives__all__c1.000000.fasta" \
    "${OUTDIR}/${DATE}_viralDB_representatives__all__c0.995000.fasta" \
    "${OUTDIR}/${DATE}_viralDB_representatives__all__c0.990000.fasta" \
    "${OUTDIR}/${DATE}_viralDB_representatives__all__c0.980000.fasta" \
    "${OUTDIR}/${DATE}_viralDB_representatives__all__c0.950000.fasta" 

  finish_manifest
fi

log "✅ Done."
