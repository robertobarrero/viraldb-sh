#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------
# Usage:
#   bash run_enrich_taxonomy_per_family.sh /path/to/families /path/to/enrich_headers_with_lineage_from_datasets_jsonl_v0.5.py
#
# Example:
#   bash run_enrich_taxonomy_per_family.sh \
#       20260307_viralDB/families \
#       /viraldb-sh/bin/enrich_headers_with_lineage_from_datasets_jsonl_v0.5.py
# --------------------------------------------------

FAMILIES_DIR="${1:-}"
PYTHON_SCRIPT="${2:-}"

if [[ -z "${FAMILIES_DIR}" || -z "${PYTHON_SCRIPT}" ]]; then
    echo "Usage:"
    echo "  bash $0 <families_dir> <python_script>"
    exit 1
fi

if [[ ! -d "${FAMILIES_DIR}" ]]; then
    echo "[ERROR] Families directory not found: ${FAMILIES_DIR}"
    exit 1
fi

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
    echo "[ERROR] Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi

echo "[INFO] Families directory : ${FAMILIES_DIR}"
echo "[INFO] Python script      : ${PYTHON_SCRIPT}"
echo

N_TOTAL=0
N_DONE=0
N_SKIP=0
N_FAIL=0

for FAMILY_PATH in "${FAMILIES_DIR}"/*; do

    [[ -d "${FAMILY_PATH}" ]] || continue

    FAMILY_NAME="$(basename "${FAMILY_PATH}")"
    N_TOTAL=$((N_TOTAL + 1))

    FASTA="${FAMILY_PATH}/ncbi_viral.fasta"
    JSONL="${FAMILY_PATH}/ncbi_viral_data_report.jsonl"
    OUT_FASTA="${FAMILY_PATH}/${FAMILY_NAME}_ncbi_viral_taxonomy.fasta"

    echo "[INFO] --------------------------------------------------"
    echo "[INFO] Family: ${FAMILY_NAME}"
    echo "[INFO] FASTA : ${FASTA}"
    echo "[INFO] JSONL : ${JSONL}"
    echo "[INFO] OUT   : ${OUT_FASTA}"

    if [[ ! -s "${FASTA}" ]]; then
        echo "[WARN] Missing or empty FASTA, skipping: ${FASTA}"
        N_SKIP=$((N_SKIP + 1))
        continue
    fi

    if [[ ! -s "${JSONL}" ]]; then
        echo "[WARN] Missing or empty JSONL, skipping: ${JSONL}"
        N_SKIP=$((N_SKIP + 1))
        continue
    fi

    if python3 "${PYTHON_SCRIPT}" \
        --fasta "${FASTA}" \
        --datasets_report "${JSONL}" \
        --out_fasta "${OUT_FASTA}" \
        --source NCBI_Virus
    then
        echo "[OK] Finished ${FAMILY_NAME}"
        N_DONE=$((N_DONE + 1))
    else
        echo "[ERROR] Failed ${FAMILY_NAME}"
        N_FAIL=$((N_FAIL + 1))
    fi

    echo
done

echo "[INFO] =================================================="
echo "[INFO] Completed family enrichment"
echo "[INFO] Total families scanned : ${N_TOTAL}"
echo "[INFO] Successful            : ${N_DONE}"
echo "[INFO] Skipped               : ${N_SKIP}"
echo "[INFO] Failed                : ${N_FAIL}"
echo "[INFO] =================================================="
