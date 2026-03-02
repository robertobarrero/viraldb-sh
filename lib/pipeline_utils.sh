#!/usr/bin/env bash
set -euo pipefail

# Nounset-safe defaults for buggy system bashrc files
: "${BASHRCSOURCED:=0}"

# -----------------------------
# Logging
# -----------------------------
ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*"; }
die() { echo "[$(ts)] [ERROR] $*" >&2; exit 1; }

# -----------------------------
# Modes
# -----------------------------
DRY_RUN=0
VALIDATE_ONLY=0

parse_modes() {
  # Accepts: --dry-run and/or --validate
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dry-run) DRY_RUN=1; shift ;;
      --validate) VALIDATE_ONLY=1; shift ;;
      *) shift ;;
    esac
  done
}

# Print + run command, respecting dry-run.
run_cmd() {
  local cmd="$*"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN] ${cmd}"
    return 0
  fi
  log "RUN: ${cmd}"
  eval "${cmd}"
}

# -----------------------------
# Validation helpers
# -----------------------------
require_file() {
  local f="$1"
  [[ -s "$f" ]] || die "Required file missing/empty: $f"
}

require_dir_writable() {
  local d="$1"
  mkdir -p "$d" 2>/dev/null || die "Cannot create directory: $d"
  [[ -w "$d" ]] || die "Directory not writable: $d"
}

require_exec() {
  local exe="$1"
  command -v "$exe" >/dev/null 2>&1 || die "Executable not found on PATH: $exe"
}

# -----------------------------
# Version stamping + checksums
# -----------------------------
init_manifest() {
  # creates ${OUTDIR}/manifest and config snapshot
  require_dir_writable "${OUTDIR}"
  mkdir -p "${OUTDIR}/manifest"

  # snapshot config
  if [[ -f "${CONFIG_FILE:-config_viralDB.txt}" ]]; then
    cp -f "${CONFIG_FILE}" "${OUTDIR}/manifest/config_viralDB.snapshot.txt" || true
  fi

  # git info if available
  if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git rev-parse HEAD > "${OUTDIR}/manifest/git_commit.txt" 2>/dev/null || true
    git status --porcelain > "${OUTDIR}/manifest/git_status.txt" 2>/dev/null || true
  else
    echo "NA" > "${OUTDIR}/manifest/git_commit.txt"
  fi

  # timestamps
  date -Iseconds > "${OUTDIR}/manifest/run_start_time.txt"
}

record_versions() {
  local out="${OUTDIR}/manifest/tool_versions.tsv"
  {
    echo -e "tool\tversion"
    # Best-effort; don't fail if a tool prints oddly.
    echo -e "python\t$(python3 --version 2>&1 | tr -d '\r')"
    echo -e "datasets\t$(datasets version 2>&1 | head -n1 | tr -d '\r')"
    echo -e "cd-hit-est\t$(cd-hit-est -h 2>&1 | head -n1 | tr -d '\r')"
  } > "${out}"
}

checksum_files() {
  # Usage: checksum_files file1 file2 ...
  local out="${OUTDIR}/manifest/sha256sums.tsv"
  local tmp="${OUTDIR}/manifest/sha256sums.tmp.tsv"
  : > "${tmp}"
  echo -e "sha256\tpath" > "${out}"

  for f in "$@"; do
    if [[ -f "$f" ]]; then
      # Use sha256sum if available, else fallback to shasum -a 256 (mac)
      if command -v sha256sum >/dev/null 2>&1; then
        local sum
        sum="$(sha256sum "$f" | awk '{print $1}')"
        echo -e "${sum}\t${f}" >> "${tmp}"
      elif command -v shasum >/dev/null 2>&1; then
        local sum
        sum="$(shasum -a 256 "$f" | awk '{print $1}')"
        echo -e "${sum}\t${f}" >> "${tmp}"
      else
        die "Neither sha256sum nor shasum found; cannot compute checksums"
      fi
    else
      log "[WARN] checksum skip (missing): $f"
    fi
  done

  cat "${tmp}" >> "${out}"
  rm -f "${tmp}"
}

finish_manifest() {
  date -Iseconds > "${OUTDIR}/manifest/run_end_time.txt"
}

#----------------------------------
# Resume pipeline steps
#----------------------------------

# Returns 0 (true) if ALL files exist and are non-empty
files_ok() {
  local f
  for f in "$@"; do
    [[ -s "$f" ]] || return 1
  done
  return 0
}

# Run a step unless outputs already exist and sentinel exists
run_step() {
  local step_name="$1"; shift
  local sentinel="$1"; shift
  local -a outputs=()
  while [[ "$1" != "--" ]]; do
    outputs+=("$1")
    shift
  done
  shift # consume --

  if [[ "${RESUME:-1}" -eq 1 ]] && [[ -f "${sentinel}" ]] && files_ok "${outputs[@]}"; then
    log "⏭️  Skip ${step_name} (resume): outputs + sentinel present"
    return 0
  fi

  log "▶️  Run ${step_name}"
  "$@"

  # Only mark done if command succeeded
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    printf "%s\t%s\n" "$(date '+%F %T')" "${step_name}" > "${sentinel}"
  fi
}


