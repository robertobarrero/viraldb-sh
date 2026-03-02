# Changelog

All notable changes to **viraldb-sh** will be documented in this file.

This project follows:
- [Keep a Changelog](https://keepachangelog.com/)
- [Semantic Versioning](https://semver.org/)

---

## [v0.1] – 2026-03-02

### Added
- Bash-based viral database construction pipeline (`viraldb-sh`)
- Support for:
  - Local execution
  - PBS Pro
  - SLURM
- Config-driven execution via `config_viralDB.txt`
- Validation mode (`--validate`)
  - Checks required executables
  - Verifies input files and writable directories
- Dry-run mode (`--dry-run`)
  - Prints commands without executing
- Automatic detection of run directory:
  - `PBS_O_WORKDIR`
  - `SLURM_SUBMIT_DIR`
  - Current directory fallback
- Centralised logging per run (`logs/viraldb_<DATE>.log`)
- Manifest generation:
  - Pipeline version
  - Tool versions
  - Configuration snapshot
- SHA256 checksums for key outputs

### Pipeline steps
1. Download NCBI Virus (GenBank + RefSeq) and ViroidDB
2. Enrich FASTA headers with taxonomic lineage
3. Merge, filter, and sort sequences
4. Cluster sequences with CD-HIT-EST
5. Summarise clusters and mixed-species membership
6. Phase 2 representative selection

### Configuration
- User-configurable clustering memory (`CDHIT_MEM_MB`)
- User-configurable CPU usage
- Identity thresholds: `1.0`, `0.995`, `0.99`
- Policy-based Phase 2 representative selection

### Notes
- This release prioritises transparency, debuggability, and HPC safety.
- No workflow engine dependency by design.

---

## [Unreleased]

### Planned
- Optional checksum validation of existing outputs
- Optional container support (Apptainer/Singularity)
