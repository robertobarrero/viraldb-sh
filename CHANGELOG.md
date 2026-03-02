# Changelog

All notable changes to **viraldb-sh** will be documented in this file.

This project follows:
- [Keep a Changelog](https://keepachangelog.com/)
- [Semantic Versioning](https://semver.org/)

---

## [v0.1] – 2026-03-02

### Added
- **viraldb-sh**: a lightweight, Bash-based pipeline for building a curated plant viral database
- Support for multiple execution environments:
  - Local / interactive shell
  - PBS Pro
  - SLURM
- Config-driven execution via `config_viralDB.txt`
- Validation mode (`--validate`)
  - Verifies required executables on `PATH`
  - Checks input files and output directory permissions
  - Confirms presence of all pipeline scripts
- Dry-run mode (`--dry-run`)
  - Prints all commands without executing them
  - Enables safe inspection and debugging prior to execution
- Automatic detection of run directory:
  - `PBS_O_WORKDIR`
  - `SLURM_SUBMIT_DIR`
  - Current working directory fallback
- Centralised logging per run:
  - Logs written to `logs/viraldb_<DATE>.log`
- Reproducibility features:
  - Manifest generation per run
  - Pipeline version stamping
  - Tool version recording
  - Configuration snapshotting
  - SHA256 checksums for key outputs

### Pipeline workflow
1. Download viral sequences from:
   - NCBI Virus (GenBank + RefSeq)
   - ViroidDB
2. Enrich FASTA headers with taxonomic lineage
3. Merge, filter, and sort sequences
   - Minimum length filtering
   - Ambiguous base (N) filtering
4. Cluster sequences using CD-HIT-EST
   - Identity thresholds: `1.0`, `0.995`, `0.99`
5. Summarise clusters
   - Cluster size statistics
   - Detection of mixed-species clusters
6. Phase 2 representative selection
   - Policy-based, hierarchical representative selection

### Configuration
- User-configurable parameters via `config_viralDB.txt`, including:
  - CPU usage
  - CD-HIT memory allocation (`CDHIT_MEM_MB`)
  - Minimum sequence length
  - Maximum ambiguous base fraction
  - Phase 2 representative selection policy
- Identity thresholds configurable at clustering time

### Notes
- This release prioritises:
  - Transparency
  - Reproducibility
  - HPC safety
  - Ease of debugging
- No workflow engine dependency by design (pure Bash + Python).

---

## [Unreleased]

### Planned
- Optional checksum verification of existing outputs
- Optional container support (Apptainer / Singularity)
- Optional resume / checkpoint support (if kept lightweight)
- Additional cluster summary levels and reporting refinements
