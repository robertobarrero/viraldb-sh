#
 viraldb-sh - A pipeline for building a harmonized plant viral database
-

viraldb-sh is a lightweight, Bash-based pipeline for building a curated viral database
from NCBI Virus (GenBank + RefSeq) and ViroidDB, with taxonomy enrichment, sequence 
filtering, clustering, and representative selection.

The pipeline is designed to be:
  •	✅ Simple & transparent (pure Bash calling Python tools)
  •	✅ HPC-safe (PBS Pro, SLURM, or local execution)
  •	✅ Reproducible (version stamping, manifests, checksums)
  •	✅ Fail-fast (validation mode)
  •	✅ Safe to test (dry-run mode)

#--------------------------------------------------------------------#
Pipeline overview
#--------------------------------------------------------------------#

The pipeline performs the following steps:

1. Download viral sequences
   • NCBI Virus (GenBank + RefSeq)
   • ViroidDB

2. Enrich FASTA headers with taxonomic lineage

3. Merge, filter, and sort sequences
   • Remove short sequences
   • Filter sequences with excessive Ns

4. Cluster sequences (CD-HIT-EST)
   • Identity thresholds: 1.0, 0.995, 0.99

5. Summarise clusters
   • Assess mixed-species clustering

6. Phase 2 representative selection
   • Policy-based representative selection

7. Generate checksums & manifest
   • For reproducibility and auditing

#--------------------------------------------------------------------#
Directory structure
#--------------------------------------------------------------------#

viraldb-sh/
├── bin/                     # Python scripts
├── lib/                     # Shared Bash utilities
├── envs/
│   └── viraldb.yaml         # Conda environment definition
├── config_viralDB.txt       # User-editable configuration
├── launch_viralDB_download.sh
├── launch_viralDB_download.pbs
├── launch_viralDB_download.slurm
└── README.md

#--------------------------------------------------------------------#
Requirements
#--------------------------------------------------------------------#

Software:
  • Bash ≥ 4
  • Conda / Mamba
  • Python ≥ 3.9
  • NCBI Datasets CLI (datasets)
  • CD-HIT (cd-hit-est)

All runtime tools are expected to be available via a Conda environment.

#--------------------------------------------------------------------#
Installation
#--------------------------------------------------------------------#

1. Clone or download the pipeline

git clone <repo-url> viraldb-sh
cd viraldb-sh

2. Create the Conda environment

conda env create -f envs/viraldb.yaml

or (recommended)
mamba env create -f envs/viraldb.yaml

#--------------------------------------------------------------------#
Configuration
#--------------------------------------------------------------------#
### Edit the main configuration file

config_viralDB.txt

### Key variables:

PIPELINE_NAME="viraldb-sh"
PIPELINE_VERSION="0.1"

CONDAENV="viralDB"

BIN_DIR="${PIPELINE_DIR}/bin"
#VIRAL_FAMILIES_TSV="/path/to/plant_virus_families.tsv"
VIRAL_FAMILIES_TSV="${BIN_DIR}/plant_virus_families.tsv"

CPUS=4
MIN_LEN=200
MAX_NS_FRACTION=3.0

CDHIT_MEM_MB=8192        # Memory passed to cd-hit-est (-M, in MB)
PHASE2_POLICY="majority"
REFSEQ_MIN_LEN_FRAC=0.90

#--------------------------------------------------------------------#
Running the pipeline
#--------------------------------------------------------------------#

# Local / interactive HPC session
bash launch_viralDB_download.sh

# PBS Pro
qsub launch_viralDB_download.pbs

# SLURM
sbatch launch_viralDB_download.slurm

# note:
The pipeline automatically detects:
  • PBS_O_WORKDIR
  • SLURM_SUBMIT_DIR
  • or defaults to the current directory

Outputs and logs are created in the run directory.

#--------------------------------------------------------------------#
Validation mode (recommended first step prior running the pipeline)
#--------------------------------------------------------------------#

Validation checks:
  • Required executables on PATH
  • Writable output directories
  • Presence of all required input files
  • Presence of pipeline scripts

# Run validation without executing the pipeline:
bash launch_viralDB_download.sh --validate

#--------------------------------------------------------------------#
Dry-run (safe preview)
#--------------------------------------------------------------------#

# Dry-run mode prints all commands without executing them:
bash launch_viralDB_download.sh --dry-run

This is ideal for:
  • Reviewing commands
  • Debugging paths
  • Testing configuration changes

#--------------------------------------------------------------------#
Logging
#--------------------------------------------------------------------#

# Logs are writtenn to:
logs/viraldb_<DATE>.log

# note:
  • Logs are not written during dry-run (clean output).
  • Each run has a unique timestamped log.

#--------------------------------------------------------------------#
Reproducibility and auditing
#--------------------------------------------------------------------#

When not in dry-run mode, the pipeline automatically:
  • Creates a run manifest
  • Records:
    • Pipeline version
    • Tool versions
    • Configuration snapshot
    • Generates checksums (sha256) for key outputs

This enables:
  • Exact reruns (resume from last successful step)
  • Provenance tracking
  • Long-term reproducibility


Feature				Default
--dry-run			off
--validate			off
Manifest + checksums		auto

#--------------------------------------------------------------------#
Output example
#--------------------------------------------------------------------#

20260302_viralDB/
├── ncbi_viral__ALL_FAMILIES.fasta
├── ncbi_viral__ALL_FAMILIES__taxonomy.fasta
├── unclassified_viroids_taxonomy.fasta
├── ncbi_viral_unclassifiedViroids__ALL_FAMILIES__taxonomy_filtered.fasta
├── reps_all/
├── clusters_summary_c1000.tsv
├── clusters_mixed_species_members_c1000.tsv
├── clusters_c1000.tsv
├── representatives__all__c1.000000.fasta
├── representatives__all__c0.995000.fasta
├── representatives__all__c0.990000.fasta
└── manifest/

#--------------------------------------------------------------------#
Why Bash (and not a worflow engine)?
#--------------------------------------------------------------------#

This pipeline intentionally avoids heavy workflow engines to:
  • Reduce cognitive overhead
  • Keep execution transparent
  • Match real HPC usage patterns
  • Enable easy debugging and maintenance

It can be run:
  • Interactively
  • As a batch job
  • Inside larger meta-workflows if needed

#--------------------------------------------------------------------#
Author
#--------------------------------------------------------------------#

Roberto A. Barrero
eResearch, Research Infrastructure – Academic Division
Queensland University of Technology (QUT)


