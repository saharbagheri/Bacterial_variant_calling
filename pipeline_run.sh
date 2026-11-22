#!/bin/bash
set -o pipefail
cd "$(dirname "$0")"                       # always run from the project root

# Make sure apptainer is available (HPC: load the module)
command -v apptainer >/dev/null || module load apptainer 2>/dev/null
command -v apptainer >/dev/null || { echo "apptainer not found: run this inside Lima or on HPC"; exit 1; }

RUN_DIR=$(python3 -c 'import yaml; print(yaml.safe_load(open("config/config.yaml"))["run_dir"])')
[ -n "$RUN_DIR" ] || { echo "could not read run_dir from config"; exit 1; }

SAMPLES=$(python3 -c 'from pathlib import Path; import yaml; p=Path(yaml.safe_load(open("config/config.yaml"))["samples"]); print(p if p.is_absolute() else (Path.cwd()/p).resolve())')
[ -n "$SAMPLES" ] || { echo "could not read samples path from config"; exit 1; }

mkdir -p "$RUN_DIR/logs"

LOG="$RUN_DIR/logs/run_$(date +%Y%m%d_%H%M%S).log"

echo "Start: $(date)" | tee "$LOG"

snakemake all \
  --workflow-profile workflow/profiles/default \
  --software-deployment-method apptainer \
  --apptainer-prefix "$PWD/.snakemake/apptainer" \
  --snakefile "$PWD/Snakefile" --directory "$RUN_DIR" \
  --configfile "$PWD/config/config.yaml" \
  --config samples="$SAMPLES" \
  --rerun-triggers mtime --cores 4 2>&1 | tee -a "$LOG"

STATUS=$?
echo "End: $(date) | Exit code: $STATUS" | tee -a "$LOG"
exit "$STATUS"



