#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Step 3: run GENESPACE (OrthoFinder + DIAMOND + MCScanX) on the five genomes from 07a
# Submit from the project folder:  sbatch scripts/08_comparative_genomics/07c_run_genespace.sh

set -euo pipefail

WORKDIR=/data/users/ksales/annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
SCRIPTDIR="$WORKDIR/scripts/08_comparative_genomics"
GENESPACE_DIR="$WORKDIR/output/comparisons/GENESPACE"
GENESPACE_SIF="$COURSEDIR/containers/genespace_latest.sif"

# Host-backed temporary storage for R's parallel workers
TMPDIR="$WORKDIR/output/comparisons/GENESPACE_R_tmp"
mkdir -p "$TMPDIR"

[[ $(ls "$GENESPACE_DIR"/bed/*.bed | wc -l) -eq 5 ]] || { echo "ERROR: expected 5 BED files; run 07a first" >&2; exit 1; }

cd "$WORKDIR"
apptainer exec --cleanenv \
    --bind /data --bind "$WORKDIR" \
    --bind "$TMPDIR":/temp --env TMPDIR=/temp --pwd "$WORKDIR" \
    "$GENESPACE_SIF" Rscript "$SCRIPTDIR/07b_genespace.R" "$GENESPACE_DIR"

ls -lh "$GENESPACE_DIR/pangenome_matrix.rds" "$GENESPACE_DIR/pangenes" "$GENESPACE_DIR/riparian" 2>/dev/null || true
