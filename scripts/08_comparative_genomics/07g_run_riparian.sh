#!/bin/bash
#SBATCH --job-name=riparian
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Step 6: draw the riparian plots with 07f_riparian.R
# Submit from the project folder after 07c has finished:
#   sbatch scripts/08_comparative_genomics/07g_run_riparian.sh

set -euo pipefail

WORKDIR=/data/users/ksales/annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
SCRIPTDIR="$WORKDIR/scripts/08_comparative_genomics"
GENESPACE_SIF="$COURSEDIR/containers/genespace_latest.sif"

TMPDIR="$WORKDIR/output/comparisons/GENESPACE_R_tmp"
mkdir -p "$TMPDIR"

cd "$WORKDIR"
apptainer exec --cleanenv \
    --bind /data --bind "$WORKDIR" \
    --bind "$TMPDIR":/temp --env TMPDIR=/temp --pwd "$WORKDIR" \
    "$GENESPACE_SIF" Rscript "$SCRIPTDIR/07f_riparian.R"

ls -lh "$WORKDIR/output/comparisons/GENESPACE/processed_riparian/"
