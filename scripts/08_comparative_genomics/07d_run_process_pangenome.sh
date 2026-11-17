#!/bin/bash
#SBATCH --job-name=process_pangenome
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Steps 4 and 5: pan-gene presence frequency (core / accessory / specific) and copy-variable core entries
# Copies the course process_pangenome.R, sets wd and focal_genome, and runs it.
# Submit from the project folder after 07c has finished:
#   sbatch scripts/08_comparative_genomics/07d_run_process_pangenome.sh

set -euo pipefail

WORKDIR=/data/users/ksales/annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
SCRIPTDIR="$WORKDIR/scripts/08_comparative_genomics"
GENESPACE_DIR="$WORKDIR/output/comparisons/GENESPACE"
Accession=Pyl_1
ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
MY_SCRIPT="$SCRIPTDIR/07e_process_pangenome.R"

# Copy the course script once, then set wd and focal_genome
if [[ ! -f "$MY_SCRIPT" ]]; then
    SRC=$(ls "$COURSEDIR"/scripts/*process_pangenome*.R | head -n 1)
    cp "$SRC" "$MY_SCRIPT"
    sed -i -E "s|^([[:space:]]*wd[[:space:]]*<-).*|\1 \"$GENESPACE_DIR\"|" "$MY_SCRIPT"
    sed -i -E "s|^([[:space:]]*focal_genome[[:space:]]*<-).*|\1 \"$Accession\"|" "$MY_SCRIPT"
fi

if ! grep -q "\"$GENESPACE_DIR\"" "$MY_SCRIPT" || ! grep -q "\"$Accession\"" "$MY_SCRIPT"; then
    echo "ERROR: could not set wd/focal_genome automatically. Edit $MY_SCRIPT by hand:" >&2
    grep -nE "wd|focal_genome" "$MY_SCRIPT" | head >&2
    exit 1
fi
echo "Settings used:"; grep -nE "^[[:space:]]*(wd|focal_genome)[[:space:]]*<-" "$MY_SCRIPT"

TMPDIR="$WORKDIR/output/comparisons/Annotation_R_tmp"
mkdir -p "$TMPDIR"

apptainer exec --cleanenv \
    --bind /data --bind "$WORKDIR" \
    --bind "$TMPDIR":/temp --env TMPDIR=/temp --pwd "$WORKDIR" \
    "$ANALYSIS_SIF" Rscript "$MY_SCRIPT"
