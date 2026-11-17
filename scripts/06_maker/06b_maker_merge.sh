#!/bin/bash
#SBATCH --job-name=maker_merge
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/ksales/annotation_course
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
OUTDIR="$WORKDIR/output/MAKER"

mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"

# Find the MAKER datastore index (there should be exactly one)
shopt -s nullglob
INDEXES=(*.maker.output/*_master_datastore_index.log)
if [[ ${#INDEXES[@]} -ne 1 ]]; then
    echo "ERROR: expected 1 datastore index, found ${#INDEXES[@]}: ${INDEXES[*]:-none}" >&2
    exit 1
fi
INDEX="${INDEXES[0]}"
echo "Using index: $INDEX"

run_maker() {
    apptainer exec --cleanenv --env TMPDIR=/tmp \
        --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" \
        --pwd "$OUTDIR" "$MAKER_SIF" "$@"
}

# Merge all contigs into one GFF (with sequences), one without sequences, and FASTA files
run_maker gff3_merge -d "$INDEX" -o assembly.all.maker.gff
run_maker gff3_merge -n -d "$INDEX" -o assembly.all.maker.noseq.gff
run_maker fasta_merge -d "$INDEX" -o assembly

ls -lh assembly.all.maker*