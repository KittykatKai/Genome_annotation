#!/bin/bash
#SBATCH --job-name=maker_gene_annotation
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/ksales/annotation_course

cd "$WORKDIR/output/EDTA"

# Genome filename (the TE library is named after it)
genome=Pyl-1.bp.p_ctg.fa

# Step 1: extract Copia and Gypsy sequences from the TE library
apptainer exec -C --bind "$WORKDIR" --pwd "$PWD" "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    seqkit grep -r -p "Copia" "$genome.mod.EDTA.TElib.fa" > Copia_sequences.fa

apptainer exec -C --bind "$WORKDIR" --pwd "$PWD" "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    seqkit grep -r -p "Gypsy" "$genome.mod.EDTA.TElib.fa" > Gypsy_sequences.fa

# Step 2: classify each set into clades with TEsorter
apptainer exec --bind "$WORKDIR" "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
    TEsorter Copia_sequences.fa -db rexdb-plant -p "$SLURM_CPUS_PER_TASK"

apptainer exec --bind "$WORKDIR" "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
    TEsorter Gypsy_sequences.fa -db rexdb-plant -p "$SLURM_CPUS_PER_TASK"