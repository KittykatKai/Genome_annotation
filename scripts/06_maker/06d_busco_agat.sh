#!/bin/bash
#SBATCH --job-name=busco_agat
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Quality assessment of the filtered MAKER annotation
#   1. Longest protein and transcript per gene
#   2. BUSCO (brassicales_odb10) on proteins and transcripts
#   3. AGAT annotation statistics
# Run after 06c_filter_refine_annotation.sh has finished.

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/ksales/annotation_course
FINALDIR="$WORKDIR/output/MAKER/final"

ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
AGAT_SIF="$COURSEDIR/containers/AGAT_v1.4.1.sif"

protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"

cd "$FINALDIR"

# ---- 1. Longest isoform per gene ----
# Gene name = record ID without the isoform suffix (-RA, -RB, ...).
# One record per gene is kept; ties keep the first one seen.
longest_per_gene() {
    local in="$1" out="$2"

    # record ID and sequence length
    awk '/^>/ {if (id) print id "\t" len; id = substr($1, 2); len = 0; next}
         {len += length($0)}
         END {if (id) print id "\t" len}' "$in" > "${out}.lengths.tsv"

    # longest record per gene
    awk -F'\t' '{gene = $1; sub(/[-.]R[A-Z]+$/, "", gene)
                 if (!(gene in best) || $2 > best[gene]) {best[gene] = $2; pick[gene] = $1}}
                 END {for (g in pick) print pick[g]}' "${out}.lengths.tsv" > "${out}.ids.txt"

    apptainer exec --cleanenv --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$ANALYSIS_SIF" \
        seqkit grep -f "${out}.ids.txt" "$in" -o "$out"

    echo "$in: $(grep -c '^>' "$in") records -> $(grep -c '^>' "$out") genes in $out"
}

longest_per_gene "${protein}.renamed.filtered.fasta"    "${protein}.renamed.longest.fasta"
longest_per_gene "${transcript}.renamed.filtered.fasta" "${transcript}.renamed.longest.fasta"

# ---- 2. BUSCO ----
# module scripts can trip over unset variables, so relax -u while loading
set +u
module load BUSCO/5.4.2-foss-2021a
set -u
export PYTHONNOUSERSITE=1

# -f overwrites a previous run's output folder if you resubmit
busco -i "${protein}.renamed.longest.fasta" \
      -l brassicales_odb10 -o busco_proteins -m proteins -c "$SLURM_CPUS_PER_TASK" -f

busco -i "${transcript}.renamed.longest.fasta" \
      -l brassicales_odb10 -o busco_transcripts -m transcriptome -c "$SLURM_CPUS_PER_TASK" -f

echo
echo "BUSCO proteins:";    grep -h "C:" busco_proteins/short_summary*.txt    || true
echo "BUSCO transcripts:"; grep -h "C:" busco_transcripts/short_summary*.txt || true

# ---- 3. AGAT statistics ----
apptainer exec --cleanenv --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$AGAT_SIF" \
    agat_sp_statistics.pl --gff filtered.genes.renamed.gff3 --output annotation.stat

echo
echo "Done. Outputs in $FINALDIR:"
echo "  ${protein}.renamed.longest.fasta, ${transcript}.renamed.longest.fasta"
echo "  busco_proteins/, busco_transcripts/"
echo "  annotation.stat"