#!/bin/bash
#SBATCH --job-name=blastp_functional
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch


set -euo pipefail

WORKDIR=/data/users/ksales/annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
FINALDIR="$WORKDIR/output/MAKER/final"
OUTDIR="$WORKDIR/output/comparisons"

PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
GENES="$FINALDIR/filtered.genes.renamed.gff3"
UNIPROT="$COURSEDIR/data/uniprot_viridiplantae_reviewed.fa"
TAIR10="$COURSEDIR/data/TAIR10_pep_20110103_representative_gene_model"
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"

mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"

# The reference databases are already indexed with makeblastdb in the course folder.

# ---- 1. blastp, both searches in parallel ----
set +u
module load BLAST+/2.15.0-gompi-2021a
set -u

blastp -query "$PROTEINS" -db "$UNIPROT" -num_threads 10 -outfmt 6 \
    -evalue 1e-5 -max_target_seqs 10 -out "$OUTDIR/blastp_uniprot.tsv" &
uniprot_job=$!

blastp -query "$PROTEINS" -db "$TAIR10" -num_threads 10 -outfmt 6 \
    -evalue 1e-5 -max_target_seqs 10 -out "$OUTDIR/blastp_TAIR10.tsv" &
tair_job=$!

# wait returns each search's exit code, so set -e stops the job if either fails
wait "$uniprot_job"
wait "$tair_job"

# ---- 2. Best hit per protein (highest bitscore, then lowest e-value) ----
for reference in uniprot TAIR10; do
    sort -k1,1 -k12,12gr -k11,11g "$OUTDIR/blastp_${reference}.tsv" | \
        awk -F '\t' '!seen[$1]++' > "$OUTDIR/blastp_${reference}.besthits.tsv"
done

# ---- 2b. Restore full UniProt IDs (sp|ACCESSION|NAME) ----
# BLAST reports bare accessions, but maker_functional looks descriptions up by the full header ID
awk 'BEGIN{OFS="\t"}
     NR==FNR { if (/^>/) { id = substr($1, 2); n = split(id, a, "|"); acc = (n >= 2 ? a[2] : id); full[acc] = id } next }
     { if ($2 in full) $2 = full[$2]; print }' \
    "$UNIPROT" "$OUTDIR/blastp_uniprot.besthits.tsv" > "$OUTDIR/blastp_uniprot.besthits.fullid.tsv"

# ---- 3. Add UniProt functions to copies of the proteins and genes ----
maker_run() {
    apptainer exec --cleanenv --env TMPDIR=/tmp \
        --bind /data --bind "$WORKDIR" --bind "$OUTDIR/tmp":/tmp \
        --pwd "$OUTDIR" "$MAKER_SIF" "$@"
}

maker_run maker_functional_fasta "$UNIPROT" "$OUTDIR/blastp_uniprot.besthits.fullid.tsv" "$PROTEINS" \
    > "$OUTDIR/maker_proteins.filtered.Uniprot.fasta"

maker_run maker_functional_gff "$UNIPROT" "$OUTDIR/blastp_uniprot.besthits.fullid.tsv" "$GENES" \
    > "$OUTDIR/filtered.genes.Uniprot.gff3"
    
# ---- Summary ----
total=$(grep -c '^>' "$PROTEINS")
uni=$(wc -l < "$OUTDIR/blastp_uniprot.besthits.tsv")
tair=$(wc -l < "$OUTDIR/blastp_TAIR10.besthits.tsv")

echo
echo "Filtered proteins:          $total"
echo "With a UniProt hit:         $uni ($(awk -v a="$uni"  -v b="$total" 'BEGIN{printf "%.1f", 100*a/b}')%)"
echo "With a TAIR10 hit:          $tair ($(awk -v a="$tair" -v b="$total" 'BEGIN{printf "%.1f", 100*a/b}')%)"
echo
echo "Outputs in $OUTDIR:"
ls -lh "$OUTDIR"/blastp_*.tsv "$OUTDIR"/maker_proteins.filtered.Uniprot.fasta "$OUTDIR"/filtered.genes.Uniprot.gff3