#!/bin/bash
#SBATCH --job-name=maker_filter_refine
#SBATCH --partition=pibu_el8
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Filtering and refining the MAKER gene annotation
#   1. Copy merged files to final/ and rename gene IDs
#   2. InterProScan (Pfam) on the proteins
#   3. Add InterProScan results to the GFF
#   4. AED distribution
#   5. Quality filter (AED < 1 and/or Pfam domain)
#   6. Keep gene features only
#   7. Filter transcript and protein FASTA to the retained models

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/ksales/annotation_course
OUTDIR="$WORKDIR/output/MAKER"
FINALDIR="$OUTDIR/final"

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
INTERPROSCAN_SIF="$COURSEDIR/containers/InterProScan_5.70.sif"
ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"

protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

prefix=PYL   # 3 to 4 letter accession prefix: gene IDs become PYL0000001, PYL0000002, ...

mkdir -p "$OUTDIR/tmp" "$FINALDIR/tmp"

# MAKER tools inside the container
maker_run() {
    apptainer exec --cleanenv --env TMPDIR=/tmp \
        --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" \
        --pwd "$FINALDIR" "$MAKER_SIF" "$@"
}

# ---- 1. Copy merged files and rename IDs ----
# Fresh copies every run, so the originals stay unchanged and IDs are never renamed twice
cp "$OUTDIR/$gff"        "$FINALDIR/${gff}.renamed.gff"
cp "$OUTDIR/$protein"    "$FINALDIR/${protein}.renamed.fasta"
cp "$OUTDIR/$transcript" "$FINALDIR/${transcript}.renamed.fasta"
cd "$FINALDIR"

echo "IDs before renaming:"
cut -f9 "${gff}.renamed.gff" | cut -d ';' -f1 | sed 's/ID=//g' | head || true

maker_run maker_map_ids --prefix "$prefix" --justify 7 "${gff}.renamed.gff" > id.map
maker_run map_gff_ids   id.map "${gff}.renamed.gff"
maker_run map_fasta_ids id.map "${protein}.renamed.fasta"
maker_run map_fasta_ids id.map "${transcript}.renamed.fasta"

echo "IDs after renaming:"
cut -f9 "${gff}.renamed.gff" | cut -d ';' -f1 | sed 's/ID=//g' | head || true

# ---- 2. InterProScan (Pfam) ----
apptainer exec --cleanenv --pwd "$FINALDIR" \
    --bind "$COURSEDIR/data/interproscan-5.70-102.0/data":/opt/interproscan/data \
    --bind "$COURSEDIR" \
    --bind "$WORKDIR" \
    --bind "$FINALDIR/tmp":/temp \
    "$INTERPROSCAN_SIF" \
    /opt/interproscan/interproscan.sh \
    -appl pfam --disable-precalc -f TSV \
    --goterms --iprlookup --seqtype p \
    -cpu "$SLURM_CPUS_PER_TASK" \
    -i "${protein}.renamed.fasta" -o output.iprscan

# ---- 3. Add InterProScan results to the GFF ----
maker_run ipr_update_gff "${gff}.renamed.gff" output.iprscan > "${gff}.renamed.iprscan.gff"

# ---- 4. AED distribution ----
maker_run perl /usr/local/bin/AED_cdf_generator.pl -b 0.025 "${gff}.renamed.gff" \
    > assembly.all.maker.renamed.gff.AED.txt

echo "Cumulative fraction of gene models with AED <= 0.5:"
awk 'NR > 1 && ($1 + 0) == 0.5 {print $2}' assembly.all.maker.renamed.gff.AED.txt

# ---- 5. Quality filter: AED < 1 and/or Pfam domain ----
maker_run perl /usr/local/bin/quality_filter.pl -s "${gff}.renamed.iprscan.gff" \
    > "${gff}_iprscan_quality_filtered.gff"

# ---- 6. Keep gene features and subfeatures only ----
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" \
    "${gff}_iprscan_quality_filtered.gff" > filtered.genes.renamed.gff3
sed -i '1i##gff-version 3' filtered.genes.renamed.gff3

echo "Feature types kept:"
cut -f3 filtered.genes.renamed.gff3 | sort | uniq

# ---- 7. Filter transcript and protein FASTA to retained mRNAs ----
grep -P "\tmRNA\t" filtered.genes.renamed.gff3 | awk '{print $9}' | cut -d ';' -f1 | sed 's/ID=//g' > list.txt

apptainer exec --cleanenv --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$ANALYSIS_SIF" \
    seqkit grep -f list.txt "${transcript}.renamed.fasta" -o "${transcript}.renamed.filtered.fasta"
apptainer exec --cleanenv --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$ANALYSIS_SIF" \
    seqkit grep -f list.txt "${protein}.renamed.fasta" -o "${protein}.renamed.filtered.fasta"

# ---- Summary ----
echo
echo "Genes before filtering: $(grep -cP '\tgene\t' "${gff}.renamed.gff")"
echo "Genes after filtering:  $(grep -cP '\tgene\t' filtered.genes.renamed.gff3)"
echo "mRNAs after filtering:  $(wc -l < list.txt)"
echo
echo "Final outputs in $FINALDIR:"
ls -lh filtered.genes.renamed.gff3 \
       "${protein}.renamed.filtered.fasta" \
       "${transcript}.renamed.filtered.fasta" \
       assembly.all.maker.renamed.gff.AED.txt