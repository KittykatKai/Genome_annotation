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
# 07a_prepare_genespace.sh
# Step 1 of the comparative genomics practical: prepare BED + protein FASTA pairs for five genomes
#   Pyl_1 (your MAKER annotation), TAIR10, and three accessions from Lian et al. 2024
# Quick (seconds); run on the login node from the project folder:
#   cd /data/users/ksales/annotation_course && bash scripts/08_comparative_genomics/07a_prepare_genespace.sh

set -euo pipefail

WORKDIR=/data/users/ksales/annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026

Accession=Pyl_1                 # GENESPACE does not allow "-" in names
OTHERS=(Kar_1 Etna_2 Altai_5)       # EDIT: three accessions from Lian et al., not your own

FINALDIR="$WORKDIR/output/MAKER/final"
GENES="$FINALDIR/filtered.genes.renamed.gff3"
LONGEST_PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.longest.fasta"
GENESPACE_DIR="$WORKDIR/output/comparisons/GENESPACE"
LIAN="$COURSEDIR/data/Lian_et_al/GENESPACE"

# ---- Check the chosen accessions exist before copying anything ---------------
problem=0
for other in "${OTHERS[@]}"; do
    if [[ "$other" == "$Accession" ]]; then
        echo "ERROR: $other is your own accession; choose another." >&2; problem=1
    elif [[ ! -s "$LIAN/bed/${other}.bed" || ! -s "$LIAN/peptide/${other}.fa" ]]; then
        echo "ERROR: no BED/FASTA pair for $other in $LIAN" >&2; problem=1
    fi
done
if (( problem )); then
    echo "Available accessions:" >&2
    ls "$LIAN/bed" | sed 's/\.bed$//' | tr '\n' ' ' >&2; echo >&2
    exit 1
fi

mkdir -p "$GENESPACE_DIR/bed" "$GENESPACE_DIR/peptide" "$GENESPACE_DIR/tmp"
cd "$GENESPACE_DIR"
rm -f bed/*.bed peptide/*.fa   # start clean, so files from an earlier choice of accessions don't linger

# ---- Your accession -----------------------------------------------------------
# BED: scaffold, 0-based start, end, gene ID (first attribute of column 9, ID=...)
awk -F'\t' 'BEGIN{OFS="\t"} $3 == "gene" {split($9, a, ";"); split(a[1], b, "="); print $1, $4 - 1, $5, b[2]}' \
    "$GENES" > "tmp/${Accession}.bed"

# FASTA: one (longest) protein per gene; drop descriptions and the -RA/-RB suffix
sed '/^>/s/ .*//' "$LONGEST_PROTEINS" | sed '/^>/s/-R[A-Z]*$//' > "tmp/${Accession}.fa"

# Keep only genes present in both files, so IDs match exactly
grep '^>' "tmp/${Accession}.fa" | sed 's/^>//' | sort -u > tmp/fa_ids.txt
cut -f4 "tmp/${Accession}.bed" | sort -u > tmp/bed_ids.txt
comm -12 tmp/bed_ids.txt tmp/fa_ids.txt > tmp/shared_ids.txt

awk 'NR == FNR {keep[$1]; next} ($4 in keep)' tmp/shared_ids.txt "tmp/${Accession}.bed" \
    > "bed/${Accession}.bed"
awk 'NR == FNR {keep[$1]; next} /^>/ {p = (substr($1, 2) in keep)} p' tmp/shared_ids.txt "tmp/${Accession}.fa" \
    > "peptide/${Accession}.fa"

echo "$Accession: $(wc -l < tmp/bed_ids.txt) genes in GFF, $(wc -l < tmp/fa_ids.txt) proteins, $(wc -l < tmp/shared_ids.txt) kept"

# ---- TAIR10 and three other accessions ---------------------------------------
cp "$COURSEDIR/data/TAIR10.bed" bed/TAIR10.bed
cp "$COURSEDIR/data/TAIR10.fa"  peptide/TAIR10.fa
for other in "${OTHERS[@]}"; do
    cp "$LIAN/bed/${other}.bed"    "bed/${other}.bed"
    cp "$LIAN/peptide/${other}.fa" "peptide/${other}.fa"
done

# ---- Validate: no duplicates, BED and FASTA IDs identical --------------------
echo
printf "%-10s %8s %8s %8s %12s %12s\n" genome genes dupBED dupFASTA missingInFA extraInFA
bad=0
for bedf in bed/*.bed; do
    g=$(basename "$bedf" .bed)
    fa="peptide/${g}.fa"
    cut -f4 "$bedf" | sort > "tmp/${g}.bed.ids"
    grep '^>' "$fa" | sed 's/^>//; s/[[:space:]].*//' | sort > "tmp/${g}.fa.ids"
    dup_bed=$(uniq -d "tmp/${g}.bed.ids" | wc -l)
    dup_fa=$(uniq -d "tmp/${g}.fa.ids" | wc -l)
    missing=$(comm -23 <(sort -u "tmp/${g}.bed.ids") <(sort -u "tmp/${g}.fa.ids") | wc -l)
    extra=$(comm -13 <(sort -u "tmp/${g}.bed.ids") <(sort -u "tmp/${g}.fa.ids") | wc -l)
    printf "%-10s %8s %8s %8s %12s %12s\n" "$g" "$(wc -l < "$bedf")" "$dup_bed" "$dup_fa" "$missing" "$extra"
    (( dup_bed + dup_fa + missing + extra == 0 )) || bad=1
done

echo
if (( bad )); then
    echo "WARNING: some genomes have mismatched or duplicated IDs (see table). Fix before running GENESPACE." >&2
    exit 1
fi
echo "All five genomes OK. Next: sbatch scripts/08_comparative_genomics/07c_run_genespace.sh"
