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
genome=Pyl-1.bp.p_ctg.fa

R_CONTAINER="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
ANNODIR="$WORKDIR/output/EDTA/$genome.mod.EDTA.anno"
PLOT_SCRIPT="$WORKDIR/scripts/06-plot_div.R"
LTR_GFF="$WORKDIR/output/LTR_full_length/intact_LTR.gff3"
OUTDIR="$WORKDIR/output/TE_age"

RATE=8.22e-9   # substitutions per synonymous site per year

mkdir -p "$OUTDIR"

# Stop early (before the slow parseRM step) if the plot script is missing
[[ -f "$PLOT_SCRIPT" ]] || { echo "ERROR: $PLOT_SCRIPT not found. Save the edited 06-plot_div.R there." >&2; exit 1; }

cd "$ANNODIR"

# ---- Steps 1 and 2: drop 0% divergence hits, then bin with parseRM.pl ----
awk 'NR == 1 || $2 != 0' "$genome.mod.EDTA.RM.out" > "$genome.mod.EDTA.RM.nonzero.out"

apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$PWD" "$R_CONTAINER" \
    perl "$COURSEDIR/scripts/05-parseRM.pl" \
    -i "$genome.mod.EDTA.RM.nonzero.out" \
    -l 50,1 -v

echo "parseRM output files:"
ls -lh "$genome.mod.EDTA.RM.nonzero.out"*

# ---- Step 3: landscape plots (runs straight after parseRM) ----
mkdir -p Plots
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$PWD" "$R_CONTAINER" \
    Rscript --vanilla "$PLOT_SCRIPT"

# ---- Step 4a: age for each 1% divergence bin, T = K / 2r ----
awk -v r="$RATE" 'BEGIN{
    OFS = "\t"
    print "divergence_pct", "age_years", "age_Mya"
    for (d = 1; d <= 50; d++) {
        t = (d / 100) / (2 * r)
        print d, sprintf("%.0f", t), sprintf("%.2f", t / 1e6)
    }
}' > "$OUTDIR/divergence_bin_ages.tsv"

# ---- Step 4b: age of each intact LTR retrotransposon, K = 1 - LTR identity ----
awk -F'\t' -v r="$RATE" 'BEGIN{
    OFS = "\t"
    print "seqid", "start", "end", "classification", "ltr_identity", "K", "age_years", "age_Mya"
}
$0 !~ /^#/ && $9 ~ /ltr_identity=/ {
    cls = "NA"; id = ""
    n = split($9, attr, ";")
    for (i = 1; i <= n; i++) {
        split(attr[i], kv, "=")
        if (kv[1] == "Classification") cls = kv[2]
        if (kv[1] == "ltr_identity")   id  = kv[2]
    }
    if (id > 1) id = id / 100          # in case identity is given as a percentage
    k = 1 - id
    t = k / (2 * r)
    print $1, $4, $5, cls, id, sprintf("%.4f", k), sprintf("%.0f", t), sprintf("%.3f", t / 1e6)
}' "$LTR_GFF" > "$OUTDIR/intact_LTR_ages.tsv"

# Short summary in the log
awk -F'\t' 'NR > 1 {n++; s += $8; if ($8 < 1) young++} END{
    printf "Intact LTR-RTs: %d\nMean age: %.2f Mya\nYounger than 1 Mya: %d (%.1f%%)\n",
           n, s / n, young, 100 * young / n
}' "$OUTDIR/intact_LTR_ages.tsv"

echo "Done."
echo "Landscape plot:  $ANNODIR/Plots/"
echo "Age tables:      $OUTDIR/"