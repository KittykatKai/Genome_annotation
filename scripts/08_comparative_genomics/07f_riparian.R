# 07f_riparian.R
# Step 6: riparian plots (synteny braids) from the finished GENESPACE run
# Called by 07g_run_riparian.sh inside genespace_latest.sif

library(GENESPACE)

wd <- "/data/users/ksales/annotation_course/output/comparisons/GENESPACE"
focal_genome <- "Pyl_1"

load(file.path(wd, "results", "gsParams.rda"))   # loads gsParam
outdir <- file.path(wd, "processed_riparian")
dir.create(outdir, showWarnings = FALSE)

# Scaffolds assembled in reverse orientation show as twisted braids.
# After looking at plot 1, list them here, e.g. c("ptg000006l", "ptg000009l"), and rerun.
flip_scaffolds <- c("ptg000006l", "ptg000001l", "ptg000002l", "ptg000012l", "ptg000011l", "ptg000026l", "ptg000005l", "ptg000017l")
flips <- if (length(flip_scaffolds) > 0) data.frame(genome = focal_genome, chr = flip_scaffolds) else NULL

# 1. Pyl_1 vs TAIR10, scaffolds ordered by synteny
plot_riparian(
  gsParam,
  genomeIDs = c("TAIR10", focal_genome),
  refGenome = "TAIR10",
  reorderBySynteny = TRUE,
  pdfFile = file.path(outdir, "Pyl_1_vs_TAIR10.pdf")
)

# 2. Same, with reversed scaffolds flipped for display
plot_riparian(
  gsParam,
  genomeIDs = c("TAIR10", focal_genome),
  refGenome = "TAIR10",
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  invertTheseChrs = flips,
  pdfFile = file.path(outdir, "Pyl_1_vs_TAIR10_oriented.pdf")
)

# 3. All five genomes
plot_riparian(
  gsParam,
  refGenome = "TAIR10",
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  invertTheseChrs = flips,
  scalePlotWidth = 1.5,
  pdfFile = file.path(outdir, "all_genomes_TAIR10.pdf")
)
