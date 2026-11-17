# 06f_plot_AED.R
# Step 6: cumulative distribution of AED (Annotation Edit Distance) of the MAKER gene models
#
# Input   Data/MAKER/final/assembly.all.maker.renamed.gff.AED.txt  (AED_cdf_generator.pl output)
# Output  Plots/06f_AED_cumulative.pdf
#
# Open annotation_course.Rproj first, so the working directory is the annotation_course folder.

dir.create("Plots", showWarnings = FALSE)

# Draw a plot in the RStudio Plots pane, then draw it again into a PDF
show_and_save <- function(draw, file, width, height) {
  draw()
  pdf(file, width = width, height = height)
  draw()
  invisible(dev.off())
}

aed <- read.table("Data/MAKER/final/assembly.all.maker.renamed.gff.AED.txt",
                  header = TRUE, sep = "\t")
colnames(aed) <- c("AED", "cum_frac")

frac_05 <- aed$cum_frac[which.min(abs(aed$AED - 0.5))]
cat("Fraction of gene models with AED <= 0.5:", frac_05, "\n")

show_and_save(function() {
  par(mar = c(4.5, 4.5, 3, 1))
  plot(aed$AED, aed$cum_frac, type = "l", lwd = 2.5, col = "darkgreen",
       xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = "AED (0 = full evidence support, 1 = no support)",
       ylab = "Cumulative fraction of gene models",
       main = "AED distribution of Pyl-1 gene models")
  grid(col = "grey90", lty = 1)
  lines(aed$AED, aed$cum_frac, lwd = 2.5, col = "darkgreen")
  segments(0.5, 0, 0.5, frac_05, lty = 2, col = "grey40")
  segments(0, frac_05, 0.5, frac_05, lty = 2, col = "grey40")
  points(0.5, frac_05, pch = 16, col = "darkgreen")
  text(0.52, frac_05 - 0.06, paste0(round(100 * frac_05, 1), "% of models have AED <= 0.5"),
       adj = 0, cex = 0.9)
}, "Plots/06f_AED_cumulative.pdf", width = 7, height = 5.5)
