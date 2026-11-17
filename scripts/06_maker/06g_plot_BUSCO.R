# 06g_plot_BUSCO.R
# Step 7: BUSCO completeness of the annotation (longest protein and transcript per gene)
#
# Inputs  Data/MAKER/final/busco_proteins/short_summary*.txt
#         Data/MAKER/final/busco_transcripts/short_summary*.txt
# Output  Plots/06g_BUSCO.pdf
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

read_busco <- function(dir) {
  file <- list.files(dir, pattern = "^short_summary.*\\.txt$", full.names = TRUE)[1]
  line <- trimws(grep("C:[0-9.]+%\\[S:", readLines(file), value = TRUE)[1])
  num  <- function(tag) as.numeric(sub(paste0(".*", tag, ":([0-9.]+).*"), "\\1", line))
  list(summary = line,
       values  = c(S = num("S"), D = num("D"), F = num("F"), M = num("M")))
}

busco <- list(
  Proteins    = read_busco("Data/MAKER/final/busco_proteins"),
  Transcripts = read_busco("Data/MAKER/final/busco_transcripts")
)
for (b in names(busco)) cat(b, ":", busco[[b]]$summary, "\n")

mat  <- sapply(busco, `[[`, "values")
cols <- c(S = "#56B4E9", D = "#3492C7", F = "#F0E442", M = "#F04442")
cats <- c("Complete, single copy", "Complete, duplicated", "Fragmented", "Missing")

show_and_save(function() {
  par(mar = c(8, 7, 3, 1))
  bp <- barplot(mat, horiz = TRUE, col = cols, border = "white", las = 1, xlim = c(0, 100),
                main = "BUSCO completeness of the Pyl-1 annotation")
  title(xlab = "% of BUSCOs (brassicales_odb10, n = 4596)", line = 2.3)
  text(2, bp, sub(",n:.*", "", sapply(busco, `[[`, "summary")), adj = 0, cex = 0.8, font = 2)
  legend("bottom", inset = c(0, -0.62), xpd = TRUE, horiz = TRUE, bty = "n",
         legend = cats, fill = cols, cex = 0.8)
}, "Plots/06g_BUSCO.pdf", width = 8, height = 4)
