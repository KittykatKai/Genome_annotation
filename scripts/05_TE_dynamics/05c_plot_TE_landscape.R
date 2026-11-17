# 05c_plot_TE_landscape.R
# Step 4: TE landscape, sequence (Mbp) per 1% divergence bin, stacked by superfamily,
# with estimated age (T = K / 2r) on the bottom axis and % divergence on the top axis.
# Redraw of the course 06-plot_div.R plot with every superfamily coloured and Helitrons removed.
#
# Input   Data/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.anno/
#           Pyl-1.bp.p_ctg.fa.mod.EDTA.RM.nonzero.out.landscape.Div.Rname.tab   (parseRM.pl output)
# Output  Plots/05c_TE_landscape_age.pdf
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

r       <- 8.22e-9   # substitutions per site per year (Kagale et al. 2014)
max_bin <- 40        # show divergence up to 40% (about 24 Mya); hardly anything lies beyond
include_unclassified_LTR <- TRUE   # FALSE drops LTR/unknown, as the course 06-plot_div.R does

# ---- Read the parseRM table (first line is a title, second the header) -----------
tab <- read.delim("Data/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.anno/Pyl-1.bp.p_ctg.fa.mod.EDTA.RM.nonzero.out.landscape.Div.Rname.tab",
                  skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
tab <- tab[!is.na(tab$Rclass) & tab$Rclass != "", ]
bins <- as.matrix(tab[, 4:53])             # bp in each 1% divergence bin, [0;1[ ... [49;50[
bins[is.na(bins)] <- 0

# ---- Group families into readable superfamilies ------------------------------------
fam <- paste(tab$Rclass, tab$Rfam, sep = "/")
group <- ifelse(fam == "LTR/Gypsy",   "Gypsy (LTR)",
         ifelse(fam == "LTR/Copia",   "Copia (LTR)",
         ifelse(fam == "LTR/unknown", "Unclassified LTR",
         ifelse(fam == "DNA/DTC",     "CACTA (TIR)",
         ifelse(fam %in% c("DNA/DTM", "DNA/MULE-MuDR"), "Mutator (TIR)",
         ifelse(fam == "DNA/DTH",     "PIF-Harbinger (TIR)",
         ifelse(tab$Rclass == "DNA" & !grepl("Helitron", fam), "Other TIR (hAT, Tc1, EnSpm)",
         ifelse(tab$Rclass == "MITE", "MITEs",
         ifelse(tab$Rclass == "LINE", "LINE (L1)", NA)))))))))
# Helitrons are removed (EDTA over-calls them), as are unclassified repeats

if (!include_unclassified_LTR) group[group == "Unclassified LTR"] <- NA
keep <- !is.na(group)
mbp  <- rowsum(bins[keep, 1:max_bin], group[keep]) / 1e6   # superfamily x bin, in Mbp

order_groups <- c("Gypsy (LTR)", "Copia (LTR)", "Unclassified LTR", "LINE (L1)",
                  "CACTA (TIR)", "Mutator (TIR)", "PIF-Harbinger (TIR)",
                  "Other TIR (hAT, Tc1, EnSpm)", "MITEs")
order_groups <- order_groups[order_groups %in% rownames(mbp)]
mbp <- mbp[order_groups, , drop = FALSE]

cols <- c("Gypsy (LTR)"                 = "#1F4E9C",
          "Copia (LTR)"                 = "#7FB3E0",
          "Unclassified LTR"            = "#B8C4D6",
          "LINE (L1)"                   = "#8C6239",
          "CACTA (TIR)"                 = "#2E9E5B",
          "Mutator (TIR)"               = "#D7301F",
          "PIF-Harbinger (TIR)"         = "#F4A3A8",
          "Other TIR (hAT, Tc1, EnSpm)" = "#F28E2B",
          "MITEs"                       = "#7B5EA7")[order_groups]

# Bin midpoints: divergence K and age T = K / 2r
k_mid   <- ((1:max_bin) - 0.5) / 100
age_mya <- k_mid / (2 * r) / 1e6

cat("Sequence per superfamily (Mbp, bins 0-", max_bin, "%):\n", sep = "")
print(round(rowSums(mbp), 2))
peak <- which.max(colSums(mbp))
cat(sprintf("Largest bin: %d-%d%% divergence, about %.1f Mya\n", peak - 1, peak, age_mya[peak]))

# ---- Plot ---------------------------------------------------------------------------
show_and_save(function() {
  par(mar = c(4.5, 4.5, 6, 15), xpd = FALSE)
  bp <- barplot(mbp, col = cols, border = NA, space = 0.15, names.arg = rep("", max_bin),
                ylim = c(0, max(colSums(mbp)) * 1.08), las = 1,
                ylab = "Sequence (Mbp)")
  title("TE landscape of Pyl-1", line = 4.2)
  # bottom axis: age in Mya, ticks every 2 Myr
  ages <- seq(0, floor(max(age_mya)), by = 2)
  axis(1, at = approx(age_mya, bp, xout = ages, rule = 2)$y, labels = ages)
  mtext("Estimated age (Mya), T = K / 2r", side = 1, line = 2.8)
  # top axis: % divergence, ticks every 5%
  divs <- seq(0, max_bin, by = 5)
  axis(3, at = approx(k_mid * 100, bp, xout = divs, rule = 2)$y, labels = paste0(divs, "%"),
       cex.axis = 0.85, col.axis = "grey30")
  mtext("Divergence from consensus", side = 3, line = 2.3, cex = 0.85, col = "grey30")
  par(xpd = TRUE)
  legend(x = max(bp) + 1.5, y = max(colSums(mbp)) * 1.05, legend = rev(names(cols)),
         fill = rev(cols), border = NA, bty = "n", cex = 0.85, title = "Superfamily", title.adj = 0)
}, "Plots/05c_TE_landscape_age.pdf", width = 10, height = 5)
