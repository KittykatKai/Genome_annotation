# 06h_circos_genes_TEs.R
# Step 7 guiding question: how does gene density compare with TE density along the scaffolds?
#   Plot 1: circos plot with rings for genes, all TEs, Gypsy, Copia and the centromeric clades
#   Plot 2: genes vs TE content per 100 kb window, with Spearman correlation
#
# Inputs  Data/data/Pyl-1.bp.p_ctg.fa.fai
#         Data/MAKER/final/filtered.genes.renamed.gff3
#         Data/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3
#         Data/TEsorter_TElib/TElib.fa.rexdb-plant.cls.tsv
# Outputs Plots/06h_circos_genes_TEs.pdf
#         Plots/06h_genes_vs_TEs_windows.pdf
#
# Open annotation_course.Rproj first, so the working directory is the annotation_course folder.

library(circlize)

dir.create("Plots", showWarnings = FALSE)
window <- 1e5   # 100 kb windows

# Draw a plot in the RStudio Plots pane, then draw it again into a PDF
show_and_save <- function(draw, file, width, height) {
  draw()
  pdf(file, width = width, height = height)
  draw()
  invisible(dev.off())
}

# ---- Ideogram: 10 largest scaffolds ---------------------------------------
fai <- read.table("Data/data/Pyl-1.bp.p_ctg.fa.fai")[, 1:2]
fai <- head(fai[order(-fai$V2), ], 10)
ideogram <- data.frame(chr = fai$V1, start = 0, end = fai$V2)

# ---- Genes ----------------------------------------------------------------
genes <- read.table("Data/MAKER/final/filtered.genes.renamed.gff3", sep = "\t", quote = "",
                    colClasses = c("character", "NULL", "character", "integer", "integer",
                                   "NULL", "NULL", "NULL", "NULL"))
colnames(genes) <- c("chr", "type", "start", "end")
genes <- genes[genes$type == "gene" & genes$chr %in% ideogram$chr, ]

# ---- TEs ------------------------------------------------------------------
tes <- read.table("Data/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3", sep = "\t", quote = "")
tes <- data.frame(chr = tes$V1, type = tes$V3, start = tes$V4, end = tes$V5,
                  name = sub(".*Name=([^;]+).*", "\\1", tes$V9))
tes <- tes[tes$chr %in% ideogram$chr &
           !tes$type %in% c("repeat_region", "long_terminal_repeat",
                            "target_site_duplication", "repeat_fragment"), ]

# Centromeric Gypsy clades (Athila, CRM) from TEsorter on the TE library
cls <- read.delim("Data/TEsorter_TElib/TElib.fa.rexdb-plant.cls.tsv")
cls$family <- sub("_(INT|LTR)$", "", sub("#.*", "", cls[, 1]))
tes$family <- sub("_(INT|LTR)$", "", sub("#.*", "", tes$name))
centromeric <- tes[tes$family %in% cls$family[cls$Clade %in% c("Athila", "CRM")], ]

cat("Genes on plotted scaffolds:", nrow(genes), "\n")
cat("TEs on plotted scaffolds:  ", nrow(tes), "\n")

# ---- Plot 1: circos -------------------------------------------------------
rings <- list(
  "Genes"             = genes,
  "All TEs"           = tes,
  "Gypsy"             = tes[tes$type == "Gypsy_LTR_retrotransposon", ],
  "Copia"             = tes[tes$type == "Copia_LTR_retrotransposon", ],
  "Athila + CRM"      = centromeric
)
ring_cols <- c("#1b9e77", "grey35", "darkorange", "steelblue", "purple")

show_and_save(function() {
  circos.par(canvas.xlim = c(-1.7, 1.05), canvas.ylim = c(-1.05, 1.05),
             start.degree = 90, gap.after = 1)
  circos.genomicInitialize(ideogram)
  for (i in seq_along(rings)) {
    circos.genomicDensity(rings[[i]][, c("chr", "start", "end")], col = ring_cols[i],
                          track.height = 0.07, window.size = window)
  }
  circos.clear()
  legend("left", legend = names(rings), fill = ring_cols, bty = "n", cex = 0.85,
         inset = 0.01, xpd = TRUE, title = "Density, 100 kb windows")
}, "Plots/06h_circos_genes_TEs.pdf", width = 9, height = 9)

# ---- Plot 2: genes vs TEs per window --------------------------------------
# Each feature is assigned to the window containing its midpoint
win_key <- function(df) paste(df$chr, ((df$start + df$end) %/% 2) %/% window)

windows <- do.call(rbind, lapply(seq_len(nrow(ideogram)), function(i) {
  n <- ceiling(ideogram$end[i] / window)
  data.frame(key = paste(ideogram$chr[i], 0:(n - 1)))
}))
gene_n  <- table(win_key(genes))
te_bp   <- tapply(tes$end - tes$start + 1, win_key(tes), sum)

windows$genes  <- as.numeric(gene_n[windows$key]);  windows$genes[is.na(windows$genes)]   <- 0
windows$te_pct <- as.numeric(te_bp[windows$key]) / window * 100
windows$te_pct[is.na(windows$te_pct)] <- 0
windows$te_pct <- pmin(windows$te_pct, 100)   # overlapping annotations can exceed 100%

ct <- cor.test(windows$te_pct, windows$genes, method = "spearman", exact = FALSE)
cat(sprintf("Spearman rho (TE %% vs genes per 100 kb): %.2f, p = %.2g, windows = %d\n",
            ct$estimate, ct$p.value, nrow(windows)))

show_and_save(function() {
  par(mar = c(4.5, 4.5, 3, 1))
  plot(windows$te_pct, jitter(windows$genes, 0.6), pch = 16, cex = 0.6,
       col = adjustcolor("grey20", 0.35),
       xlab = "TE content of window (%)", ylab = "Genes per 100 kb window",
       main = "Gene density vs TE density, Pyl-1 (10 largest scaffolds)")
  lines(lowess(windows$te_pct, windows$genes), col = "#1b9e77", lwd = 3)
  legend("topright", bty = "n",
         legend = sprintf("Spearman rho = %.2f, %d windows", ct$estimate, nrow(windows)))
}, "Plots/06h_genes_vs_TEs_windows.pdf", width = 7, height = 5.5)
