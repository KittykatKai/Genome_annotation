# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# Scaffold lengths: keep the 10 largest
fai <- read.table("Pyl-1.bp.p_ctg.fa.fai")[, 1:2]
fai <- head(fai[order(-fai$V2), ], 10)
ideogram <- data.frame(chr = fai$V1, start = 0, end = fai$V2)

# TE annotation: chr, type, start, end, name
gff <- read.table("Pyl-1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3",
                  sep = "\t", quote = "")
gff <- data.frame(chr = gff$V1, type = gff$V3, start = gff$V4, end = gff$V5,
                  name = sub(".*Name=([^;]+).*", "\\1", gff$V9))
gff <- gff[gff$chr %in% ideogram$chr, ]
gff <- gff[!gff$type %in% c("repeat_region", "long_terminal_repeat", "target_site_duplication"), ]

# Function: one circos plot with one density ring per group
plot_circos <- function(groups, cols, file) {
  circos.par(canvas.xlim = c(-1.7, 1.05), canvas.ylim = c(-1.05, 1.05))
  circos.genomicInitialize(ideogram)
  for (i in seq_along(groups)) {
    circos.genomicDensity(groups[[i]][, c("chr", "start", "end")],
                          col = cols[i], track.height = 0.08, window.size = 1e5)
  }
  circos.clear()
  legend("left", legend = names(groups), fill = cols, bty = "n", cex = 0.8,
         inset = 0.01, xpd = TRUE)
  dev.copy(pdf, file, width = 9, height = 9)
  dev.off()
}

# The 6 most abundant TE superfamilies (by bp)
gff <- gff[gff$type != "repeat_fragment", ]
top <- names(head(sort(tapply(gff$end - gff$start, gff$type, sum), decreasing = TRUE), 6))
plot_circos(split(gff, gff$type)[top], rainbow(6), "TE_superfamilies.pdf")

# Athila and CRM clades
cls <- read.delim("TElib.fa.rexdb-plant.cls.tsv")
cls$family <- sub("_(INT|LTR)$", "", sub("#.*", "", cls[, 1]))
gff$family <- sub("_(INT|LTR)$", "", gff$name)

clades <- list(
  Athila = gff[gff$family %in% cls$family[cls$Clade == "Athila"], ],
  CRM    = gff[gff$family %in% cls$family[cls$Clade == "CRM"], ]
)
plot_circos(clades, c("pink", "purple"), "Athila_CRM.pdf")