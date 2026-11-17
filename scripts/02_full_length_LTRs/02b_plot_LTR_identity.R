# 02b_plot_LTR_identity.R
# Step 1 / 3: LTR identity and insertion age of intact LTR retrotransposons, per clade
#
# Inputs  Data/LTR_full_length/intact_LTR.gff3                   (ltr_identity per element)
#         Data/LTR_full_length/intact_LTR.fa.rexdb-plant.cls.tsv (TEsorter clade per element)
# Outputs Plots/02b_LTR_identity_by_clade.pdf
#         Plots/02b_LTR_age_Copia_Gypsy.pdf
#         Tables/02b_LTR_identity_by_clade.tsv
#
# Open annotation_course.Rproj first, so the working directory is the annotation_course folder.

dir.create("Plots", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)

# Draw a plot in the RStudio Plots pane, then draw it again into a PDF
show_and_save <- function(draw, file, width, height) {
  draw()
  pdf(file, width = width, height = height)
  draw()
  invisible(dev.off())
}

r <- 8.22e-9   # substitutions per site per year (Kagale et al. 2014)

# ---- Intact elements and their LTR identity -------------------------------
gff <- read.table("Data/LTR_full_length/intact_LTR.gff3", sep = "\t", quote = "",
                  stringsAsFactors = FALSE)
gff <- gff[gff$V3 == "repeat_region", ]   # one row per element

ltr <- data.frame(
  name     = sub(".*Name=([^;]+).*", "\\1", gff$V9),
  identity = as.numeric(sub(".*ltr_identity=([^;]+).*", "\\1", gff$V9))
)

# ---- TEsorter clades ------------------------------------------------------
cls <- read.delim("Data/LTR_full_length/intact_LTR.fa.rexdb-plant.cls.tsv", check.names = FALSE)
colnames(cls)[1] <- "TE"

# IDs look like ptg000001l:89522..94244 in the GFF, but minus strand elements are
# written end..start by TEsorter, so both sides are turned into chr:min..max
norm_id <- function(x) {
  x   <- sub("#.*", "", x)
  chr <- sub(":.*", "", x)
  pos <- strsplit(sub(".*:", "", x), "..", fixed = TRUE)
  a   <- as.numeric(sapply(pos, `[`, 1))
  b   <- as.numeric(sapply(pos, `[`, 2))
  paste0(chr, ":", pmin(a, b), "..", pmax(a, b))
}
ltr$key <- norm_id(ltr$name)
cls$key <- norm_id(cls$TE)

ltr <- merge(ltr, cls[, c("key", "Superfamily", "Clade")], by = "key")
ltr <- ltr[ltr$Superfamily %in% c("Copia", "Gypsy") & !ltr$Clade %in% c("unknown", "mixture"), ]
ltr$age_Mya <- (1 - ltr$identity) / (2 * r) / 1e6

cat("Intact LTR-RTs with a clade:", nrow(ltr), "\n")

# ---- Summary table per clade ----------------------------------------------
clade_tab <- do.call(rbind, lapply(split(ltr, ltr$Clade), function(d) data.frame(
  Superfamily       = d$Superfamily[1],
  Clade             = d$Clade[1],
  n                 = nrow(d),
  median_identity   = round(median(d$identity), 4),
  median_age_Mya    = round(median(d$age_Mya), 2),
  pct_younger_1Mya  = round(100 * mean(d$age_Mya < 1), 1)
)))
clade_tab <- clade_tab[order(clade_tab$Superfamily, -clade_tab$n), ]
rownames(clade_tab) <- NULL
print(clade_tab)
write.table(clade_tab, "Tables/02b_LTR_identity_by_clade.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Plot 1: LTR identity per clade ----------------------------------------
sf_cols <- c(Copia = "steelblue", Gypsy = "darkorange")
clade_order <- clade_tab$Clade                               # Copia clades, then Gypsy clades
ltr$Clade <- factor(ltr$Clade, levels = rev(clade_order))
box_cols <- sf_cols[clade_tab$Superfamily[match(levels(ltr$Clade), clade_tab$Clade)]]

show_and_save(function() {
  par(mar = c(4.5, 8, 3, 4))
  boxplot(identity ~ Clade, data = ltr, horizontal = TRUE, las = 1, outline = FALSE,
          col = adjustcolor(box_cols, 0.35), border = box_cols,
          xlab = "Identity between the two LTRs", ylab = "",
          main = "Intact LTR retrotransposons in Pyl-1")
  stripchart(identity ~ Clade, data = ltr, add = TRUE, method = "jitter", jitter = 0.2,
             pch = 16, cex = 0.6, col = adjustcolor(box_cols, 0.7))
  labels_n <- table(ltr$Clade)
  axis(4, at = seq_along(labels_n), labels = paste0("n=", labels_n), las = 1,
       tick = FALSE, line = -0.5, cex.axis = 0.75)
  legend("bottomleft", legend = names(sf_cols), fill = sf_cols, bty = "n", cex = 0.85)
}, "Plots/02b_LTR_identity_by_clade.pdf", width = 8, height = 6)

# ---- Plot 2: insertion age, Copia vs Gypsy ---------------------------------
breaks <- seq(0, ceiling(max(ltr$age_Mya) * 2) / 2, by = 0.25)
counts <- rbind(
  Copia = hist(ltr$age_Mya[ltr$Superfamily == "Copia"], breaks = breaks, plot = FALSE)$counts,
  Gypsy = hist(ltr$age_Mya[ltr$Superfamily == "Gypsy"], breaks = breaks, plot = FALSE)$counts
)
colnames(counts) <- breaks[-length(breaks)]

show_and_save(function() {
  par(mar = c(4.5, 4.5, 3, 1))
  bp <- barplot(counts, beside = TRUE, col = sf_cols[rownames(counts)], border = NA,
                names.arg = rep("", ncol(counts)), space = c(0, 0.4),
                xlab = "Estimated insertion age (Mya), T = K / 2r", ylab = "Intact elements",
                main = "Insertion age of intact LTR retrotransposons")
  whole <- which(breaks[-length(breaks)] %% 1 == 0)          # label every 1 Mya
  axis(1, at = colMeans(bp)[whole] - diff(colMeans(bp))[1] / 2, labels = breaks[whole])
  legend("topright", legend = rownames(counts), fill = sf_cols[rownames(counts)], bty = "n")
}, "Plots/02b_LTR_age_Copia_Gypsy.pdf", width = 8, height = 5)
