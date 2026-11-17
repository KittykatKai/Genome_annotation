dir.create("output/TE_overview", showWarnings = FALSE)

# 1. TEsorter clades for the Copia and Gypsy library sequences
cls <- rbind(
  read.delim("Copia_sequences.fa.rexdb-plant.cls.tsv"),
  read.delim("Gypsy_sequences.fa.rexdb-plant.cls.tsv")
)

# TE_00000123_INT#LTR/Gypsy  ->  TE_00000123  (drop the classification and _INT/_LTR)
cls$family <- sub("_(INT|LTR)$", "", sub("#.*", "", cls[, 1]))

# One clade per family (the _INT part carries the protein domains)
cls <- cls[cls$Clade != "unknown", ]
cls <- cls[!duplicated(cls$family), c("family", "Superfamily", "Clade")]

# 2. Genome-wide TE annotation from EDTA
gff <- read.table("Pyl-1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3",
                  sep = "\t", quote = "")
gff <- data.frame(type = gff$V3, bp = gff$V5 - gff$V4 + 1,
                  name = sub(".*Name=([^;]+).*", "\\1", gff$V9))
gff <- gff[!gff$type %in% c("repeat_region", "long_terminal_repeat", "target_site_duplication"), ]
gff$family <- sub("_(INT|LTR)$", "", sub("#.*", "", gff$name))

# 3. Match annotations to clades
hits <- merge(gff, cls, by = "family")

# 4. Abundance per clade
genome_bp <- sum(read.table("Pyl-1.bp.p_ctg.fa.fai")$V2)

abundance <- aggregate(bp ~ Superfamily + Clade, data = hits, FUN = sum)
abundance$copies     <- aggregate(bp ~ Superfamily + Clade, data = hits, FUN = length)$bp
abundance$Mb         <- round(abundance$bp / 1e6, 2)
abundance$pct_genome <- round(100 * abundance$bp / genome_bp, 2)
abundance <- abundance[order(abundance$Superfamily, -abundance$bp), ]

print(abundance)
write.table(abundance, "Copia_Gypsy_clade_abundance.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Bar plot: % of genome per clade, coloured by superfamily
cols <- ifelse(abundance$Superfamily == "Gypsy", "orange", "steelblue")
par(mar = c(4, 8, 2, 1))
barplot(rev(abundance$pct_genome), names.arg = rev(abundance$Clade), col = rev(cols),
        horiz = TRUE, las = 1, xlab = "% of genome", main = "Copia and Gypsy clades, Pyl-1")
legend("bottomright", legend = c("Gypsy", "Copia"), fill = c("orange", "steelblue"), bty = "n")
dev.copy(pdf, "Copia_Gypsy_clade_abundance.pdf", width = 8, height = 6)
dev.off()