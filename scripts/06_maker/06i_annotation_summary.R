# 06i_annotation_summary.R
# One table of the headline numbers for the presentation: assembly, TEs, genes, AED, BUSCO, BLAST
#
# Inputs  Data/data/Pyl-1.bp.p_ctg.fa.fai
#         Data/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.TEanno.sum
#         Data/MAKER/final/filtered.genes.renamed.gff3
#         Data/MAKER/final/annotation.stat                         (AGAT)
#         Data/MAKER/final/assembly.all.maker.renamed.gff.AED.txt
#         Data/MAKER/final/busco_*/short_summary*.txt
#         Data/comparisons/blastp_uniprot.besthits.tsv, blastp_TAIR10.besthits.tsv
# Output  Tables/06i_annotation_summary.tsv
#
# Open annotation_course.Rproj first, so the working directory is the annotation_course folder.

dir.create("Tables", showWarnings = FALSE)

# ---- Assembly -------------------------------------------------------------
fai <- read.table("Data/data/Pyl-1.bp.p_ctg.fa.fai")
assembly_bp <- sum(fai$V2)

# ---- TEs (total line of TEanno.sum) ---------------------------------------
te_sum  <- readLines("Data/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.TEanno.sum")
te_line <- strsplit(trimws(grep("total interspersed", te_sum, value = TRUE)[1]), "\\s+")[[1]]
te_bp   <- as.numeric(te_line[4])
te_pct  <- te_line[5]

# ---- Genes and mRNAs (filtered GFF) ---------------------------------------
types <- read.table("Data/MAKER/final/filtered.genes.renamed.gff3", sep = "\t", quote = "",
                    colClasses = c("NULL", "NULL", "character", rep("NULL", 6)))$V3
n_genes <- sum(types == "gene")
n_mrna  <- sum(types == "mRNA")   # = number of proteins searched with BLAST

# ---- AGAT statistics (first section of annotation.stat) -------------------
stat_lines <- trimws(readLines("Data/MAKER/final/annotation.stat"))
agat <- function(label) {
  hit <- stat_lines[startsWith(stat_lines, label) &
                    grepl("^[0-9.]+$", trimws(substring(stat_lines, nchar(label) + 1)))][1]
  as.numeric(trimws(substring(hit, nchar(label) + 1)))
}

# ---- AED ------------------------------------------------------------------
aed <- read.table("Data/MAKER/final/assembly.all.maker.renamed.gff.AED.txt", header = TRUE, sep = "\t")
aed_05 <- aed[which.min(abs(aed[, 1] - 0.5)), 2]

# ---- BUSCO ----------------------------------------------------------------
busco_line <- function(dir) {
  f <- list.files(dir, pattern = "^short_summary.*\\.txt$", full.names = TRUE)[1]
  trimws(grep("C:[0-9.]+%\\[S:", readLines(f), value = TRUE)[1])
}

# ---- BLAST ----------------------------------------------------------------
n_uniprot <- length(readLines("Data/comparisons/blastp_uniprot.besthits.tsv"))
n_tair    <- length(readLines("Data/comparisons/blastp_TAIR10.besthits.tsv"))
pct <- function(a, b) paste0(a, " (", round(100 * a / b, 1), "%)")

# ---- Table ----------------------------------------------------------------
summary_tab <- data.frame(
  Metric = c(
    "Assembly size (Mb)", "Contigs", "Longest contig (Mb)",
    "TE content (Mb)", "TE content (% of assembly)",
    "Genes", "mRNAs", "mRNAs per gene (mean)",
    "Single exon genes", "Gene length, mean (bp)", "Gene length, median (bp)",
    "Exons per mRNA (mean)", "CDS length, mean (bp)",
    "Gene models with AED <= 0.5",
    "BUSCO proteins", "BUSCO transcripts",
    "Proteins with UniProt hit", "Proteins with TAIR10 hit"
  ),
  Value = c(
    round(assembly_bp / 1e6, 1), nrow(fai), round(max(fai$V2) / 1e6, 2),
    round(te_bp / 1e6, 2), te_pct,
    n_genes, n_mrna, agat("mean mrnas per gene"),
    pct(agat("Number of single exon gene"), n_genes),
    agat("mean gene length (bp)"), agat("median gene length (bp)"),
    agat("mean exons per mrna"), agat("mean cds length (bp)"),
    paste0(round(100 * aed_05, 1), "%"),
    busco_line("Data/MAKER/final/busco_proteins"),
    busco_line("Data/MAKER/final/busco_transcripts"),
    pct(n_uniprot, n_mrna), pct(n_tair, n_mrna)
  )
)

print(summary_tab, right = FALSE, row.names = FALSE)
write.table(summary_tab, "Tables/06i_annotation_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
