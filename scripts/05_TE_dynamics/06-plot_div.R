library(reshape2)
# library(hrbrthemes)
library(tidyverse)
library(data.table)


# get data from parameter

data <- "/data/users/ksales/annotation_course/output/EDTA/Pyl-1.bp.p_ctg.fa.mod.EDTA.anno/Pyl-1.bp.p_ctg.fa.mod.EDTA.RM.nonzero.out.landscape.Div.Rname.tab"

rep_table <- fread(data, header = FALSE, sep = "\t")
rep_table %>% head()
# How does the data look like?

colnames(rep_table) <- c("Rname", "Rclass", "Rfam", 1:50)
rep_table <- rep_table %>% filter(Rfam != "unknown")
rep_table$fam <- paste(rep_table$Rclass, rep_table$Rfam, sep = "/")

table(rep_table$fam)
# How many elements are there in each Superfamily?

rep_table.m <- melt(rep_table)


# Arrange the data so that they are in the following order:
# LTR/Copia, LTR/Gypsy, all types of DNA transposons (TIR transposons), DNA/Helitron, all types of MITES
fam_order <- c(
  "LTR/Copia", "LTR/Gypsy", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron",
  "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM"
)
# Any superfamilies in the data that are not listed above are added at the end
fam_order <- c(fam_order, setdiff(unique(rep_table.m$fam), fam_order))
rep_table.m$fam <- factor(rep_table.m$fam, levels = fam_order)

rep_table.m$distance <- as.numeric(rep_table.m$variable) / 100 # as it is percent divergence

# Question:
# T = K / 2r, with K = divergence and r = 8.22e-9 substitutions per site per year (Kagale et al. 2014)
r <- 8.22e-9
rep_table.m$age <- rep_table.m$distance / (2 * r)   # years; 1% divergence is about 608,000 years


# options(scipen = 999)

# remove helitrons as EDTA is not able to annotate them properly (https://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A)
rep_table.m <- rep_table.m %>% filter(fam != "DNA/Helitron")

ggplot(rep_table.m, aes(fill = fam, x = distance, weight = value / 1000000)) +
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired") +
  xlab("Distance") +
  ylab("Sequence (Mbp)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust = 0.5))

ggsave(filename = "Plots/TE_landscape_divergence.pdf", width = 10, height = 5, useDingbats = F)

# Same landscape with estimated age (million years) on the x axis
ggplot(rep_table.m, aes(fill = fam, x = age / 1e6, weight = value / 1000000)) +
  geom_bar(width = 0.55) +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired") +
  xlab("Estimated age (Mya)") +
  ylab("Sequence (Mbp)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "Plots/TE_landscape_age.pdf", width = 10, height = 5, useDingbats = F)



# Why is it important to have this plot in Mbp instead of counts?
# Hint:
# Consider a scenario where there is a lot of small fragments of TEs due to nested insertions and deletions.
# How would that affect the plot if you used counts instead of Mbp?
#
