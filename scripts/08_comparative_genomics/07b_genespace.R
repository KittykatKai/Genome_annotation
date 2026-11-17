# 07b_genespace.R
# Step 2: initialise and run GENESPACE (called by 07c_run_genespace.sh inside genespace_latest.sif)
# Usage: Rscript 07b_genespace.R <GENESPACE working directory>

library(GENESPACE)

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]   # folder that contains bed/ and peptide/

gpar <- init_genespace(wd = wd, path2mcscanx = "/usr/local/bin", nCores = 20)

out <- run_genespace(gpar)

pangenome <- query_pangenes(out, bed = NULL, refGenome = "TAIR10", transform = TRUE,
                            showArrayMem = TRUE, showNSOrtho = TRUE, maxMem2Show = Inf)

saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))
