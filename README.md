#Genome annotation of Arabidopsis thaliana Pyl-1

Transposable element (TE) annotation of the Pyl-1 accession, carried out for the Genome and Transcriptome Annotation course (MSc Bioinformatics and Computational Biology, Universities of Bern and Fribourg).

The input is the primary hifiasm assembly of Pyl-1 PacBio HiFi reads, produced in the preceding assembly course. The pipeline identifies and classifies TEs, describes their abundance and distribution along the genome, and estimates when they were active.

All jobs were run on the IBU cluster (Slurm, partition pibu_el8) using the Apptainer containers provided by the course.

## Repository structure
 
```
scripts/
├── 00_setup/                 project folders, GFA to FASTA, genome index
├── 01_TE_annotation_EDTA/    de novo TE library and whole genome TE annotation
├── 02_full_length_LTRs/      clade classification and LTR identity of intact LTR-RTs
├── 03_TE_visualisation/      TE abundance and circos density plots
├── 04_Copia_Gypsy_clades/    genome wide abundance of Copia and Gypsy clades
└── 05_TE_dynamics/           divergence landscape and TE age estimation
```

## Key results so far
 
| | |
|---|---|
| Assembly size | 161.7 Mb in 514 contigs |
| Total TE content | 20.0 Mb (12.4% of the assembly) |
| Gypsy | 6.53 Mb (4.04%) |
| Helitron | 4.43 Mb (2.74%) |
| Unclassified LTR | 2.42 Mb (1.50%) |
| Copia | 2.12 Mb (1.31%) |
| CACTA | 1.23 Mb (0.76%) |
| Mutator | 1.17 Mb (0.72%) |
 
* LTR retrotransposons make up over half of all TE sequence, and Gypsy outweighs Copia roughly threefold.
* Gypsy elements form dense clusters on several scaffolds, which most likely correspond to centromeric regions.
* The divergence landscape shows two waves of TE activity: a large, Gypsy dominated burst around 4 to 6 million years ago and an older wave around 12 to 20 million years ago, in which DNA transposons are more prominent.
* Helitron numbers should be treated with caution, since EDTA is known to overannotate them; they are excluded from the landscape plot.
* Ages are approximate: no CpG correction was applied, and the substitution rate is derived from genes rather than TEs.
## Tools
 
EDTA 2.3.0, RepeatMasker (within EDTA), TEsorter 1.3.0, SeqKit, SAMtools 1.13, parseRM.pl, R 4.4 (circlize, tidyverse, data.table, reshape2, cowplot).
 
## References
 
* Ou S. et al. (2019) Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline. *Genome Biology* 20: 275.
* Zhang R.G. et al. (2022) TEsorter: an accurate and fast method to classify LTR retrotransposons in plant genomes. *Horticulture Research* 9: uhac017.
* Kapusta A., Suh A. and Feschotte C. (2017) Dynamics of genome size evolution in birds and mammals. *PNAS* 114: E1460–E1469.
* Kagale S. et al. (2014) Polyploid evolution of the Brassicaceae during the Cenozoic era. *The Plant Cell* 26: 2777–2791.