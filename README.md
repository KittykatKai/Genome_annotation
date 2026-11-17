# Genome annotation of *Arabidopsis thaliana* Pyl-1
 
Transposable element (TE) and gene annotation of the Pyl-1 accession, carried out for the Genome and Transcriptome Annotation course (MSc Bioinformatics and Computational Biology, Universities of Bern and Fribourg).
 
The input is the primary hifiasm assembly of Pyl-1 PacBio HiFi reads, produced in the preceding assembly course. The pipeline identifies, classifies and dates TEs; predicts, filters and functionally annotates protein-coding genes; checks annotation quality; and compares Pyl-1 with TAIR10 and three other accessions from the *A. thaliana* pangenome (Lian et al. 2024).
 
Cluster jobs were run on the IBU cluster (Slurm, partition `pibu_el8`) with the Apptainer containers provided by the course. The R plotting scripts were run locally in RStudio (see [R plotting scripts](#r-plotting-scripts)).
 
## Repository structure
 
```
scripts/
├── 00_setup/                  project folders, GFA to FASTA, genome index
├── 01_TE_annotation_EDTA/     de novo TE library and whole genome TE annotation
├── 02_full_length_LTRs/       clade classification and LTR identity of intact LTR-RTs
├── 03_TE_visualisation/       TE abundance and circos density plots
├── 04_Copia_Gypsy_clades/     genome wide abundance of Copia and Gypsy clades
├── 05_TE_dynamics/            divergence landscape and TE age estimation
├── 06_gene_annotation/        MAKER, ID renaming, InterProScan, AED filtering, BUSCO, AGAT
├── 07_blastp/                 homology to UniProt and TAIR10, functional annotation
└── 08_comparative_genomics/   GENESPACE: orthology, pangenome, synteny (riparian plots)
```
 
Scripts are numbered in the order they are run. `.sh` files are Slurm jobs, submitted from the project root, for example:
 
```bash
cd /data/users/ksales/annotation_course
sbatch scripts/01_TE_annotation_EDTA/01a_run_EDTA.sh
```
 
Data, outputs and logs stay on the cluster and are not tracked in this repository.
 
## Pipeline
 
### 00 Setup
The hifiasm primary contigs are converted from GFA to FASTA and indexed with `samtools faidx`. Only the primary assembly is used, since adding the haplotype assemblies would count each TE and gene several times.
 
### 01 TE annotation with EDTA
EDTA (v2.3.0) builds a non redundant TE library by detecting LTR retrotransposons, TIR transposons, Helitrons, LINEs and SINEs from their structural features, removing gene sequences with the TAIR10 coding sequences. RepeatMasker, run within EDTA, then annotates all TE copies in the genome. The annotation stage was rerun on its own (`--step anno`) with more memory after the first run ran out of memory.
 
### 02 Full length LTR retrotransposons
Intact LTR retrotransposons are classified into superfamilies and clades with TEsorter (v1.3.0, REXdb plant). The identity between the two LTRs of each element is used to estimate insertion age.
 
### 03 TE abundance and distribution
Genomic occupancy per superfamily comes from `TEanno.sum`. Circos plots (circlize) show TE density in 100 kb windows on the largest scaffolds, for the most abundant superfamilies and for the centromeric Gypsy clades Athila and CRM.
 
### 04 Copia and Gypsy clade abundance
Copia and Gypsy sequences are extracted from the TE library with SeqKit, classified with TEsorter and matched to all annotated copies, giving the genome wide abundance of each clade.
 
### 05 TE dynamics and age estimation
RepeatMasker hits with 0% divergence are removed, divergence is binned in 1% classes with parseRM.pl, and the landscape is plotted per superfamily. Divergence K is converted to age with T = K / 2r, using r = 8.22 × 10⁻⁹ substitutions per site per year (Kagale et al. 2014).
 
### 06 Gene annotation (MAKER) and quality
MAKER 3.01.03 combines transcript evidence (Trinity assembly), protein homology (TAIR10 and reviewed UniProt Viridiplantae) and AUGUSTUS ab initio prediction (arabidopsis model) on the genome masked with the EDTA TE library. Outputs are merged, genes are renamed (prefix `PYL`), Pfam domains and GO terms are added with InterProScan 5.70, and models are kept if AED < 1 and/or they carry a Pfam domain. Completeness is checked with BUSCO 5.4.2 (brassicales_odb10) on the longest protein and transcript per gene, and gene structure statistics come from AGAT 1.4.1.
 
### 07 Functional annotation by homology
Filtered proteins are searched with BLASTP (BLAST+ 2.15, e-value 1e-5) against reviewed UniProt Viridiplantae and TAIR10 representative proteins. The best UniProt hit per protein is added to the GFF3 and FASTA with `maker_functional_gff` and `maker_functional_fasta`. BLAST reports bare UniProt accessions, so subject IDs are first mapped back to the full FASTA headers (`sp|ACCESSION|NAME`), which these tools require.
 
### 08 Comparative genomics (GENESPACE)
BED and longest protein files for Pyl-1 (`Pyl_1`), TAIR10 and three accessions from Lian et al. 2024 (Kar_1, Etna_2, Altai_5) are prepared and checked for matching IDs. GENESPACE (OrthoFinder, DIAMOND, MCScanX) builds orthogroups and syntenic blocks. Pan-gene entries are classified as core, accessory or accession specific, copy variable core entries are compared, and riparian plots show synteny. Pyl-1 contigs assembled in reverse orientation are flipped for display.
 
## R plotting scripts
 
These scripts were run locally from an RStudio project whose `Data/` folder mirrors the cluster outputs (for example `Data/MAKER/final/`, `Data/EDTA/`). Each script draws its plot on screen and saves a PDF to `Plots/`.
 
| Script | Plot |
|---|---|
| `02b_plot_LTR_identity.R` | LTR identity per clade and insertion age of intact LTR-RTs |
| `03a_circos_TE_density.R` | Course template: Gypsy and Copia density |
| `03b_circos_TE_superfamilies_clades.R` | Circos plots of the top superfamilies and of Athila and CRM |
| `04b_clade_abundance.R` | Genome wide abundance of Copia and Gypsy clades |
| `05c_plot_TE_landscape.R` | TE landscape with age and divergence axes, all superfamilies coloured |
| `06f_plot_AED.R` | Cumulative AED distribution |
| `06g_plot_BUSCO.R` | BUSCO completeness, proteins and transcripts |
| `06h_circos_genes_TEs.R` | Circos plot of genes vs TEs, and genes vs TE content per 100 kb window |
| `06i_annotation_summary.R` | Table of headline numbers |
 
## Key results
 
**Assembly and TEs**
 
| | |
|---|---|
| Assembly size | 161.7 Mb in 514 contigs (16 contigs > 1 Mb total 140.3 Mb) |
| Total TE content | 20.0 Mb (12.4% of the assembly) |
| Gypsy | 6.53 Mb (4.04%) |
| Helitron | 4.43 Mb (2.74%) |
| Unclassified LTR | 2.42 Mb (1.50%) |
| Copia | 2.12 Mb (1.31%) |
| CACTA | 1.23 Mb (0.76%) |
| Mutator | 1.17 Mb (0.72%) |
 
* LTR retrotransposons make up over half of all TE sequence, and Gypsy outweighs Copia roughly threefold.
* Gypsy, and especially the Athila and CRM clades, forms one dense cluster per large scaffold, marking the centromeres.
* The divergence landscape shows a large Gypsy burst around 5 million years ago and an older wave around 12 to 20 million years ago, in which DNA transposons are more prominent. The youngest bins are dominated by unclassified LTR elements, and most intact LTR-RTs are under 1 million years old.
* Helitrons are excluded from the landscape plot, since EDTA is known to overannotate them. Ages are approximate: no CpG correction was applied, and the substitution rate is derived from genes rather than TEs.
**Genes**
 
| | |
|---|---|
| Genes / mRNAs | 36,732 / 45,331 |
| Single exon genes | 12,037 (32.8%) |
| Mean / median gene length | 1,957 / 1,500 bp |
| Mean exons per mRNA | 4.7 |
| Gene models with AED ≤ 0.5 | 94.7% |
| BUSCO, proteins | C:95.1% [S:93.6%, D:1.5%], F:0.4%, M:4.5%, n:4596 |
| BUSCO, transcripts | C:96.2% [S:92.2%, D:4.0%], F:0.4%, M:3.4%, n:4596 |
| Proteins with a UniProt hit | 36,395 (80.3%) |
| Proteins with a TAIR10 hit | 43,182 (95.3%) |
 
* Gene density decreases as TE content increases (Spearman rho = −0.13 over 100 kb windows): genes lie on the chromosome arms, TEs at the centromeres.
* FLC (AT5G10140) is annotated as PYL0005381 (99% identity to TAIR10, AED 0.08, MADS-box and K-box domains), with no TE insertion within 10 kb.
* The gene count is higher than the 27,000 to 29,000 protein coding genes per accession reported by Lian et al. 2024, consistent with over prediction and the extra contigs below.
**Comparative genomics**
 
| | |
|---|---|
| Pan-gene entries (5 genomes) | 44,096 |
| Core (all 5 genomes) | 23,665 |
| Accessory (2 to 4 genomes) | 11,933 |
| Accession specific | 8,498, of which 6,156 in Pyl-1 |
 
* Every large Pyl-1 contig is fully collinear with one TAIR10 chromosome; there are no large rearrangements. Chromosome 5 is a single contig, and chromosome 2 (which carries rDNA arrays) is the most fragmented.
* About 160 small contigs (6.8 Mb) are mostly chloroplast sequence, carrying up to 150 copies of single chloroplast genes. Most Pyl-1 specific pan-genes lie on small contigs. Without the small contigs, the assembly (140 Mb) matches the published Pyl-1 assembly (about 143 Mb), so purging organelle and duplicate contigs is the main improvement for future work.
## Tools
 
EDTA 2.3.0, RepeatMasker (within EDTA), TEsorter 1.3.0, SeqKit, SAMtools 1.13, parseRM.pl, MAKER 3.01.03, AUGUSTUS, InterProScan 5.70, BUSCO 5.4.2, AGAT 1.4.1, BLAST+ 2.15, GENESPACE (OrthoFinder, DIAMOND, MCScanX), R 4.4 (circlize, GENESPACE, tidyverse, data.table, reshape2, cowplot).
 
## References
 
* Ou S. et al. (2019) Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline. *Genome Biology* 20: 275.
* Zhang R.G. et al. (2022) TEsorter: an accurate and fast method to classify LTR retrotransposons in plant genomes. *Horticulture Research* 9: uhac017.
* Kapusta A., Suh A. and Feschotte C. (2017) Dynamics of genome size evolution in birds and mammals. *PNAS* 114: E1460–E1469.
* Kagale S. et al. (2014) Polyploid evolution of the Brassicaceae during the Cenozoic era. *The Plant Cell* 26: 2777–2791.
* Lian Q. et al. (2024) A pan-genome of 69 *Arabidopsis thaliana* accessions reveals a conserved genome structure throughout the global species range. *Nature Genetics* 56: 982–991.
* Lovell J.T. et al. (2022) GENESPACE tracks regions of interest and gene copy number variation across multiple genomes. *eLife* 11: e78526.
 