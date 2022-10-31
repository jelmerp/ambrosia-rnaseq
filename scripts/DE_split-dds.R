# SETUP ------------------------------------------------------------------------
## Load packages
if (!require(pacman)) install.packages("pacman")
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(DESeq2)) BiocManager::install("DESeq2")
packages <- c("tidyverse", "DESeq2", "here")
pacman::p_load(char = packages, install = TRUE)

## Define input files
dds_file <- here("results/kallisto/beetle/deseq_object.rds")
tax_file <- here("results/trinotate/blast_tax.txt")
annot_file <- here("results/trinotate/blast_res.txt")

## Define output files
outdir <- here("results/DE")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dds_insect_file <- here(outdir, "dds_insect_preDE.rds")
dds_fungal_file <- here(outdir, "dds_fungal_preDE.rds")
dds_nonfungal_file <- here(outdir, "dds_nonfungal_preDE.rds")

## Settings
tissue_lvls <- c("DM_S", "IM_S", "IM")
tissue_plot_lvls <- c("DM+S", "IM+S", "IM")


# PREP DATA --------------------------------------------------------------------
## Load taxonomic and gene info
tax <- read_tsv(tax_file, show_col_types = FALSE)
annot <- read_tsv(annot_file, show_col_types = FALSE)

## Import DESeq2 object
dds <- readRDS(dds_file)
colData(dds)$tissue <- factor(colData(dds)$tissue, levels = tissue_lvls)
colData(dds)$tissue_plot <- factor(colData(dds)$tissue_plot, levels = tissue_plot_lvls)

# CREATE SEPARATE DESEQ OBJECTS ------------------------------------------------
## Only select genes with an insect annotation
insect_genes <- filter(tax, taxlev6 == "Insecta") %>% pull(gene_id) %>% unique()
fungal_genes <- tax %>% filter(taxlev2 == "Fungi") %>% pull(gene_id) %>% unique()
nonfungal_genes <- tax %>% filter(taxlev2 != "Fungi") %>% pull(gene_id) %>% unique()

## Create 3 different DESeq objects
dds_insect <- dds[insect_genes, ]
dds_fungal <- dds[fungal_genes, ]
dds_nonfungal <- dds[nonfungal_genes, ]

## Write to file
saveRDS(dds_insect, dds_insect_file)
saveRDS(dds_fungal, dds_fungal_file)
saveRDS(dds_nonfungal, dds_nonfungal_file)
