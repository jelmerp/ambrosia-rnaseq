# SETUP ------------------------------------------------------------------------
## Load packages
if (!require(pacman)) install.packages("pacman")
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(DESeq2)) BiocManager::install("DESeq2")
packages <- c("tidyverse", "DESeq2", "here")
pacman::p_load(char = packages, install = TRUE)

## Source script with functions
source(here("scripts/DE_funs.R"))

## Define input files
dds_file <- here("results/DE/dds_insect_preDE.rds")
annot_file <- here("results/trinotate/blast_res.txt")

## Define output files
outdir <- here("results/DE/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
DE_file <- here(outdir, "DE_insectgenes.tsv")

## Settings
p_tres <- 0.05
mean_tres <- 1
lfc_tres <- 1

tissue_lvls <- c("DM_S", "IM_S", "IM")
tissue_plot_lvls <- c("DM+S", "IM+S", "IM")
complist <- list(c("DM_S", "IM_S"), c("DM_S", "IM"), c("IM_S", "IM"))


# PREP DATA --------------------------------------------------------------------
## Import DESeq2 object
dds <- readRDS(dds_file)
colData(dds)$tissue <- factor(colData(dds)$tissue, levels = tissue_lvls)
colData(dds)$tissue_plot <- factor(colData(dds)$tissue_plot, levels = tissue_plot_lvls)

## Read annotation data
annot <- read_tsv(annot_file, show_col_types = FALSE)

## Get normalized count data (for computing means)
dds <- estimateSizeFactors(dds)
counts_raw <- dds@assays@data$counts
count_mat <- sweep(counts_raw, 2, sizeFactors(dds), "/")
count_df <- as.data.frame(count_mat) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "id_short", values_to = "count") %>%
  left_join(as.data.frame(colData(dds)), by = "id_short")


# DE ANALYSIS ------------------------------------------------------------------
## DE analysis
design(dds) <- formula(~ tissue)
dds <- DESeq(dds)

## Get DE results for each pairwise comparison
DE_padj <- map_dfr(.x = complist, .f = extract_DE,
                   dds = dds, count_df = count_df, annot = annot,
                   p_tres = 0.05, mean_tres = mean_tres, lfc_tres = lfc_tres) %>%
  rename(isDE_padj = isDE) %>%
  select(-log2FoldChange) # We'll keep the shrunken LFC instead

## Shrunken LFC values + LFC test
DE_lfc <- map_dfr(complist, shrink_lfc,
                  dds = dds, fac = "tissue", lfc_threshold = 1) %>%
  select(gene_id, contrast, log2FoldChange, svalue, isDE_lfc = isDE)

## Combine regular DE analysis with LFC shrinkage
DE <- left_join(DE_padj, DE_lfc, by = c("gene_id", "contrast")) %>%
  mutate(isDE = ifelse(isDE_padj & isDE_lfc, TRUE, FALSE)) %>%
  select(-isDE_lfc, -isDE_padj) %>%
  mutate(group1 = sub("_", "+", group1),
         group2 = sub("_", "+", group2),
         contrast = paste0(group1, " v. ", group2))

## Check numbers of DEGs
DE %>% filter(isDE == TRUE) %>% count(isDE, contrast)

## Write to file
write_tsv(DE, DE_file)
