# SETUP ------------------------------------------------------------------------
## Load packages
if (!require(pacman)) install.packages("pacman")
if (!require(BiocManager)) install.packages("BiocManager")
packages <- c("tidyverse", "goseq", "here")
pacman::p_load(char = packages, install = TRUE)

## Source script with functions
source("mcic-scripts/rnaseq/rfuns/enrich_funs.R")

## Define input files
GO_map_file <- here("results/GO/GO_map.txt")
gene_len_file <- here("results/GO/gene_lengths.txt")
DE_file <- here("results/DE/DE_insectgenes.tsv")

## Define output files
GO_outfile <- "results/GO/GO_insectgenes.tsv"

## Read input files
GO_map <- read_tsv(GO_map_file, show_col_types = FALSE) %>%
  as.data.frame() # goseq fails with a tibble
gene_len_df <- read_tsv(gene_len_file, show_col_types = FALSE)
DE <- read_tsv(DE_file, show_col_types = FALSE)


# RUN --------------------------------------------------------------------------
GO_all <- run_GO_all(unique(DE$contrast),
                     DE, GO_map, gene_len_df,
                     use_sig_column = "isDE")

write_tsv(GO_all, GO_outfile)
