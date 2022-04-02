#TODO: SOME GENES ARE NOT PRESENT IN THE GFF FILE!

# SET-UP --------------------------------------------------------------
## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "ape")             # To read GFF
pacman::p_load(char = packages)

## Input files
annot_dir <- here("results/trinotate/")
GO_info_file <- here(annot_dir, "Trinotate.xls.gene_ontology")
gff_file <- here(annot_dir, "trinity_cl0.95.fa.transdecoder.gff3")

## Output files
GO_dir <- here("results", "GO")
if(!dir.exists(GO_dir)) dir.create(GO_dir, recursive = TRUE)
GO_map_file <- here(GO_dir, "GO_map.txt")
gene_len_file <- here(GO_dir, "gene_lengths.txt")

## Read annot files
GO_info <- read_tsv(GO_info_file, col_names = c("gene_id", "GO_terms"))
gff <- read.gff(gff_file)


# CREATE GO MAP ----------------------------------------------------------------
## Define columns for different GO terms for 1 gene:
nmax <- max(str_count(GO_info$GO_terms, ","), na.rm = TRUE) + 1
term_cols <- paste0("term", seq_len(nmax))

GO_map <- GO_info %>%
  separate(GO_terms, into = term_cols, sep = ",", fill = "right") %>%
  pivot_longer(-gene_id, names_to = "term_idx", values_to = "GO_term") %>%
  drop_na() %>%
  select(gene_id, GO_term) %>%
  ## Merge GO terms from different transcripts!
  separate(gene_id, sep = "_i", into = c("gene_id", NA)) %>%
  distinct(gene_id, GO_term, .keep_all = TRUE)

write_tsv(GO_map, GO_map_file)


# CREATE GENE-LENGTH DF FOR GO ANALYSIS ----------------------------------------
gene_len <- gff %>%
  filter(type == "gene") %>%
  mutate(length = end - start + 1) %>%
  select(gene_id = seqid, length) %>%
  arrange(gene_id) %>%
  ## Taking the mean length across transcripts!
  separate(gene_id, sep = "_i", into = c("gene_id", "transcript_nr")) %>%
  group_by(gene_id) %>%
  summarize(length = mean(length))

write_tsv(gene_len, gene_len_file)

## Check if all IDs from the GO map have gene lengths:
#all(GO_map$gene_id %in% gene_len$gene_id)
sum(!unique(GO_map$gene_id) %in% gene_len$gene_id)
sum(unique(GO_map$gene_id) %in% gene_len$gene_id)
