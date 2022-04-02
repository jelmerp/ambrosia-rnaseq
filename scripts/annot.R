#TODO - ALSO CHECK BLASTP RESULTS
#TODO - INCLUDE RRNA ANNOTS FROM RNAMMER GFF?

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
packages <- c("tidyverse", "here")
pacman::p_load(char = packages, install = TRUE)

## Define input files
annot_file <- "results/trinotate/Trinotate.xls"

## Define output files
blax_file <- "results/trinotate/blastx_res.txt"
blax_tax_file <- "results/trinotate/blastx_tax.txt"
blap_file <- "results/trinotate/blastp_res.txt"
blap_tax_file <- "results/trinotate/blastp_tax.txt"
blast_file <- "results/trinotate/blast_res.txt"
blast_tax_file <- "results/trinotate/blast_tax.txt"

## Read input files
annot_raw <- read_tsv(annot_file)

## BLAST column names
blast_colns <- c("gene_name", "gene_name2", "QH", "pct_idty", "eval",
                 "description", "tax")

## Pre-process annotation df
annot <- annot_raw %>%
  rename(gene_id = `#gene_id`,
         blastx = sprot_Top_BLASTX_hit,
         blastp = sprot_Top_BLASTP_hit,
         GO_blastp = gene_ontology_BLASTP,
         GO_blastx = gene_ontology_BLASTX,
         GO_pfam = gene_ontology_Pfam) %>%
  ## Replace missing-data-periods by NA:
  mutate(across(everything(), function(x) sub("^\\.$", NA, x)))


# BLASTX -----------------------------------------------------------------------
## Create a dataframe with the Blastx results
blax <- annot %>%
  select(gene_id, transcript_id, blastx) %>%
  drop_na() %>%
  separate(blastx, sep = "`", fill = "right",
           into = paste0("section", seq_len(max(str_count(.$blastx, "`"), na.rm = TRUE) + 1))) %>%
  pivot_longer(cols = contains("section"), names_to = NULL, values_to = "blastx") %>%
  drop_na() %>%
  separate(blastx, sep = "\\^", fill = "right", into = blast_colns) %>%
  mutate(QH = gsub("[QH]:", "", QH)) %>%
  separate(QH, sep = ",", into = c("q_loc", "h_loc")) %>%
  separate(q_loc, sep = "-", into = c("q_start", "q_end"), remove = FALSE) %>%
  mutate(q_start = as.integer(q_start), q_end = as.integer(q_end)) %>%
  select(-gene_name2) %>%
  mutate(description = sub("RecName: Full=", "", description),
         description = sub(";", "", description),
         pct_idty = as.numeric(sub("%ID", "", pct_idty)),
         eval = as.numeric(sub("E:", "", eval))) %>%
  # TAKE BEST BLAST HIT PER GENE
  group_by(gene_id) %>%
  filter(eval == min(eval)) %>%
  filter(pct_idty == max(pct_idty)) %>%
  filter(q_end - q_start == max(q_end - q_start)) %>%
  select(-q_end, -q_start) %>%
  # REMOVE MULTIPLE ENTRIES FOR SAME GENE
  distinct(gene_id, .keep_all = TRUE) %>%
  select(-transcript_id)

## Create a separate df with just the taxonomic info
blax_tax <- blax %>%
  select(gene_id, tax) %>%
  separate(tax, sep = "; ", fill = "right",
           into = paste0("taxlev", seq_len(max(str_count(.$tax, ";"), na.rm = TRUE) + 1)))

## Check how many genes are fungi / not-fungi / insects
blax_tax %>% filter(taxlev2 == "Fungi") %>% pull(gene_id) %>% unique() %>% length()   # 6,844
blax_tax %>% filter(taxlev2 != "Fungi") %>% pull(gene_id) %>% unique() %>% length()   # 13,062
blax_tax %>% filter(taxlev6 == "Insecta") %>% pull(gene_id) %>% unique() %>% length() # 4,333

## Remove taxonomy from main blast results df
blax <- blax %>% select(-tax)

## Write blastx and blastx-tax df's
write_tsv(blax, blax_file)
write_tsv(blax_tax, blax_tax_file)


# BLASTP -----------------------------------------------------------------------
## Create a dataframe with the Blastp results
blap <- annot %>%
  select(gene_id, transcript_id, blastp) %>%
  drop_na() %>%
  separate(blastp, sep = "`", fill = "right",
           into = paste0("section", seq_len(max(str_count(.$blastp, "`"), na.rm = TRUE) + 1))) %>%
  pivot_longer(cols = contains("section"), names_to = NULL, values_to = "blastp") %>%
  drop_na() %>%
  separate(blastp, sep = "\\^", fill = "right", into = blast_colns) %>%
  mutate(QH = gsub("[QH]:", "", QH)) %>%
  separate(QH, sep = ",", into = c("q_loc", "h_loc")) %>%
  separate(q_loc, sep = "-", into = c("q_start", "q_end"), remove = FALSE) %>%
  mutate(q_start = as.integer(q_start), q_end = as.integer(q_end)) %>%
  select(-gene_name2) %>%
  mutate(description = sub("RecName: Full=", "", description),
         description = sub(";", "", description),
         pct_idty = as.numeric(sub("%ID", "", pct_idty)),
         eval = as.numeric(sub("E:", "", eval))) %>%
  # TAKE BEST BLAST HIT PER GENE
  group_by(gene_id) %>%
  filter(eval == min(eval)) %>%
  filter(pct_idty == max(pct_idty)) %>%
  filter(q_end - q_start == max(q_end - q_start)) %>%
  select(-q_end, -q_start) %>%
  # REMOVE MULTIPLE ENTRIES FOR SAME GENE
  distinct(gene_id, .keep_all = TRUE) %>%
  select(-transcript_id)

## Create a separate df with just the taxonomic info
blap_tax <- blap %>%
  select(gene_id, tax) %>%
  separate(tax, sep = "; ", fill = "right",
           into = paste0("taxlev", seq_len(max(str_count(.$tax, ";"), na.rm = TRUE) + 1)))

## Check how many genes are fungi / not-fungi / insects
blap_tax %>% filter(taxlev2 == "Fungi") %>% pull(gene_id) %>% unique() %>% length()   # 5,066
blap_tax %>% filter(taxlev2 != "Fungi") %>% pull(gene_id) %>% unique() %>% length()   # 10,039
blap_tax %>% filter(taxlev6 == "Insecta") %>% pull(gene_id) %>% unique() %>% length() # 3,186

## Remove taxonomy from main blast results df
blap <- blap %>% select(-tax)

## Write blastx and blastx-tax df's
write_tsv(blap, blap_file)
write_tsv(blap_tax, blap_tax_file)


# COMBINE BLAST RESULTS --------------------------------------------------------
## Take BLASTP results _and_ for genes without a BLASTP result, take the BLASTX result
blast <- rbind(blap, anti_join(blax, blap, by = "gene_id"))
blast_tax <- rbind(blap_tax, anti_join(blax_tax, blap_tax, by = "gene_id"))

## Write to file
write_tsv(blast, blast_file)
write_tsv(blast_tax, blast_tax_file)
