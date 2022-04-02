## Create named vector of DE genes
get_DE_vec <- function(contrast_id, DE_res) {
  DE_foc <- DE_res %>%
    filter(contrast == contrast_id,
           !is.na(padj)) %>%   # Exclude genes with NA adj-p-val - not tested
    mutate(sig = ifelse(padj < 0.05, 1, 0))

  contrast_vector <- DE_foc$sig
  names(contrast_vector) <- DE_foc$gene_id

  return(contrast_vector)
}

GO_run <- function(contrast_id, DE_vec, GO_map, gene_lens) {

  if(sum(DE_vec > 0)) {
    ## Remove rows from gene length df not in the DE_vec
    fgene_lens <- gene_lens %>% filter(gene_id %in% names(DE_vec))

    ## Remove elements from DE_vec not among the gene lengths
    fDE_vec <- DE_vec[names(DE_vec) %in% fgene_lens$gene_id]

    ## Probability weighting function based on gene lengths
    pwf <- nullp(DEgenes = fDE_vec,
                 bias.data = fgene_lens$length,
                 plot.fit = FALSE)

    ## Run GO test
    GO_df <- goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")

    ## Process GO results
    GO_df <- GO_df %>%
      filter(numDEInCat > 0, # P-adjustment only for genes that were actually tested
             !is.na(ontology)) %>% # Remove GO terms without an ontology and desciption
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             contrast = contrast_id) %>%
      select(contrast,
             p = over_represented_pvalue, padj,
             numDEInCat, numInCat,
             category, ontology, description = term)

    n_sig <- nrow(filter(GO_df, padj < 0.05))
    message("Contrast: ", contrast_id,
            "     Nr DE: ", sum(DE_vec),
            "     Nr sign. GO: ", n_sig)

    return(GO_df)
  } else {
    message("Contrast:", contrast_id, " -- No significant DE genes")
  }
}

GO_wrap <- function(contrast_id, DE_res, GO_map, gene_lens) {
  DE_vec <- get_DE_vec(contrast_id, DE_res)
  GO_df <- GO_run(contrast_id, DE_vec, GO_map, gene_lens)
}

## Prep df for GO plot
prep_goplot <- function(df, contrasts) {

  ## Subset to focal contrasts and tissues
  df <- df %>% filter(contrast %in% contrasts)

  vars <- c("category", "ontology", "description", "tissue", "contrast", "padj")

  df %>%
    select(any_of(vars)) %>%
    ## Pivot wider and then longer to include all terms in all contrasts
    pivot_wider(names_from = contrast, values_from = padj) %>%
    pivot_longer(cols = any_of(contrasts),
                 names_to = "contrast", values_to = "padj") %>%
    left_join(df %>% select(contrast, category, numDEInCat),
              by = c("contrast", "category")) %>%
    ## No labels if not significant
    mutate(numDEInCat = ifelse(padj >= 0.05, NA, numDEInCat)) %>%
    mutate(contrast = sub("padj_", "", contrast),
           padj = ifelse(padj >= 0.05, NA, padj),
           padj_log = -log10(padj)) %>%
    ## Only take GO categories with at least one significant contrast
    filter(category %in% (filter(., padj < 0.05) %>% pull(category))) %>%
    ## Only take contrast with at least one significant category
    filter(contrast %in% (filter(., padj < 0.05) %>% pull(contrast))) %>%
    arrange(padj_log) %>%
    mutate(description = str_trunc(description, width = 45),
           description = ifelse(is.na(description), category, description),
           description = fct_inorder(description))
}

## Heatmap-style plot for significant GO categories
goplot <- function(df, x_var = "contrast",
                   title = NULL, xlabs = NULL,
                   ylabsize = 9, cntlabsize = 1.5,
                   repeat_labs = FALSE) {

  ## Whether to repeat the contrast axis labels for each ontology (i.e. in each facet)
  if (repeat_labs == TRUE) my_scale <- "free" else my_scale <- "free_y"

  p <- df %>%
    ggplot(aes_string(x_var, "description", fill = "padj_log")) +
    geom_tile(stat = "identity", size = 0.25, color = "grey80") +
    geom_label(aes(label = numDEInCat),
               fill = "grey95", size = cntlabsize) +
    scale_fill_viridis_c(option = "D", na.value = "grey95") +
    labs(fill = "-log10\n(adj. p)", title = title) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(legend.position = "left",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = ylabsize))

  ## ggforce::facet_col will keep tile heights constant
  p <- p +
    ggforce::facet_col(vars(ontology),
                       scales = my_scale, space = "free") +
    theme(strip.text = element_text(size = 10, face = "bold"))

  if (!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)

  return(p)
}

## GO dotplot
godotplot <- function(df, type = "GO") {

  #if (type == "GO") group_by <- "ontology" else group_by <- NULL

  p <- ggdotchart(df,
                  x = "description", y = "padj_log",
                  color = "padj_log",
                  sorting = "descending",                       # Sort value in descending order
                  add = "segments",                             # Add segments from y = 0 to dots
                  rotate = TRUE,                                # Rotate vertically
                  #group = group_by,                             # Order by groups
                  dot.size = 5,                                 # Large dot size
                  label = "numDEInCat",                         # Add nr DE genes as dot labels
                  font.label = list(color = "white", size = 9, vjust = 0.5),
                  ggtheme = theme_bw()) +                       # ggplot2 theme
    labs(y = "-log10(adj. p-value)", x = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_color_viridis_c(option = "D", na.value = "grey95") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
          plot.title = element_text(size = 15, face = "bold"),
          strip.text.y = element_text(angle = 270, face = "bold"),
          strip.placement = "outside",
          axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank())

  if (type == "GO") {
    p <- p + facet_grid(ontology~contrast, space = "free", scales = "free")
  } else if (type == "KEGG") {
    p <- p + facet_wrap(vars(contrast), scales = "free_x", nrow = 1)
  }

  print(p)
}

