## Get DE results
get_DE <- function(comp, dds, count_df,
                   annot = NULL, sig_only = FALSE,
                   p_tres = 0.05, mean_tres = 1, lfc_tres = 0.5) {

  ## Get DEseq results
  res <- results(dds, contrast = c("tissue", comp))

  res <- res %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    dplyr::rename(lfc = log2FoldChange) %>%
    mutate(group1 = comp[1],
           group2 = comp[2],
           contrast = paste0(comp, collapse = "_")) %>%
    select(-baseMean) %>%
    arrange(padj)

  ## Include mean normalized counts
  fcount_df <- count_df %>%
      filter(tissue %in% comp) %>%
      group_by(gene_id, tissue) %>%
      summarize(mean = mean(count)) %>%
      pivot_wider(id_cols = gene_id, values_from = mean, names_from = tissue)
  colnames(fcount_df)[2:3] <- c("mean_A", "mean_B")
  fcount_df <- fcount_df %>% mutate(mean = (mean_A + mean_B) / 2)
  res <- res %>% left_join(fcount_df, by = "gene_id")

  ## Determine whether a gene is DE
  res <- res %>% mutate(
    isDE = ifelse(padj < p_tres & abs(lfc) > lfc_tres & (mean_A > mean_tres | mean_B > mean_tres),
                  TRUE, FALSE)
    )
  nsig <- sum(res$isDE, na.rm = TRUE)
  message(comp[1], " vs ", comp[2], " - Nr signif.: ", nsig)

  ## Only take significant genes
  if (sig_only == TRUE) res <- res %>% filter(isDE == TRUE)

  ## Add gene annotation
  if (!is.null(annot)) {
    res <- res %>%
      left_join(annot %>% select(gene_id, gene_name, description),
                res, by = "gene_id") %>%
      arrange(padj)
  }

  return(res)
}

## Volcano plot
pvolc <- function(DE_df,
                  sig_only = TRUE,
                  fcontrast = NULL,
                  cols = NULL,
                  interactive = FALSE) {

  DE_df <- DE_df %>% filter(!is.na(padj))

  if (sig_only == TRUE) DE_df <- DE_df %>% filter(isDE == TRUE)
  if (!is.null(fcontrast)) DE_df <- DE_df %>% filter(contrast %in% fcontrast)

  p <- ggplot(DE_df) +
    aes(x = lfc, y = -log10(padj), color = contrast,
        text = glue("Trinity gene ID: {gene_id}
                    Gene name: {gene_name}
                    Description: {description}")) +
    geom_point() +
    geom_vline(xintercept = 0, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Log-fold change (>0: 1st treatment has higher counts)",
         y = "-log10(adjusted p-value)") +
    guides(color = "none")

  if (is.null(cols)) {
    p <- p + scale_color_brewer(palette = "Dark2")
  } else {
    p <- p + scale_color_manual(values = cols)
  }

  ## When no focal contrast is specified, show all with a facet
  if (is.null(fcontrast)) p <- p + facet_wrap(vars(contrast))

  if (interactive == TRUE) ggplotly(p, tooltip = "text") else print(p)
}


## Function to provide shrunken LFC estimates
shrink_lfc <- function(dds, fac, comp, lfc_threshold = 1) {

  coef <- paste0(fac, "_", paste0(comp, collapse = "_vs_"))
  message("\nCoefficient: ", coef)

  if(! coef %in% resultsNames(dds)) {
    ## For a comparison not in resultsNames
    dds[[fac]] <- relevel(dds[[fac]], ref = comp[2])
    design(dds) <- as.formula(paste("~", fac))
    dds <- DESeq(dds)
  }

  res <- lfcShrink(dds,
                   coef = coef,
                   type = "apeglm",
                   lfcThreshold = lfc_threshold) %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    dplyr::rename(mean = baseMean,
                  lfc = log2FoldChange) %>%
    mutate(group1 = comp[1],
           group2 = comp[2],
           contrast = paste0(comp, collapse = "_"),
           isDE = ifelse(svalue < 0.005, TRUE, FALSE))

  return(res)
}

## Function to make an MA plot
#> Input df should be DE results from get_DE()
pMA <- function(df,
                fcontrast,
                rm_padj_na = TRUE,
                interactive = FALSE) {

  df <- df %>% filter(contrast %in% fcontrast)
  if (rm_padj_na == TRUE) df <- df %>% filter(!is.na(padj))

  p <- ggplot(df) +
    aes(x = mean, y = lfc, color = isDE,
        text = glue("Trinity gene ID: {gene_id}
                    Gene name: {gene_name}
                    Description: {description}")) +
    geom_point(size = 0.5) +
    scale_x_log10(labels = scales::comma_format(accuracy = 1)) +
    scale_color_manual(values = c("grey50", "blue")) +
    guides(color = "none") +
    labs(x = "Mean of normalized counts",
         y = "Mean log-fold change")

  if (length(fcontrast) > 1) {
    p <- p + facet_wrap(vars(contrast))
  }

  if (interactive == TRUE) ggplotly(p, tooltip = "text") else print(p)
}

## Heatmap plot showing abundances
pheat <- function(IDs, count_mat, meta_df,
                  groups = "tissue",
                  show_rownames = TRUE, show_colnames = FALSE,
                  cluster_rows = TRUE,
                  logtrans = TRUE,
                  annotation_colors = NULL,
                  id_labsize = 10, ...) {

  ## Select groups and IDs
  fmeta <- meta_df[, groups, drop = FALSE]
  #fmeta <- fmeta %>% mutate(across(everything(), as.character))

  ## Arrange metadata according to the columns with included factors
  if (length(groups) == 1) fmeta <- fmeta %>% arrange(.data[[groups[1]]])
  if (length(groups) == 2) fmeta <- fmeta %>% arrange(.data[[groups[1]]],
                                                      .data[[groups[2]]])
  if (length(groups) == 3) fmeta <- fmeta %>% arrange(.data[[groups[1]]],
                                                      .data[[groups[2]]],
                                                      .data[[groups[3]]])

  ## Select and arrange count matrix
  fcount_mat <- count_mat[match(IDs, rownames(count_mat)),
                          match(rownames(fmeta), colnames(count_mat))]
  fcount_mat <- as.matrix(fcount_mat)

  ## Log-transform
  if (logtrans == TRUE) {
    fcount_mat <- log10(fcount_mat)
    fcount_mat[fcount_mat == -Inf] <- 0
  }

  ## If few features are included, reduce the cell (row) height
  cellheight <- ifelse(length(IDs) > 20, NA, 20)
  id_labsize <- ifelse(length(IDs) > 40, 6, id_labsize)

  ## Truncate long taxon names
  row.names(fcount_mat) <- str_trunc(row.names(fcount_mat),
                                     width = 20, ellipsis = "")

  ## Function to create the plot
  pheatmap(fcount_mat, annotation_col = fmeta,
           cluster_rows = cluster_rows, cluster_cols = FALSE,
           show_rownames = show_rownames, show_colnames = show_colnames,
           annotation_colors = annotation_colors,
           cellheight = cellheight,
           fontsize = 9, fontsize_row = id_labsize, cex = 1,
           ...)
}

pbox <- function(ID, count_df,
                 annot = NULL, cols = NULL,
                 x_by = "tissue", col_by = "tissue",
                 printplot = TRUE, theme_base_size = 13) {

  if (!is.null(annot)) {
    g_descrip <- annot %>%
      filter(gene_id == ID) %>%
      pull(description) %>%
      str_trunc(width = 50)

    g_name <- annot %>%
      filter(gene_id == ID) %>%
      pull(gene_name)

    psub <- paste0(g_descrip, "\n", ID)
    ptitle <- g_name
  } else {
    ptitle <- ID
    psub <- NULL
  }

  count_df <- count_df %>% filter(gene_id %in% ID)

  p <- ggplot(count_df) +
    aes(y = count, x = .data[[x_by]], color = .data[[x_by]]) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(w = 0.15, h = 0),
               size = 1.5) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.001, 0.2))) +
    labs(title = ptitle,
         subtitle = psub,
         x = NULL,
         y = "Normalized count") +
    guides(color = "none") +
    theme_bw(base_size = theme_base_size) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold.italic"),
          plot.subtitle = element_text(hjust = 0.5, face = "plain"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm"))

  if (is.null(cols)) {
    p <- p + scale_color_brewer(palette = "Set1")
  } else {
    p <- p + scale_color_manual(values = cols)
  }

  if (printplot == TRUE) print(p) else return(p)
}

p4box <- function(IDs, count_df,
                  annot = NULL, cols = NULL,
                  x_by = "tissue", col_by = "tissue") {

  p1 <- pbox(IDs[1], count_df, annot, cols,
             printplot = FALSE, theme_base_size = 10) +
    theme(axis.title.x = element_blank())
          #axis.title.y = element_text(size = 11),
          #plot.title = element_text(size = 11))

  p2 <- pbox(IDs[2], count_df, annot, cols,
             printplot = FALSE, theme_base_size = 10) +
    theme(axis.title = element_blank())
          #plot.title = element_text(size = 11))

  p3 <- pbox(IDs[3], count_df, annot, cols,
             printplot = FALSE, theme_base_size = 10)
    #theme(axis.title = element_text(size = 11),
    #      plot.title = element_text(size = 11))

  p4 <- pbox(IDs[4], count_df, annot, cols,
             printplot = FALSE, theme_base_size = 10) +
    theme(axis.title.y = element_blank())
          #axis.title.x = element_text(size = 11),
          #plot.title = element_text(size = 11))

  p <- (p1 + p2) / (p3 + p4)
  print(p)
}
