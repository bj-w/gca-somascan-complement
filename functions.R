# COLOURS -----------------------------------------------------------------
col.dx <- c("non-GCA" = "#66C2A5", "GCA" = "#FC8D62")

# Plot theme --------------------------------------------------------------
figTheme <- function(layout = "line", fontsize = 10, font = "Helvetica") {
  require(ggtext)
  line_theme <- theme_classic(base_size = fontsize) +
    theme(
      plot.title = ggtext::element_markdown(
        face = "plain",
        hjust = 0.5,
        size = fontsize + 1
      ),
      plot.subtitle = ggtext::element_markdown(size = fontsize),
      plot.tag = element_text(size = 14, face = "bold", family = font),
      # legend.title = ggtext::element_markdown(),
      axis.line = element_line(colour = "black"),
      axis.title.x = ggtext::element_markdown(size = fontsize),
      axis.title.y = ggtext::element_markdown(size = fontsize),
      axis.text = element_text(colour = "black", size = fontsize),
      strip.text = ggtext::element_markdown(size = fontsize),
      strip.background = element_rect(fill = NA, colour = NA),
      legend.text = element_text(size = fontsize - 1),
      legend.title = ggtext::element_markdown(size = fontsize)
    )
  if (layout == "box") {
    box_theme <- line_theme +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_blank()
      )
    return(box_theme)
  } else {
    return(line_theme)
  }
}


# Format P-value ----------------------------------------------------------
formatPvalLabel <- function(p, digits = 1, includeP = TRUE) {
  if (p < 2.2e-16) {
    if (includeP == TRUE) {
      paste0("<i>P</i> < 2.2x10<sup>-16</sup>")
    } else {
      paste0("< 2.2x10<sup>-16</sup>")
    }
  } else if (p < 0.001) {
    if (includeP == TRUE) {
      sprintf(
        paste0("<i>P</i> = %.", digits, "fx10<sup>%d</sup>"),
        p / 10^floor(log10(abs(p))),
        floor(log10(abs(p)))
      )
    } else {
      sprintf(
        paste0("%.", digits, "fx10<sup>%d</sup>"),
        p / 10^floor(log10(abs(p))),
        floor(log10(abs(p)))
      )
    }
  } else {
    if (includeP == TRUE) {
      paste0("<i>P</i> = ", signif(p, digits = digits))
    } else {
      paste0(signif(p, digits = digits))
    }
  }
}


# QC ----------------------------------------------------------------------

outlierPCA <- function(
  prcompObject,
  sdLimit = 3,
  fillGroup = outlier,
  outlierLabelSize = 3
) {
  require(ggrepel)
  # needed for plotting required group
  fillGroup <- enquo(fillGroup)

  # calculate upper and lower limits for PCs
  calculate_limits <- function(values) {
    mean_val <- mean(values)
    sd_val <- sd(values)
    lower_limit <- mean_val - sdLimit * sd_val
    upper_limit <- mean_val + sdLimit * sd_val
    list(lower = lower_limit, upper = upper_limit)
  }
  pc1_limits <- calculate_limits(prcompObject$x[, 1])
  pc2_limits <- calculate_limits(prcompObject$x[, 2])

  # calculate % variance explained by each PC
  perc_variance_explained <- summary(prcompObject)$importance[2, ] * 100

  plot_df <- as.data.frame(prcompObject$x) %>%
    rownames_to_column("study_id") %>%
    left_join(meta, by = "study_id") %>%
    mutate(
      outlier = case_when(
        (PC1 > pc1_limits$upper | PC1 < pc1_limits$lower) |
          (PC2 > pc2_limits$upper | PC2 < pc2_limits$lower) ~ "Outlier",
        .default = "Not outlier"
      )
    )

  outlier_samples <- plot_df %>%
    filter(outlier == "Outlier") %>%
    pull(study_id)

  p <- ggplot(plot_df, aes(x = PC1, y = PC2, fill = {{ fillGroup }})) +
    geom_blank() + # adds both 'Outlier' and 'Not outlier' to legend if even no outliers
    geom_hline(
      yintercept = c(unlist(pc2_limits)),
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = c(unlist(pc1_limits)),
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_point(shape = 21, colour = "black", size = 2) +
    # geom_text_repel(
    #   data = filter(plot_df, outlier == "Outlier"),
    #   aes(label = study_id),
    #   bg.colour = "white",
    #   size = outlierLabelSize,
    #   min.segment.length = 0
    # ) +
    labs(
      x = paste0("PC1 (", round(perc_variance_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(perc_variance_explained[2], 2), "%)")
    ) +
    scale_fill_manual(
      values = c("Outlier" = "red", "Not outlier" = "grey75"),
      drop = FALSE
    )

  output <- list(data = plot_df, outliers = outlier_samples, plot = p)
  return(output)
}

# PCA biplot
plotPCABiplot <- function(
  prcompObject,
  metadata,
  color_var,
  x_pc,
  y_pc,
  sample_id_col = "SampleID",
  marginalDensity = TRUE
) {
  require(ggside)
  # For ggplot
  color_var_sym <- sym(color_var)
  x_pc_sym <- sym(x_pc)
  y_pc_sym <- sym(y_pc)

  # calculate proportion of variance explained
  prop_explained <- data.frame(
    prop_explained = summary(prcompObject)$importance[2, ]
  ) %>%
    rownames_to_column("PC")
  x_prop <- prop_explained %>%
    filter(PC == x_pc) %>%
    pull(prop_explained)
  y_prop <- prop_explained %>%
    filter(PC == y_pc) %>%
    pull(prop_explained)
  # plot
  plot <- prcompObject$x %>%
    as.data.frame() %>%
    rownames_to_column(sample_id_col) %>%
    left_join(
      unique(metadata %>% select(all_of(c(sample_id_col, color_var)))),
      by = sample_id_col
    ) %>%
    ggplot(aes(x = !!x_pc_sym, y = !!y_pc_sym)) +
    geom_point(
      aes(fill = !!color_var_sym),
      shape = 21,
      size = 2.5,
      colour = "black"
    ) +
    labs(
      x = paste0(x_pc, " (", round(x_prop * 100, 2), "% variance)"),
      y = paste0(y_pc, " (", round(y_prop * 100, 2), "% variance)")
    )
  if (marginalDensity == TRUE) {
    plot <- plot +
      geom_xsidedensity(aes(fill = !!color_var_sym), alpha = 0.75) +
      geom_ysidedensity(aes(fill = !!color_var_sym), alpha = 0.75) +
      theme_ggside_void() +
      theme(ggside.axis.line = element_blank())
  }

  return(plot)
}

# outlier plot based on median and IQR
outlierMedianIQR <- function(
  counts_df,
  fillGroup = outlier,
  metadata,
  sdLimit,
  outlierLabelSize = 3
) {
  fillGroup <- enquo(fillGroup)
  require(ggrepel)
  outlier <- counts_df %>%
    summarise(
      median = median(log2_RFU),
      iqr = IQR(log2_RFU),
      .by = study_id
    )

  # calculate median and IQR RFU values based on 3 SDs from the mean
  calculate_limits <- function(values) {
    mean_val <- mean(values)
    sd_val <- sd(values)
    lower_limit <- mean_val - sdLimit * sd_val
    upper_limit <- mean_val + sdLimit * sd_val
    list(lower = lower_limit, upper = upper_limit)
  }
  median_limits <- calculate_limits(outlier$median)
  iqr_limits <- calculate_limits(outlier$iqr)

  # plot
  plot_df <- outlier %>%
    mutate(
      outlier = case_when(
        median > median_limits$upper | median < median_limits$lower ~ "Outlier",
        iqr > iqr_limits$upper | iqr < iqr_limits$lower ~ "Outlier",
        .default = "Not outlier"
      )
    ) %>%
    left_join(metadata, by = "study_id")

  outlier_samples <- plot_df %>%
    filter(outlier == "Outlier") %>%
    pull(study_id)

  p <- plot_df %>%
    ggplot(aes(x = median, y = iqr, fill = {{ fillGroup }})) +
    geom_blank() + # adds both 'Outlier' and 'Not outlier' to legend if even no outliers
    geom_point(shape = 21, colour = "black", size = 2) +
    geom_vline(
      xintercept = c(median_limits$lower, median_limits$upper),
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = c(iqr_limits$lower, iqr_limits$upper),
      linetype = "dashed",
      alpha = 0.5
    ) +
    # geom_text_repel(
    #   data = plot_df %>% filter(outlier == "Outlier"),
    #   aes(label = study_id),
    #   bg.colour = "white",
    #   min.segment.length = 0,
    #   size = outlierLabelSize
    # ) +
    labs(x = "Median log<sub>2</sub> RFU", y = "IQR log<sub>2</sub> RFU") +
    scale_fill_manual(
      values = c("Outlier" = "red", "Not outlier" = "grey75"),
      drop = FALSE
    )

  return(list(plot = p, outliers = outlier_samples, data = plot_df))
}


# function to run linear regression
lmPCA <- function(metricDF, prcompObject, uniqueIdColumn, PC, title) {
  lmPCA_output <- map(1:PC, function(pc) {
    map(names(metricDF[-1]), function(column) {
      lm_df <- left_join(
        data.frame(pc = prcompObject$x[, pc]) %>%
          rownames_to_column(uniqueIdColumn),
        metricDF,
        by = uniqueIdColumn
      )
      glance(lm(pc ~ ., data = lm_df %>% select(pc, {{ column }}))) %>%
        mutate(metric = column, pc = pc)
    })
  }) %>%
    list_flatten() %>%
    list_rbind() %>%
    filter(metric != "id") %>%
    mutate(pc = paste0("PC", pc))

  prop_explained <- summary(prcompObject)[["importance"]][2, ] %>%
    data.frame(prop_explained = .) %>%
    rownames_to_column("pc")

  lmPCA_output <- lmPCA_output %>%
    left_join(prop_explained, by = "pc") %>%
    relocate(pc, metric)
  return(lmPCA_output)
}

plotlmPCA <- function(lmPCA_output, plotTitle = NULL) {
  # plot screeplot
  screeplot <- lmPCA_output %>%
    select(pc, prop_explained) %>%
    distinct() %>%
    ggplot(aes(
      y = fct_reorder(pc, prop_explained, .desc = TRUE),
      x = prop_explained * 100
    )) +
    geom_col(aes(fill = prop_explained), colour = "black") +
    scale_x_reverse(expand = expansion(mult = c(0.2, 0))) +
    scale_y_discrete(limits = rev) +
    scale_fill_gradient(
      limits = c(0, 1),
      low = "white",
      high = "steelblue",
      na.value = NA
    ) +
    labs(y = NULL, x = "Variance<br>explained<br>(%)", title = plotTitle) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  # plot heatmap
  heatmap <- lmPCA_output %>%
    mutate(r.squared = case_when(p.value > 0.05 ~ NA, .default = r.squared)) %>%
    ggplot(aes(
      y = fct_reorder(pc, prop_explained, .desc = TRUE),
      x = metric,
      fill = r.squared,
      label = round(r.squared, 2)
    )) +
    geom_tile(colour = "black") +
    geom_text(size = 2.5) +
    scale_x_discrete(
      expand = c(0, 0),
      guide = guide_axis(angle = 45),
      position = "top"
    ) +
    scale_y_discrete(expand = c(0, 0), limits = rev) +
    scale_fill_gradient(
      limits = c(0, 1),
      low = "white",
      high = "steelblue",
      na.value = NA,
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
    ) +
    labs(x = NULL, y = NULL, fill = "R<sup>2</sup>") +
    figTheme(layout = "box") +
    theme(legend.position = "right")
  # combine plots
  (screeplot | heatmap) + plot_layout(widths = c(1, 6))
}

# LIMMA -------------------------------------------------------------------
# simple generic function to run limma using log2 somascan RFU
runLimma <- function(
  metadata,
  counts,
  design,
  contrastMatrix = NULL,
  numericOutcome = NULL,
  featureInfo = NULL,
  featureId,
  outputCsvPath = NULL
) {
  require(edgeR)
  require(limma)

  if (is.null(contrastMatrix) & is.null(numericOutcome)) {
    stop("Specify at least a contrast matrix (if factor) or numeric outcome")
  }

  # remove samples with missing numeric outcomes
  if (!is.null(numericOutcome)) {
    # set up for numeric outcomes
    # identify samples missing numeric outcome
    missing_outcome <- row.names(metadata[is.na(metadata[[numericOutcome]]), ])
    message(paste0(
      "Excluding ",
      length(missing_outcome),
      " samples with missing numeric outcome data."
    ))
    # exclude from metadata
    metadata <- metadata[!row.names(metadata) %in% missing_outcome, ]
    # exclude from counts
    counts <- counts[, !colnames(counts) %in% missing_outcome]
  }

  # check metadata and counts align
  if (!identical(row.names(metadata), colnames(counts))) {
    counts <- counts[, match(row.names(metadata), colnames(counts))]
    stopifnot(identical(row.names(metadata), colnames(counts)))
  }

  # FIT LINEAR MODELS
  vfit <- lmFit(counts, design)
  # add contrast matrix if provided
  if (!is.null(contrastMatrix)) {
    vfit <- contrasts.fit(vfit, contrasts = contrastMatrix)
  }

  # EBAYES SHRINKING
  efit <- eBayes(vfit, robust = TRUE, trend = TRUE)

  # COLLECT RESULTS
  if (is.null(contrastMatrix)) {
    results <- topTable(efit, number = Inf, coef = numericOutcome) # coef is the numeric outcome
  } else {
    results <- topTable(efit, number = Inf, coef = 1) # selects first comparison from contrast matrix by default
  }
  # extract results
  results <- results %>%
    rownames_to_column(featureId) %>%
    as_tibble()
  # add on SOMAmer info if provded
  if (!is.null(featureInfo)) {
    results <- results %>% left_join(featureInfo, by = featureId)
  }
  # export csv file if path if provided
  if (!is.null(outputCsvPath)) {
    write_csv(results, outputCsvPath)
    message(paste0("Limma output: ", outputCsvPath))
  }
  # return all results objects
  return(list(vfit = vfit, efit = efit, results = results))
}


# GSEA --------------------------------------------------------------------

# prepare gene sets using msigdbr
prepareGeneSets <- function(
  geneIdentifier,
  hallmark = FALSE,
  kegg = FALSE,
  reactome = FALSE,
  go_bp = FALSE
) {
  require(dplyr)
  require(msigdbr)
  gene_sets <- list()
  if (isTRUE(hallmark)) {
    message("Preparing Hallmark gene sets...")
    # MSigDB Hallmark
    hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
      distinct(gs_name, !!sym(geneIdentifier))
    ## remove HALLMARK_ and replace underscore with space
    hallmark <- hallmark %>%
      mutate(
        gs_name = str_remove(gs_name, "HALLMARK_"),
        gs_name = str_replace_all(gs_name, "_", " ")
      )
    gene_sets$hallmark <- hallmark
  }

  # KEGG
  if (isTRUE(kegg)) {
    message("Preparing KEGG Medicus gene sets...")
    kegg <- msigdbr(
      species = "Homo sapiens",
      collection = "C2",
      subcollection = "CP:KEGG_MEDICUS"
    ) %>%
      distinct(gs_name, !!sym(geneIdentifier))
    kegg <- kegg %>%
      mutate(
        gs_name = str_remove(gs_name, "KEGG_MEDICUS_"),
        gs_name = str_replace_all(gs_name, "_", " ")
      )
    gene_sets$kegg <- kegg
  }

  # REACTOME
  if (isTRUE(reactome)) {
    message("Preparing Reactome gene sets...")
    reactome <- msigdbr(
      species = "Homo sapiens",
      collection = "C2",
      subcollection = "CP:REACTOME"
    ) %>%
      distinct(gs_name, !!sym(geneIdentifier))
    reactome <- reactome %>%
      mutate(
        gs_name = str_remove(gs_name, "REACTOME_"),
        gs_name = str_replace_all(gs_name, "_", " ")
      )
    gene_sets$reactome <- reactome
  }
  # GENE ONTOLOGY - BIOLOGICAL PROCESSES
  if (isTRUE(go_bp)) {
    message("Preparing Gene Ontology (Biological Processes) gene sets...")
    go_bp <- msigdbr(
      species = "Homo sapiens",
      collection = "C5",
      subcollection = "GO:BP"
    ) %>%
      distinct(gs_name, !!sym(geneIdentifier))
    go_bp <- go_bp %>%
      mutate(
        gs_name = str_remove(gs_name, "GOBP_"),
        gs_name = str_replace_all(gs_name, "_", " ")
      )
    gene_sets$go_bp <- go_bp
  }
  if (length(gene_sets) < 1) {
    stop(" No gene sets returned. Did you specify any gene sets?")
  } else {
    return(gene_sets)
  }
}

# add somamers with multiple gene ids to gene set if there's a consensus
addComplexEntrezIdsToGeneSet <- function(differentialAbundance, geneSet) {
  id_mapping <- differentialAbundance |> select(apt_id, EntrezGeneID)
  complex_gene_ids <- id_mapping |>
    # just get somamers with multiple gene IDs. Separated by |
    filter(str_detect(EntrezGeneID, "\\|")) |>
    separate_longer_delim(cols = EntrezGeneID, delim = "|") |>
    left_join(
      geneSet,
      by = c("EntrezGeneID" = "ncbi_gene"),
      relationship = "many-to-many"
    ) |>
    # count number of times each gene set appears for each somamer
    add_count(gs_name, .by = apt_id, name = "gs_n") |>
    # count number of unique gene IDs for each somamer
    mutate(gene_n = n_distinct(EntrezGeneID), .by = apt_id) |>
    # if number of times gene set occurs is the same as the number of gene IDs
    # for each somamer, there is a consensus
    filter(gs_n == gene_n) |>
    # formatting
    select(apt_id, gs_name) |>
    left_join(id_mapping, by = "apt_id") |>
    select(gs_name, ncbi_gene = EntrezGeneID)
  # add complex gene IDs to original gene sets
  geneSet |>
    bind_rows(complex_gene_ids) |>
    arrange(gs_name)
}

# generic function to run GSEA
runGSEA <- function(results, geneSet, featureIdentifier = "EntrezGeneID") {
  require(fgsea)
  filt <- results %>%
    # rank according to descending t stats. high = upregulated in GCA
    arrange(desc(t))
  t_stats <- filt$t
  names(t_stats) <- filt[[featureIdentifier]]
  # prepare pathways
  geneSet <- split(geneSet[[2]], geneSet[[1]])
  # GSEA
  set.seed(1997)
  gsea_obj <- fgsea(
    pathways = geneSet,
    stats = t_stats,
    eps = 0,
    nPermSimple = 10000,
    nproc = 1
  )
  gsea_obj <- arrange(gsea_obj, pval)
  return(list(pathways = geneSet, stats = t_stats, gseaOutput = gsea_obj))
}

# function to creating GSEA plot
# modified fgsea::plotEnrichment
plotGSEA_ES <- function(
  geneset,
  pathways,
  stats,
  gseaOutput,
  lineColour = "green",
  gseaParam = 1,
  ticksSize = 0.2,
  statsSize = 3.5,
  subtitleSize = 10
) {
  require(fgsea)
  require(ggtext)
  pathway <- pathways[[geneset]]
  pd <- plotEnrichmentData(
    pathway = pathway,
    stats = stats,
    gseaParam = gseaParam
  )
  # plot
  p <- with(
    pd,
    ggplot(data = curve) +
      geom_line(aes(x = rank, y = ES), color = lineColour) +
      geom_segment(
        data = ticks,
        mapping = aes(
          x = rank,
          y = -spreadES / 16,
          xend = rank,
          yend = spreadES / 16 # , colour = rank
        ),
        linewidth = ticksSize
      ) +
      # scale_colour_pretty_div("RedBlues", midpoint = median(pd$stats$rank), limits = c(min(pd$stats$rank), max(pd$stats$rank)), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      geom_hline(yintercept = posES, colour = "red", linetype = "dashed") +
      geom_hline(yintercept = negES, colour = "red", linetype = "dashed") +
      geom_hline(yintercept = 0, colour = "black") +
      labs(x = "Rank", y = "Enrichment score", subtitle = geneset) +
      theme(plot.subtitle = element_markdown(size = subtitleSize))
  )

  # add on NES and p-value
  # figure out where to put stats on the plot
  if (gseaOutput[gseaOutput$pathway == geneset, ][["NES"]] < 0) {
    # NES is negative so put stats in bottom left
    x <- -Inf
    y <- -Inf
    vjust <- -0.5
    hjust <- 0
  } else {
    # NES is positive so put stats in top right
    x <- Inf
    y <- Inf
    vjust <- 1.5
    hjust <- 1
  }
  # add to plot
  p +
    geom_richtext(
      data = gseaOutput %>% filter(pathway == {{ geneset }}),
      aes(
        x = x,
        y = y,
        label = paste0(
          "NES = ",
          round(NES, 2),
          "<br>",
          formatPvalLabel(pval)
        )
      ),
      vjust = vjust,
      hjust = hjust,
      fill = NA,
      label.color = NA,
      size = statsSize
    )
}

plotGSEA_heatmap <- function(gseaDf, filterNoneSig = FALSE, plotTitle = NULL) {
  gseaDf <- gseaDf %>% mutate(asterisk = if_else(padj < 0.05, "*", NA))
  if (!isFALSE(filterNoneSig)) {
    gseaDf <- gseaDf %>% filter(padj < 0.05)
  }
  ggplot(gseaDf, aes(x = "x", y = fct_reorder(pathway, NES, .desc = FALSE))) +
    geom_tile(aes(fill = NES)) +
    geom_text(
      aes(label = asterisk),
      colour = "white",
      fontface = "bold",
      size = 5
    ) +
    labs(x = NULL, y = NULL, fill = "NES", title = plotTitle) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), labels = scales::label_wrap(25)) +
    scale_fill_gradient2(
      low = "#00344aff",
      mid = "#f5f5f5ff",
      high = "#b53737ff",
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
    ) +
    # scale_fill_pretty_c(palette = "RedBlues", direction=-1, guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    figTheme("box") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
}

# volcano plot ------------------------------------------------------------

plotVolcano <- function(resultsDf, labelN = 0, title = NULL) {
  require(ggrepel)
  # features to label
  to_label <- resultsDf %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(P.Value) %>%
    slice_head(n = labelN) %>%
    pull(apt_id)
  df <- resultsDf %>%
    mutate(
      colour = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "up",
        adj.P.Val < 0.05 & logFC < 0 ~ "down",
        .default = "notSig"
      ),
      colour = factor(colour, levels = c("up", "down", "notSig")),
      label = if_else(apt_id %in% to_label, EntrezGeneSymbol, NA)
    )
  plot <- df %>%
    ggplot(aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(colour = colour)) +
    geom_text_repel(aes(label = label), size = 3, min.segment.length = 0) +
    annotate(
      geom = "text",
      x = Inf,
      y = Inf,
      label = nrow(filter(df, colour == "up")),
      hjust = 2,
      vjust = 1,
      colour = "#B53737",
      fontface = "bold",
      size = 4
    ) +
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      label = nrow(filter(df, colour == "down")),
      hjust = -2,
      vjust = 1,
      colour = "#00344a",
      fontface = "bold",
      size = 4
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_colour_manual(
      values = c("up" = "#B53737", down = "#00344a", notSig = "grey80")
    ) +
    labs(
      x = "log<sub>2</sub>(FC)",
      y = "-log<sub>10</sub>(<i>P</i>)",
      title = title
    ) +
    theme(legend.position = "none")
  return(plot)
}


# PROTEIN CLASSIFICATION ENRICHMENT ---------------------------------------
# function to run TissueEnrich

# format input genes and background sets
formatTissueEnrichGeneSets <- function(limmaResultDf) {
  # requires unique gene IDs, so keeping most sig SOMAmer per ensembl gene ID.
  limmaResultDf <- limmaResultDf %>%
    arrange(P.Value) %>%
    distinct(ensembl_id, .keep_all = TRUE) %>%
    drop_na(ensembl_id)
  # genes of interest are differentially abundant proteins
  gs <- GeneSet(
    geneIds = filter(limmaResultDf, adj.P.Val < 0.05)$ensembl_id,
    organism = "Homo Sapiens",
    geneIdType = ENSEMBLIdentifier()
  )
  # background is all proteins tested for differential abundance
  bg <- GeneSet(
    geneIds = limmaResultDf$ensembl_id,
    organism = "Homo Sapiens",
    geneIdType = ENSEMBLIdentifier()
  )
  return(list(inputGenes = gs, backgroundGenes = bg))
}

runTissueEnrich <- function(
  reference,
  tpmMatrix = NULL,
  inputGenes,
  backgroundGenes,
  categoryIdentifier,
  plotTitle = NULL
) {
  require(TissueEnrich)
  if (!reference %in% c("custom", 1, 2, 3)) {
    stop(
      "reference argument must be 'custom' / 1 (HPA) / 2 (GTEx) / 3 (Mouse ENCODE)"
    )
  }
  if (!reference == "custom") {
    enrich_output <- teEnrichment(
      inputGenes = inputGenes,
      rnaSeqDataset = reference,
      backgroundGenes = backgroundGenes,
      # will perform manually after removing
      # categories with no specific features
      multiHypoCorrection = FALSE
    )
    enrich_results <- enrich_output[[1]] %>%
      assay() %>%
      as.data.frame()

    enrich_categories <- map(names(enrich_output[[3]]), \(x) {
      pluck(enrich_output, 3, x) %>%
        assay() %>%
        as_tibble() %>%
        mutate(Tissue = x)
    }) %>%
      # remove empty dfs from list
      discard(\(x) nrow(x) == 0) %>%
      list_rbind()
  } else {
    # convert to summarized experiment
    se <- SummarizedExperiment(
      assays = SimpleList(as.matrix(tpmMatrix)),
      rowData = row.names(tpmMatrix),
      colData = colnames(tpmMatrix)
    )
    # enrichment categorisation
    enrich_categories <- teGeneRetrieval(se)
    # run enrichment
    enrich_results <- teEnrichmentCustom(
      inputGenes = inputGenes,
      tissueSpecificGenes = enrich_categories,
      backgroundGenes = backgroundGenes,
      # will perform manually after removing
      # categories with no specific features
      multiHypoCorrection = FALSE
    )
    # format source data
    enrich_results <- setNames(
      data.frame(
        assay(enrich_results[[1]]),
        row.names = rowData(enrich_results[[1]])[, 1]
      ),
      colData(enrich_results[[1]])[, 1]
    )
    enrich_categories <- assay(enrich_categories) %>%
      as_tibble()
  }

  # remove uninformative categories and readjust P vals
  source_data <- enrich_results %>%
    filter(!Tissue.Specific.Genes == 0) %>%
    mutate(
      PValue = 10^-Log10PValue,
      fdr = p.adjust(PValue, method = "fdr")
    ) %>%
    rownames_to_column(categoryIdentifier) %>%
    arrange(PValue) %>%
    as_tibble()

  return(list(source_data = source_data, categories = enrich_categories))
}

# RFU plot ----------------------------------------------------------------
plotDxAbundance <- function(
  aptID,
  tidyCounts,
  residualise = FALSE,
  adjustedP = FALSE
) {
  require(ggtext)
  # correct for DE covariates for plotting
  if (residualise) {
    cm <- tidyCounts %>%
      filter(apt_id == aptID) %>%
      pivot_wider(
        id_cols = apt_id,
        names_from = study_id,
        values_from = log2_RFU
      ) %>%
      column_to_rownames("apt_id") %>%
      as.matrix()
    # residualise
    stopifnot(identical(colnames(cm), meta$study_id))
    cm_resid <- resid(lm(
      t(cm) ~ poly(age_at_dx, 2) +
        sex +
        poly(iqr_log2_RFU, 2) +
        poly(gc_days, 2),
      data = meta
    ))

    filt_counts <- data.frame(log2_RFU = cm_resid) %>%
      rownames_to_column("study_id")
  } else {
    filt_counts <- tidyCounts %>%
      filter(apt_id == aptID)
  }
  # filter results for summary stats on plot
  filt_res <- res.limma %>% filter(apt_id == aptID)
  # mainly plot formatting
  filt_counts <- filt_counts %>%
    left_join(meta[, c("study_id", "ref_dx")], by = "study_id") %>%
    mutate(
      ref_dx = case_when(ref_dx == "not_cGCA" ~ "non-cGCA", .default = ref_dx),
      ref_dx = factor(ref_dx, levels = c("non-cGCA", "cGCA"))
    )
  # plot
  p <- filt_counts %>%
    ggplot(aes(x = ref_dx, y = log2_RFU)) +
    geom_violin(aes(fill = ref_dx)) +
    geom_boxplot(width = 0.25, outlier.shape = NA) +
    geom_richtext(
      data = filt_res,
      aes(
        x = 1.5,
        y = Inf,
        vjust = 0.75,
        label = paste0(
          "log<sub>2</sub>(FC) = ",
          round(logFC, 2),
          "<br>",
          if_else(
            adjustedP == TRUE,
            paste0("FDR ", formatPvalLabel(adj.P.Val)),
            formatPvalLabel(P.Value)
          )
        )
      ),
      size = 3,
      fill = NA,
      label.color = NA
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_fill_manual(values = col.dx) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = if_else(
        residualise,
        "log<sub>2</sub> RFU residuals",
        "log<sub>2</sub> RFU"
      ),
      title = paste0(filt_res$Target, " (", filt_res$SeqId, ")")
    ) +
    theme(legend.position = "none")
  p
  return(p)
}


# CLASSIFICATION ----------------------------------------------------------

tidyROC <- function(rocObject) {
  data.frame(
    threshold = rocObject[["thresholds"]],
    specificity = rocObject[["specificities"]],
    sensitivity = rocObject[["sensitivities"]]
  )
}

tidyAUC <- function(rocObject) {
  data.frame(
    conf_lower = ci.auc(rocObject)[[1]],
    auc = ci.auc(rocObject)[[2]],
    conf_higher = ci.auc(rocObject)[[3]]
  )
}


# Adjust abundance for covariates -----------------------------------------

adjustMatrixCovariates <- function(
  mat,
  metadata,
  sampleIdentifier,
  designFormula,
  tidyOutput = FALSE,
  featureIdentifier = NULL,
  valueColName = NULL
) {
  require(limma)
  # check ID orders
  isOrdered <- identical(colnames(mat), metadata[[sampleIdentifier]])
  if (isOrdered == FALSE) {
    message(
      "Sample ID order does not match in matrix and metadata: reordering.."
    )
    mat <- mat[, order(match(colnames(mat), metadata[[sampleIdentifier]]))]
    stopifnot(identical(colnames(mat), metadata[[sampleIdentifier]]))
  } else {
    message("Sample IDs in order\n")
  }
  # design
  design <- model.matrix(designFormula, data = metadata)
  # fit residuals
  resid_mat <- limma::removeBatchEffect(mat, covariates = design)
  #return(resid_mat)

  if (tidyOutput == TRUE) {
    tidy_output <- resid_mat %>%
      as.data.frame() %>%
      rownames_to_column(featureIdentifier) %>%
      pivot_longer(
        cols = !featureIdentifier,
        names_to = sampleIdentifier,
        values_to = valueColName
      ) |>
      drop_na(valueColName)
    return(tidy_output)
  } else {
    return(resid_mat)
  }
}
