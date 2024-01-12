#' @include plot.R
NULL

#' @rdname sig_gseaplot
setMethod(
  "sig_gseaplot", signature(
    data = "MArrayLM",
    sigs = "vector"
  ),
  function(data,
           sigs,
           gene_id = "SYMBOL",
           method = c("dotplot", "gseaplot"),
           col = "-log10(p.adjust)",
           size = "enrichmentScore",
           pvalue_table = FALSE,
           digits = 2,
           rank_stat = "logFC",
           ...) {
    stopifnot(is.character(gene_id), is.numeric(digits))
    method <- match.arg(method)

    ## get DEGs tables list with statistics
    ## save DE result tables into list
    tDEG <- list()
    for (i in seq_len(ncol(data))) {
      ## use limma::topTreat() to get statistics of DEA
      tDEG[[i]] <- na.omit(limma::topTreat(data, coef = i, number = Inf))
    }
    names(tDEG) <- colnames(data)

    gsets <- data.frame(SYMBOL = sigs, set = "Signature")
    ids <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
      gsets[, 1],
      columns = gene_id,
      keytype = "SYMBOL"
    )
    gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

    ## gsea
    gse <- gsea_analysis(
      tDEG = tDEG, gsets = gsets,
      gene_id = gene_id, digits = digits,
      rank_stat = rank_stat
    )

    if (method == "gseaplot") {
      p <- gsea_plot_init(gse, pvalue_table = pvalue_table)
    } else {
      p <- gsea_dotplot_init(gse, col = col, size = size)
    }

    p <- patchwork::wrap_plots(p) + patchwork::plot_layout(guides = "collect")
    return(p)
  }
)

#' @rdname sig_gseaplot
setMethod(
  "sig_gseaplot", signature(
    data = "MArrayLM",
    sigs = "list"
  ),
  function(data,
           sigs,
           gene_id = "SYMBOL",
           method = c("dotplot", "gseaplot"),
           col = "-log10(p.adjust)",
           size = "enrichmentScore",
           pvalue_table = FALSE,
           digits = 2,
           rank_stat = "logFC",
           ...) {
    stopifnot(is.character(gene_id), is.numeric(digits))
    method <- match.arg(method)

    ## get DEGs tables list with statistics
    ## save DE result tables into list
    tDEG <- list()
    for (i in seq_len(ncol(data))) {
      ## use limma::topTreat() to get statistics of DEA
      tDEG[[i]] <- na.omit(limma::topTreat(data, coef = i, number = Inf))
    }
    names(tDEG) <- colnames(data)

    if (is.null(names(sigs))) {
      names(sigs) <- seq_along(sigs)
    } ## set gene list names
    gsets <- utils::stack(sigs)
    colnames(gsets) <- c("SYMBOL", "set")
    ids <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
      gsets[, 1],
      columns = gene_id,
      keytype = "SYMBOL"
    )
    gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

    ## gsea
    gse <- gsea_analysis(
      tDEG = tDEG, gsets = gsets,
      gene_id = gene_id, digits = digits,
      rank_stat = rank_stat
    )

    if (method == "gseaplot") {
      p <- gsea_plot_init(gse, pvalue_table = pvalue_table)
    } else {
      p <- gsea_dotplot_init(gse, col = col, size = size)
    }

    p <- patchwork::wrap_plots(p) + patchwork::plot_layout(guides = "collect")
    return(p)
  }
)

#' @rdname sig_gseaplot
setMethod(
  "sig_gseaplot", signature(
    data = "DGEList",
    sigs = "ANY",
    group_col = "ANY",
    target_group = "ANY"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           gene_id = "SYMBOL",
           slot = "counts",
           method = c("dotplot", "gseaplot"),
           col = "-log10(p.adjust)",
           size = "enrichmentScore",
           pvalue_table = FALSE,
           digits = 2,
           rank_stat = "logFC",
           ...) {
    stopifnot(is.character(gene_id), is.numeric(digits))
    method <- match.arg(method)

    ## get processed data if there's no tfit in it
    if (is.null(data$tfit)) {
      data <- process_data(
        data = data, group_col = group_col,
        target_group = target_group,
        markers = Reduce(union, sigs),
        gene_id = gene_id,
        slot = slot,
        ...
      )
    }

    p <- sig_gseaplot(
      data = data$tfit, sigs = sigs,
      gene_id = gene_id,
      method = method,
      col = col,
      size = size,
      pvalue_table = pvalue_table,
      digits = digits,
      rank_stat = rank_stat,
      ...
    )
    return(p)
  }
)

#' @rdname sig_gseaplot
setMethod(
  "sig_gseaplot", signature(
    data = "ANY",
    sigs = "ANY",
    group_col = "ANY",
    target_group = "ANY"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           gene_id = "SYMBOL",
           slot = "counts",
           method = c("dotplot", "gseaplot"),
           col = "-log10(p.adjust)",
           size = "enrichmentScore",
           pvalue_table = FALSE,
           digits = 2,
           rank_stat = "logFC",
           ...) {
    stopifnot(is.character(gene_id), is.numeric(digits))
    method <- match.arg(method)

    ## get processed data
    data <- process_data(
      data = data, group_col = group_col,
      target_group = target_group,
      markers = Reduce(union, sigs),
      gene_id = gene_id,
      slot = slot,
      ...
    )

    p <- sig_gseaplot(
      data = data$tfit, sigs = sigs,
      gene_id = gene_id,
      method = method,
      col = col,
      size = size,
      pvalue_table = pvalue_table,
      digits = digits,
      rank_stat = rank_stat,
      ...
    )
    return(p)
  }
)

#' @rdname sig_gseaplot
setMethod(
  "sig_gseaplot", signature(
    data = "list",
    sigs = "ANY",
    group_col = "ANY",
    target_group = "ANY"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           gene_id = "SYMBOL",
           slot = "counts",
           method = c("dotplot", "gseaplot"),
           col = "-log10(p.adjust)",
           size = "enrichmentScore",
           pvalue_table = FALSE,
           digits = 2,
           rank_stat = "logFC",
           ...) {
    stopifnot(is.character(gene_id), is.numeric(digits))
    method <- match.arg(method)

    if (length(target_group) == 1) {
      target_group <- rep(target_group, length(data))
    }
    if (length(group_col) == 1) {
      group_col <- rep(group_col, length(data))
    }
    if (length(gene_id) == 1) {
      gene_id <- rep(gene_id, length(data))
    }
    if (length(slot) == 1) {
      slot <- rep(slot, length(data))
    }

    p <- list()
    for (i in seq_along(data)) {
      p[[i]] <- sig_gseaplot(
        data = data[[i]],
        sigs = sigs,
        group_col = group_col[[i]],
        target_group = target_group[i],
        gene_id = gene_id[i],
        slot = slot[i],
        method = method,
        col = col,
        size = size,
        pvalue_table = pvalue_table,
        digits = digits,
        rank_stat = rank_stat,
        ...
      )
      p[[i]] <- ggpubr::as_ggplot(patchwork::patchworkGrob(p[[i]]))
      if (!is.null(names(data))) {
        p[[i]] <- p[[i]] + ggplot2::labs(subtitle = names(data)[i]) +
          ggplot2::theme(plot.subtitle = element_text(hjust = 0.5))
      }
    }
    p <- patchwork::wrap_plots(p) +
      patchwork::plot_layout(guides = "collect")

    return(p)
  }
)
