#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "matrix",
    sigs = "vector",
    group_col = "vector"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- rankdensity_init(
      expr = data, sigs = sigs, by = group_col,
      aggregate = aggregate, gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "Matrix",
    sigs = "vector",
    group_col = "vector"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- rankdensity_init(
      expr = data, sigs = sigs, by = group_col,
      aggregate = aggregate, gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "data.frame",
    sigs = "vector",
    group_col = "vector"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- rankdensity_init(
      expr = as.matrix(data), sigs = sigs, by = group_col,
      aggregate = aggregate, gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "DGEList",
    sigs = "vector",
    group_col = "character"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           slot = "counts",
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- sig_rankdensity_plot(
      data = data[[slot]], sigs = sigs,
      group_col = data$samples[[group_col]],
      aggregate = aggregate, gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "ExpressionSet",
    sigs = "vector",
    group_col = "character"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- sig_rankdensity_plot(
      data = Biobase::exprs(data), sigs = sigs,
      group_col = data[[group_col]],
      aggregate = aggregate, gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "Seurat",
    sigs = "vector",
    group_col = "character"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           slot = "counts",
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- sig_rankdensity_plot(
      data = SeuratObject::GetAssayData(data, slot = slot),
      sigs = sigs, group_col = slot(data, "meta.data")[[group_col]],
      aggregate = aggregate,
      gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "SummarizedExperiment",
    sigs = "vector",
    group_col = "character"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           slot = "counts",
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    p <- sig_rankdensity_plot(
      data = SummarizedExperiment::assay(data, slot),
      sigs = sigs, group_col = SummarizedExperiment::colData(data)[[group_col]],
      aggregate = aggregate, gene_id = gene_id
    )

    return(p)
  }
)

#' @rdname sig_rankdensity_plot
setMethod(
  "sig_rankdensity_plot", signature(
    data = "list",
    sigs = "vector",
    group_col = "character"
  ),
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           slot = "counts",
           gene_id = "SYMBOL") {
    stopifnot(is.character(gene_id), is.logical(aggregate))

    if (length(group_col) == 1) {
      group_col <- rep(group_col, length(data))
    }
    if (length(slot) == 1) {
      slot <- rep(slot, length(data))
    }
    if (length(aggregate) == 1) {
      aggregate <- rep(aggregate, length(data))
    }
    if (length(gene_id) == 1) {
      gene_id <- rep(gene_id, length(data))
    }

    p <- list()
    for (i in seq_along(data)) {
      p[[i]] <- sig_rankdensity_plot(
        data = data[[i]], sigs = sigs,
        group_col = group_col[[i]],
        slot = slot[i],
        aggregate = aggregate[i],
        gene_id = gene_id[i]
      )

      if (!is.null(names(data))) {
        p[[i]] <- p[[i]] +
          patchwork::plot_annotation(
            subtitle = names(data)[i],
            theme = theme(plot.subtitle = element_text(hjust = 0.5))
          )
      }
      p[[i]] <- ggpubr::as_ggplot(patchwork::patchworkGrob(p[[i]]))
    }

    p <- patchwork::wrap_plots(p) +
      patchwork::plot_layout(guides = "collect")
    return(p)
  }
)
