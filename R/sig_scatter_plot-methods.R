#' @include plot.R
#' @import ggplot2
NULL

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "matrix",
    sigs = "vector",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    ## scale data by row
    data <- t(scale(t(data)))

    p <- scatter_plot_init(
      expr = data, sigs = sigs,
      target_group = target_group,
      by = group_col, xint = xint, yint = yint,
      gene_id = gene_id
    )
    return(p)
  }
)

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "Matrix",
    sigs = "vector",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    p <- sig_scatter_plot(
      data = as.matrix(data), sigs = sigs,
      group_col = group_col,
      target_group = target_group,
      xint = xint, yint = yint,
      gene_id = gene_id
    )
    return(p)
  }
)

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "DGEList",
    sigs = "vector",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           slot = "counts",
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    p <- sig_scatter_plot(
      data = data[[slot]], sigs = sigs,
      group_col = data$samples[[group_col]],
      target_group = target_group,
      xint = xint, yint = yint,
      gene_id = gene_id
    )
    return(p)
  }
)

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "ExpressionSet",
    sigs = "vector",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    p <- sig_scatter_plot(
      data = Biobase::exprs(data),
      sigs = sigs,
      group_col = data[[group_col]],
      target_group = target_group,
      xint = xint, yint = yint,
      gene_id = gene_id
    )
    return(p)
  }
)

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "Seurat",
    sigs = "vector",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           slot = "counts",
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    p <- sig_scatter_plot(
      data = SeuratObject::GetAssayData(data, slot = slot),
      sigs = sigs,
      group_col = slot(data, "meta.data")[[group_col]],
      target_group = target_group,
      xint = xint, yint = yint,
      gene_id = gene_id
    )
    return(p)
  }
)

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "SummarizedExperiment",
    sigs = "vector",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           slot = "counts",
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    p <- sig_scatter_plot(
      data = SummarizedExperiment::assay(data, slot),
      sigs = sigs,
      group_col = SummarizedExperiment::colData(data)[[group_col]],
      target_group = target_group,
      xint = xint, yint = yint,
      gene_id = gene_id
    )
    return(p)
  }
)

#' @rdname sig_scatter_plot
setMethod(
  "sig_scatter_plot", signature(
    data = "list",
    sigs = "vector",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           sigs,
           group_col,
           target_group,
           slot = "counts",
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    if (length(target_group) == 1) {
      target_group <- rep(target_group, length(data))
    }
    if (length(group_col) == 1) {
      group_col <- rep(group_col, length(data))
    }
    if (length(slot) == 1) {
      slot <- rep(slot, length(data))
    }
    if (length(gene_id) == 1) {
      gene_id <- rep(gene_id, length(data))
    }

    p <- list()
    for (i in seq_along(data)) {
      p[[i]] <- sig_scatter_plot(
        data = data[[i]], sigs = sigs,
        target_group = target_group[i],
        slot = slot[i], group_col = group_col[[i]],
        xint = xint, yint = yint,
        gene_id = gene_id[i]
      )
      if (!is.null(names(data))) {
        p[[i]] <- p[[i]] + labs(subtitle = names(data)[i]) +
          theme(plot.subtitle = element_text(hjust = 0.5))
      }
    }
    p <- patchwork::wrap_plots(p) +
      patchwork::plot_layout(guides = "collect")
    return(p)
  }
)
