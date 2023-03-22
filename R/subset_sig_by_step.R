#' @include DE_functions.R plot.R
NULL

#' @title process data
#'
#' @description filter low expression genes, normalize data by 'TMM' and apply
#'     [limma::voom()], [limma::lmFit()] and [limma::treat()] on normalized data
#'
#' @inheritParams de_analysis
#' @param data expression object
#' @param slot character, specify which slot to use only for DGEList, sce or
#'             seurat object, optional, default 'counts'
#' @param ... params for [voom_fit_treat()]
#'
#' @return A DGEList containing vfit by [limma::voom()] (if normalize = TRUE)
#'         and tfit by [limma::treat()]
#' @export
#'
#' @examples
#' proc_data <- process_data(mastR::im_data_6,
#'   group_col = "celltype:ch1",
#'   target_group = "NK"
#' )
setGeneric(
  "process_data",
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           slot = "counts",
           ...) {
    standardGeneric("process_data")
  }
)

#' @rdname process_data
setMethod(
  "process_data", signature(
    data = "DGEList",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           slot = "counts",
           ...) {
    stopifnot(
      is.logical(normalize),
      is.numeric(filter),
      is.numeric(lfc), is.numeric(p)
    )

    data$original_counts <- data$counts
    data$counts <- data[[slot]]
    ## filter low counts genes
    keep <- filterGenes(
      dge = data,
      group_col = group_col,
      filter = filter,
      normalize = normalize,
      markers = markers,
      gene_id = gene_id
    )
    data <- data[keep, , keep.lib.sizes = FALSE]

    ## normalization for raw counts data
    if (normalize == TRUE) {
      data <- edgeR::calcNormFactors(data, method = "TMM")
    }
    ## update group based group_col and target_group and apply voom
    data <- voom_fit_treat(
      dge = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      lfc = lfc,
      p = p,
      ...
    )
    ## remove vfit if normalize = FALSE as voom isn't applied
    if (normalize == FALSE) data$vfit <- NULL

    return(data)
  }
)

#' @rdname process_data
setMethod(
  "process_data", signature(
    data = "matrix",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           batch = NULL,
           ...) {
    data <- edgeR::DGEList(counts = data, group = group_col)
    if (!is.null(batch)) {
      data$samples <- data.frame(data$samples, batch)
      batch <- colnames(data$samples)[-(1:3)]
    }
    group_col <- "group"

    data <- process_data(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      filter = filter,
      lfc = lfc,
      p = p,
      markers = markers,
      gene_id = gene_id,
      batch = batch,
      slot = "counts",
      ...
    )
    return(data)
  }
)

#' @rdname process_data
setMethod(
  "process_data", signature(
    data = "Matrix",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           batch = NULL,
           ...) {
    data <- edgeR::DGEList(counts = data, group = group_col)
    if (!is.null(batch)) {
      data$samples <- data.frame(data$samples, batch)
      batch <- colnames(data$samples)[-(1:3)]
    }
    group_col <- "group"

    data <- process_data(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      filter = filter,
      lfc = lfc,
      p = p,
      markers = markers,
      gene_id = gene_id,
      batch = batch,
      slot = "counts",
      ...
    )
    return(data)
  }
)

#' @rdname process_data
setMethod(
  "process_data", signature(
    data = "ExpressionSet",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           batch = NULL,
           ...) {
    data <- edgeR::DGEList(
      counts = Biobase::exprs(data),
      samples = Biobase::pData(data),
      group = Biobase::pData(data)[[group_col]]
    )

    if (!is.null(batch)) batch <- make.names(batch)
    group_col <- make.names(group_col)

    data <- process_data(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      filter = filter,
      lfc = lfc,
      p = p,
      markers = markers,
      gene_id = gene_id,
      batch = batch,
      slot = "counts",
      ...
    )
    return(data)
  }
)

#' @rdname process_data
setMethod(
  "process_data", signature(
    data = "SummarizedExperiment",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           slot = "counts",
           batch = NULL,
           ...) {
    data <- edgeR::DGEList(
      counts = SummarizedExperiment::assay(data, slot),
      samples = SummarizedExperiment::colData(data),
      group = SummarizedExperiment::colData(data)[[group_col]]
    )

    if (!is.null(batch)) batch <- make.names(batch)
    group_col <- make.names(group_col)

    data <- process_data(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      filter = filter,
      lfc = lfc,
      p = p,
      markers = markers,
      gene_id = gene_id,
      batch = batch,
      slot = "counts",
      ...
    )
    return(data)
  }
)

#' @rdname process_data
setMethod(
  "process_data", signature(
    data = "Seurat",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           slot = "counts",
           batch = NULL,
           ...) {
    data <- edgeR::DGEList(
      counts = Seurat::GetAssayData(data, slot = slot),
      samples = slot(data, "meta.data"),
      group = slot(data, "meta.data")[[group_col]]
    )

    if (!is.null(batch)) batch <- make.names(batch)
    group_col <- make.names(group_col)

    data <- process_data(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      filter = filter,
      lfc = lfc,
      p = p,
      markers = markers,
      gene_id = gene_id,
      batch = batch,
      slot = "counts",
      ...
    )
    return(data)
  }
)


#' plot diagnostics before and after [process_data()]
#'
#' @param expr1 expression matrix 1 for original data
#' @param expr2 expression matrix 2 for processed data
#' @param group_col vector of group of samples
#' @param abl num, cutoff line
#'
#' @return multiple plots
#' @export
#'
#' @examples
#' dge <- edgeR::DGEList(
#'   counts = Biobase::exprs(im_data_6),
#'   samples = Biobase::pData(im_data_6)
#' )
#' dge$logCPM <- edgeR::cpm(dge, log = TRUE)
#' proc_data <- process_data(dge,
#'   group_col = "celltype.ch1",
#'   target_group = "NK"
#' )
#' plot_diagnostics(proc_data$logCPM, proc_data$vfit$E,
#'   group_col = proc_data$samples$group
#' )
plot_diagnostics <- function(expr1, expr2, group_col, abl = 2) {
  stopifnot(
    "expr1 must be matrix!" = is(expr1, "matrix") | is(expr1, "Matrix"),
    "expr2 must be matrix!" = is(expr2, "matrix") | is(expr2, "Matrix")
  )

  data1 <- tidyr::pivot_longer(as.data.frame(expr1),
    tidyr::everything(),
    names_to = "Sample",
    values_to = "logcounts"
  )
  data1 <- dplyr::left_join(data1, data.frame(
    Sample = colnames(expr1),
    Group = group_col
  ),
  by = "Sample"
  )
  data2 <- tidyr::pivot_longer(as.data.frame(expr2),
    tidyr::everything(),
    names_to = "Sample",
    values_to = "logcounts"
  )
  data2 <- dplyr::left_join(data2, data.frame(
    Sample = colnames(expr2),
    Group = group_col
  ),
  by = "Sample"
  )

  ## density plot for each sample
  p1 <- plot_density_init(data1, data2, abl)
  ## boxplot of RLE
  p2 <- plot_rle_init(expr1, group_col) + ggtitle("Original data") +
    plot_rle_init(expr2, group_col) + ggtitle("Processed data") +
    patchwork::plot_layout(guides = "collect", ncol = 1)
  ## MDS plot
  p3 <- plot_MDS_init(expr1, expr2, group_col)

  return(list("density" = p1, "RLE" = p2, "MDS" = p3))
}

#' plot Mean-variance trend after voom and after final linear fit
#'
#' @param proc_data processed data returned by [process_data()]
#' @param span num, span for [lowess()]
#'
#' @return comparison plot of mean-variance of voom and final model
#' @export
#'
#' @examples
#' proc_data <- process_data(im_data_6,
#'   group_col = "celltype:ch1",
#'   target_group = "NK"
#' )
#' plot_mean_var(proc_data)
plot_mean_var <- function(proc_data, span = 0.5) {
  if (is.null(proc_data$vfit) | is.null(proc_data$tfit)) {
    stop("vfit or tfit is not found in proc_data!")
  }
  p1 <- plot_voom(proc_data$vfit, span = span)
  p2 <- plot_sa(proc_data$tfit, title = "Final model: Mean-variance trend")

  p1 + p2
}

#' select DEGs from multiple comparisons
#'
#' @param tfit processed tfit by [limma::treat()] or
#'             processed data returned by [process_data()]
#' @param feature_selection one of "auto" (default), "rankproduct" or "none",
#'                          choose if to use rank product or not to select DEGs
#'                          from multiple comparisons of DE analysis, 'auto'
#'                          uses 'rankproduct' but change to 'none' if final
#'                          genes < 5 for both UP and DOWN
#' @param ... params for [DEGs_RP()] or [DEGs_Group()]
#'
#' @return GeneSetCollection contains UP and DOWN genesets
#' @export
#'
#' @examples
#' proc_data <- process_data(im_data_6,
#'   group_col = "celltype:ch1",
#'   target_group = "NK"
#' )
#' select_sig(proc_data$tfit)
select_sig <- function(tfit,
                       feature_selection = c("auto", "rankproduct", "none"),
                       ...) {
  if (is(tfit, "DGEList")) {
    stopifnot("No 'tfit' is found!" = "tfit" %in% names(tfit))
    tfit <- tfit$tfit
  }

  stopifnot("tfit is not MArrayLM!" = is(tfit, "MArrayLM"))
  feature_selection <- match.arg(feature_selection)

  ## assemble DEGs from comparisons by Rank Product or simply intersect/union
  if (feature_selection == "auto") {
    DEGs <- DEGs_RP(tfit = tfit, ...)
    if (all(lengths(DEGs) < 5)) {
      DEGs <- DEGs_Group(tfit = tfit, ...)
    }
  } else if (feature_selection == "rankproduct") {
    DEGs <- DEGs_RP(tfit = tfit, ...)
  } else {
    DEGs <- DEGs_Group(tfit = tfit, ...)
  }

  DEGs <- gls2gsc(DEGs) ## convert to GeneSetCollection object
  return(DEGs)
}
