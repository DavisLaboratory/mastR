#' @include DE_functions.R plot.R
NULL

#' @rdname get_degs
setMethod(
  "get_degs", signature(
    data = "DGEList",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           slot = "counts",
           batch = NULL,
           ...) {
    feature_selection <- match.arg(feature_selection)
    stopifnot(
      "Please provide column names as batch!" =
        is.null(batch) | is.vector(batch),
      "slot must be character!" = is.character(slot)
    )

    ## process data, filtration, normalization, fit
    proc_data <- process_data(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      slot = slot,
      batch = batch,
      ...
    )

    ## plot diagnostics before and after process_data()
    if (!is.null(list(...)$plot) && list(...)$plot == TRUE) {
      ## generate expr1 adn expr2 to compare before and after process
      if (slot == "counts") {
        expr1 <- proc_data$original_counts
      } else {
        expr1 <- proc_data[[slot]]
      }
      if (normalize == TRUE) {
        expr1 <- edgeR::cpm(expr1, log = TRUE)
        expr2 <- proc_data$vfit$E
      } else {
        expr2 <- proc_data$counts
      }

      plot_diagnostics(expr1, expr2,
        group_col = proc_data$samples[[group_col]]
      )
      if (normalize == TRUE) {
        plot_mean_var(proc_data)
      } else {
        limma::plotSA(proc_data$tfit,
          main = "Final model: Mean-variance trend"
        )
      }
    }

    ## select DEGs from multiple comparisons
    DEGs <- select_sig(
      tfit = proc_data$tfit,
      feature_selection = feature_selection,
      ...
    )

    DEGs <- list(DEGs = DEGs, proc_data = proc_data)

    return(DEGs)
  }
)

#' @rdname get_degs
setMethod(
  "get_degs", signature(
    data = "matrix",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           batch = NULL,
           ...) {
    DGE <- edgeR::DGEList(counts = data, group = group_col)
    if (!is.null(batch)) {
      DGE$samples <- data.frame(DGE$samples, batch)
      batch <- colnames(DGE$samples)[-seq_len(3)]
    }
    group_col <- "group"
    rm(data)

    DEGs <- get_degs(
      data = DGE, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      batch = batch, slot = "counts",
      ...
    )
    return(DEGs)
  }
)

#' @rdname get_degs
setMethod(
  "get_degs", signature(
    data = "Matrix",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           batch = NULL,
           ...) {
    DGE <- edgeR::DGEList(counts = data, group = group_col)
    if (!is.null(batch)) {
      DGE$samples <- data.frame(DGE$samples, batch)
      batch <- colnames(DGE$samples)[-seq_len(3)]
    }
    group_col <- "group"
    rm(data)

    DEGs <- get_degs(
      data = DGE, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      batch = batch, slot = "counts",
      ...
    )
    return(DEGs)
  }
)

#' @rdname get_degs
setMethod(
  "get_degs", signature(
    data = "ExpressionSet",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           batch = NULL,
           ...) {
    expr <- Biobase::exprs(data)
    coldata <- Biobase::pData(data)

    DGE <- edgeR::DGEList(
      counts = expr,
      samples = coldata,
      group = coldata[[group_col]]
    )
    if (!is.null(batch)) batch <- make.names(batch)
    group_col <- make.names(group_col)
    rm(data, expr, coldata)

    DEGs <- get_degs(
      data = DGE, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      batch = batch, slot = "counts",
      ...
    )
    return(DEGs)
  }
)

#' @rdname get_degs
setMethod(
  "get_degs", signature(
    data = "SummarizedExperiment",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           slot = "counts",
           batch = NULL,
           ...) {
    expr <- SummarizedExperiment::assay(data, slot)
    coldata <- SummarizedExperiment::colData(data)

    DGE <- edgeR::DGEList(
      counts = expr,
      samples = coldata,
      group = coldata[[group_col]]
    )
    if (!is.null(batch)) batch <- make.names(batch)
    group_col <- make.names(group_col)
    rm(data, expr, coldata)

    DEGs <- get_degs(
      data = DGE, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      batch = batch, slot = "counts",
      ...
    )
    return(DEGs)
  }
)

#' @rdname get_degs
setMethod(
  "get_degs", signature(
    data = "Seurat",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           slot = "counts",
           batch = NULL,
           ...) {
    expr <- SeuratObject::GetAssayData(data, slot = slot)
    coldata <- slot(data, "meta.data")

    DGE <- edgeR::DGEList(
      counts = expr,
      samples = coldata,
      group = coldata[[group_col]]
    )
    if (!is.null(batch)) batch <- make.names(batch)
    group_col <- make.names(group_col)
    rm(data, expr, coldata)

    DEGs <- get_degs(
      data = DGE, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      batch = batch, slot = "counts",
      ...
    )
    return(DEGs)
  }
)


#' @title DE analysis pipeline
#'
#' @description Standard DE analysis by using edgeR and limma::voom pipeline
#'
#' @param dge DGEList object for DE analysis, including expr and samples info
#' @param group_col character, column name of coldata to specify the DE comparisons
#' @param target_group pattern, specify the group of interest, e.g. NK
#' @param normalize logical, if the expr in data is raw counts needs to be normalized
#' @param group logical, TRUE to separate samples into only 2 groups:
#'              `target_group`` and 'Others'; FALSE to set each level as a group
#' @param filter a vector of 2 numbers, filter condition to remove low expression
#'               genes, the 1st for min.counts (if normalize = TRUE) or CPM/TPM
#'               (if normalize = FALSE), the 2nd for samples size 'large.n'
#' @param plot logical, if to make plots to show QC before and after filtration
#' @param lfc num, cutoff of logFC for DE analysis
#' @param p num, cutoff of p value for DE analysis and permutation test if
#'          feature_selection = "rankproduct"
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                kept anyway after filtration. Default 'NULL' means no special
#'                genes need to be kept.
#' @param gene_id character, specify the gene ID target_group of rownames of expression data
#'                when markers is not NULL, could be one of 'ENSEMBL', 'SYMBOL',
#'                'ENTREZ'..., default 'SYMBOL'
#' @param slot character, specify which slot to use for DGEList,
#'             default 'counts'
#' @param batch vector of character, column name(s) of coldata to be treated as
#'              batch effect factor, default NULL
#' @param summary logical, if to show the summary of DE analysis
#' @param ... omitted
#'
#' @return MArrayLM object generated by [limma::treat()]
#'
#' @examples
#' data("im_data_6")
#' dge <- edgeR::DGEList(
#'   counts = Biobase::exprs(im_data_6),
#'   samples = Biobase::pData(im_data_6)
#' )
#' de_analysis(dge, group_col = "celltype.ch1", target_group = "NK")
#'
#' @export
de_analysis <- function(dge,
                        group_col,
                        target_group,
                        normalize = TRUE,
                        group = FALSE,
                        filter = c(10, 10),
                        plot = FALSE,
                        lfc = 0,
                        p = 0.05,
                        markers = NULL,
                        gene_id = "SYMBOL",
                        slot = "counts",
                        batch = NULL,
                        summary = TRUE,
                        ...) {
  stopifnot(
    is.logical(normalize),
    is.logical(group), is.numeric(filter),
    is.logical(plot), is.numeric(lfc),
    is.numeric(p), is.character(gene_id)
  )

  stopifnot(
    "Please provide column names as batch!" =
      is.null(batch) | is.vector(batch),
    "slot must be character!" = is.character(slot)
  )

  ## process data, filtration, normalization, fit
  proc_data <- process_data(
    data = dge,
    group_col = group_col,
    target_group = target_group,
    normalize = normalize,
    group = group,
    filter = filter,
    lfc = lfc,
    p = p,
    markers = markers,
    gene_id = gene_id,
    slot = slot,
    batch = batch,
    summary = summary,
    ...
  )

  ## plot diagnostics before and after process_data()
  if (plot == TRUE) {
    ## generate expr1 adn expr2 to compare before and after process
    if (slot == "counts") {
      expr1 <- proc_data$original_counts
    } else {
      expr1 <- proc_data[[slot]]
    }
    if (normalize == TRUE) {
      expr1 <- edgeR::cpm(expr1, log = TRUE)
      expr2 <- proc_data$vfit$E
    } else {
      expr2 <- proc_data$counts
    }

    plot_diagnostics(expr1, expr2,
      group_col = proc_data$samples[[group_col]]
    )
    if (normalize == TRUE) {
      plot_mean_var(proc_data)
    } else {
      limma::plotSA(proc_data$tfit,
        main = "Final model: Mean-variance trend"
      )
    }
  }

  return(proc_data)
}
