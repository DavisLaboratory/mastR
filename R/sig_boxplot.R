#' @include plot.R
#' @import ggplot2
NULL

#' Boxplot of median expression or scores of signature
#'
#' Make boxplot and show expression or score level of signature across subsets.
#'
#' @param data expression data, can be matrix, DGEList, eSet, seurat, sce...
#' @param sigs a vector of signature (Symbols)
#' @param group_col character or vector, specify the column name to compare in coldata
#' @param target_group pattern, specify the group of interest as reference
#' @param plot.score logical, if to plot score instead of expression of signature
#' @param normalize logical, TRUE indicates raw counts data to normalize and
#'                  calculate cpm
#' @param method a character string indicating which method to be used for
#'               `stat_compare_means()` to compare the means across groups,
#'               could be "t.test", 'wilcox.test', 'anova'..., default "t.test"
#' @param slot character, indicate which slot used for seurat/sce object
#' @param gene_id character, indicate the ID type of rowname of expression data's ,
#'                could be one of 'ENSEMBL', 'SYMBOL', ... default 'SYMBOL'
#'
#' @return patchwork or ggplot of boxplot
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' p <- sig_boxplot(
#'   im_data_6, sigs = NK_markers$HGNC_Symbol[1:30],
#'   group_col = "celltype:ch1", target_group = "NK",
#'   gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric("sig_boxplot",
           function(data,
                    sigs,
                    group_col,
                    target_group,
                    plot.score = TRUE,
                    normalize = TRUE,
                    method = "t.test",
                    slot = "counts",
                    gene_id = "SYMBOL")
             standardGeneric("sig_boxplot"))

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "matrix",
  sigs = "vector",
  group_col = "vector",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  stopifnot(is.logical(normalize), is.character(method), is.character(gene_id))

  if(plot.score == TRUE) {
    ## plot scores of signature
    ## calculate scores
    scores <- singscore_init(expr = data, sigs = sigs, by = group_col,
                             normalize = normalize, gene_id = gene_id)
    ## boxplot of scores
    p <- score_boxplot_init(scores = scores, by = group_col,
                            target_group = target_group, method = method)
  } else {
    ## plot median expression of signature
    p <- exp_boxplot_init(expr = data, sigs = sigs,
                          target_group = target_group,
                          by = group_col, normalize = normalize,
                          method = method,
                          gene_id = gene_id)
  }

  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "Matrix",
  sigs = "vector",
  group_col = "vector",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  stopifnot(is.logical(normalize), is.character(method), is.character(gene_id))

  if(plot.score == TRUE) {
    ## plot scores of signature
    ## calculate scores
    scores <- singscore_init(expr = as.matrix(data),
                             sigs = sigs,
                             by = group_col,
                             normalize = normalize,
                             gene_id = gene_id)
    ## boxplot of scores
    p <- score_boxplot_init(scores = scores, by = group_col,
                            target_group = target_group, method = method)
  } else {
    ## plot median expression of signature
    p <- exp_boxplot_init(expr = data, sigs = sigs,
                          target_group = target_group,
                          by = group_col, normalize = normalize,
                          method = method,
                          gene_id = gene_id)
  }

  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "data.frame",
  sigs = "vector",
  group_col = "vector",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  stopifnot(is.logical(normalize), is.character(method), is.character(gene_id))

  p <- sig_boxplot(data = as.matrix(data),
                   sigs = sigs,
                   group_col = group_col,
                   target_group = target_group,
                   plot.score = plot.score,
                   normalize = normalize,
                   method = method,
                   gene_id = gene_id)

  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "DGEList",
  sigs = "vector",
  group_col = "character",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = data$counts,
                   sigs = sigs,
                   group_col = data$samples[[group_col]],
                   target_group = target_group,
                   plot.score = plot.score,
                   normalize = normalize,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "ExpressionSet",
  sigs = "vector",
  group_col = "character",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = Biobase::exprs(data),
                   sigs = sigs,
                   group_col = Biobase::pData(data)[[group_col]],
                   target_group = target_group,
                   plot.score = plot.score,
                   normalize = normalize,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "Seurat",
  sigs = "vector",
  group_col = "character",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         slot = "counts",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = Seurat::GetAssayData(data, slot = slot),
                   sigs = sigs,
                   group_col = data@meta.data[[group_col]],
                   target_group = target_group,
                   plot.score = plot.score,
                   normalize = normalize,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "SummarizedExperiment",
  sigs = "vector",
  group_col = "character",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         slot = "counts",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = SummarizedExperiment::assay(data, slot),
                   sigs = sigs,
                   group_col = data@colData[[group_col]],
                   target_group = target_group,
                   plot.score = plot.score,
                   normalize = normalize,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "list",
  sigs = "vector",
  group_col = "character",
  target_group = "character"
),
function(data,
         sigs,
         group_col,
         target_group,
         plot.score = TRUE,
         normalize = TRUE,
         method = "t.test",
         slot = "counts",
         gene_id = "SYMBOL") {

  if(length(normalize) == 1)
    normalize <- rep(normalize, length(data))
  if(length(target_group) == 1)
    target_group <- rep(target_group, length(data))
  if(length(group_col) == 1)
    group_col <- rep(group_col, length(data))
  if(length(slot) == 1)
    slot <- rep(slot, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))
  if(length(plot.score) == 1)
    plot.score <- rep(plot.score, length(data))

  p <- list()
  for (i in seq_along(data)) {
    p[[i]] <- sig_boxplot(data = data[[i]],
                          sigs = sigs,
                          group_col = group_col[[i]],
                          target_group = target_group[i],
                          plot.score = plot.score[i],
                          normalize = normalize[i],
                          method = method,
                          slot = slot[i],
                          gene_id = gene_id[i])
    if(!is.null(names(data)))
      p[[i]] <- p[[i]] + ggplot2::labs(subtitle = names(data)[i]) +
        ggplot2::theme(plot.subtitle = element_text(hjust = 0.5))
  }
  p <- patchwork::wrap_plots(p) +
    patchwork::plot_layout(guides = "collect")
  return(p)
})
