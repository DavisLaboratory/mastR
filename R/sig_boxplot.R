#' @include plot.R
#' @import ggplot2
NULL

#' Boxplot of median expression or scores of signature
#'
#' Make boxplot and show expression or score level of signature across subsets.
#'
#' @param data expression data, can be matrix, DGEList, eSet, seurat, sce...
#' @param sigs a vector of signature (Symbols)
#' @param ID character or vector, specify the column name to compare in coldata
#' @param type pattern, specify the group of interest as reference
#' @param plot.score logical, if to plot score instead of expression of signature
#' @param counts logical, indicate if data is raw counts data
#' @param method a character string indicating which method to be used for
#'               `stat_compare_means()` to compare the means across types,
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
#'   ID = "celltype:ch1", type = "NK",
#'   gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric("sig_boxplot",
           function(data,
                    sigs,
                    ID,
                    type,
                    plot.score = TRUE,
                    counts = TRUE,
                    method = "t.test",
                    slot = "counts",
                    gene_id = "SYMBOL")
             standardGeneric("sig_boxplot"))

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "matrix",
  sigs = "vector",
  ID = "vector",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  stopifnot(is.logical(counts), is.character(method), is.character(gene_id))

  if(plot.score == TRUE) {
    ## plot scores of signature
    ## calculate scores
    scores <- singscore_init(expr = data, sigs = sigs, by = ID,
                             counts = counts, gene_id = gene_id)
    ## boxplot of scores
    p <- score_boxplot_init(scores = scores, by = ID,
                            type = type, method = method)
  } else {
    ## plot median expression of signature
    p <- exp_boxplot_init(expr = data, sigs = sigs, type = type,
                          by = ID, counts = counts,
                          method = method,
                          gene_id = gene_id)
  }

  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "Matrix",
  sigs = "vector",
  ID = "vector",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  stopifnot(is.logical(counts), is.character(method), is.character(gene_id))

  if(plot.score == TRUE) {
    ## plot scores of signature
    ## calculate scores
    scores <- singscore_init(expr = as.matrix(data), sigs = sigs, by = ID,
                             counts = counts, gene_id = gene_id)
    ## boxplot of scores
    p <- score_boxplot_init(scores = scores, by = ID,
                            type = type, method = method)
  } else {
    ## plot median expression of signature
    p <- exp_boxplot_init(expr = data, sigs = sigs, type = type,
                          by = ID, counts = counts,
                          method = method,
                          gene_id = gene_id)
  }

  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "data.frame",
  sigs = "vector",
  ID = "vector",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  stopifnot(is.logical(counts), is.character(method), is.character(gene_id))

  p <- sig_boxplot(data = as.matrix(data),
                   sigs = sigs,
                   ID = ID,
                   type = type,
                   plot.score = plot.score,
                   counts = counts,
                   method = method,
                   gene_id = gene_id)

  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "DGEList",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = data$counts,
                   sigs = sigs,
                   ID = data$samples[[ID]],
                   type = type,
                   plot.score = plot.score,
                   counts = counts,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "ExpressionSet",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = Biobase::exprs(data),
                   sigs = sigs,
                   ID = Biobase::pData(data)[[ID]],
                   type = type,
                   plot.score = plot.score,
                   counts = counts,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "Seurat",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         slot = "counts",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = Seurat::GetAssayData(data, slot = slot),
                   sigs = sigs,
                   ID = data@meta.data[[ID]],
                   type = type,
                   plot.score = plot.score,
                   counts = counts,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "SummarizedExperiment",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         slot = "counts",
         gene_id = "SYMBOL") {

  p <- sig_boxplot(data = SummarizedExperiment::assay(data, slot),
                   sigs = sigs,
                   ID = data@colData[[ID]],
                   type = type,
                   plot.score = plot.score,
                   counts = counts,
                   method = method,
                   gene_id = gene_id)
  return(p)
})

#' @rdname sig_boxplot
setMethod("sig_boxplot", signature(
  data = "list",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         plot.score = TRUE,
         counts = TRUE,
         method = "t.test",
         slot = "counts",
         gene_id = "SYMBOL") {

  if(length(counts) == 1)
    counts <- rep(counts, length(data))
  if(length(type) == 1)
    type <- rep(type, length(data))
  if(length(ID) == 1)
    ID <- rep(ID, length(data))
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
                          ID = ID[[i]],
                          type = type[i],
                          plot.score = plot.score[i],
                          counts = counts[i],
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
