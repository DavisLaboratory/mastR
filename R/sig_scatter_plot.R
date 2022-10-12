#' @include plot.R
#' @import ggplot2
NULL

#' @title Biplot of signature for specific subset vs others
#'
#' @description  Biplot depicts mean expression for each signature gene in
#'   specific subset against other cell types.
#'
#' @param data expression data objects, can be DGEList, eSet, Seurat, sce...
#' @param sigs a vector of signature (Symbols)
#' @param ID a vector of chr, specify the column name to compare in coldata
#' @param type pattern, specify the group of interest as reference
#' @param counts logical, indicate if data is raw counts data
#' @param slot chr, indicate which slot used for seurat/sce object
#' @param xint intercept of vertical dashed line
#' @param yint intercept of horizontal dashed line
#' @param gene_id chr, indicate the ID type of rowname of expression data's ,
#'                could be one of 'ENSEMBL', 'SYMBOL', ... default 'SYMBOL'
#'
#' @return patchwork or ggplot of scatter plot of median expression
#'
#' @examples
#' sig_scatter_plot(
#'   sigs = NK_markers$HGNC_Symbol, data = im_data_6,
#'   ID = "celltype:ch1", type = "NK",
#'   gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric("sig_scatter_plot",
           function(data,
                    sigs,
                    ID,
                    type,
                    counts = TRUE,
                    slot = "counts",
                    xint = 1,
                    yint = 1,
                    gene_id = "SYMBOL")
             standardGeneric("sig_scatter_plot"))

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "matrix",
  sigs = "vector",
  ID = "vector",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         xint = 1,
         yint = 1,
         gene_id = "SYMBOL") {

  p <- scatter_plot_init(expr = data, sigs = sigs, counts = counts,
                         type = type, by = ID, xint = xint, yint = yint,
                         gene_id = gene_id)
  return(p)
})

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "Matrix",
  sigs = "vector",
  ID = "vector",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         xint = 1,
         yint = 1,
         gene_id = "SYMBOL") {

  p <- scatter_plot_init(expr = data, sigs = sigs, counts = counts,
                         type = type, by = ID, xint = xint, yint = yint,
                         gene_id = gene_id)
  return(p)
})

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "DGEList",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         xint = 1,
         yint = 1,
         gene_id = "SYMBOL") {

  p <- scatter_plot_init(expr = data$counts, sigs = sigs, counts = counts,
                         type = type, by = data$samples[[ID]],
                         xint = xint, yint = yint,
                         gene_id = gene_id)
  return(p)
})

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "ExpressionSet",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         xint = 1,
         yint = 1,
         gene_id = "SYMBOL") {

  p <- scatter_plot_init(expr = Biobase::exprs(data), sigs = sigs, counts = counts,
                         type = type, by = data[[ID]], xint = xint, yint = yint,
                         gene_id = gene_id)
  return(p)
})

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "Seurat",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         slot = "counts",
         xint = 1,
         yint = 1,
         gene_id = "SYMBOL") {

  p <- scatter_plot_init(expr = Seurat::GetAssayData(data, slot = slot),
                         sigs = sigs, counts = counts, type = type,
                         by = data@meta.data[[ID]], xint = xint, yint = yint,
                         gene_id = gene_id)
  return(p)
})

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "SummarizedExperiment",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         slot = "counts",
         xint = 1,
         yint = 1,
         gene_id = "SYMBOL") {

  p <- scatter_plot_init(expr = SummarizedExperiment::assay(data, slot),
                         sigs = sigs, counts = counts, type = type,
                         by = data@colData[[ID]], xint = xint, yint = yint,
                         gene_id = gene_id)
  return(p)
})

#' @rdname sig_scatter_plot
setMethod("sig_scatter_plot", signature(
  data = "list",
  sigs = "vector",
  ID = "character",
  type = "character"
),
function(data,
         sigs,
         ID,
         type,
         counts = TRUE,
         slot = "counts",
         xint = 1,
         yint = 1,
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

  p <- list()
  for (i in 1:length(data)) {
    p[[i]] <- sig_scatter_plot(data = data[[i]], sigs = sigs, counts = counts[i],
                               type = type[i], slot = slot[i], ID = ID[i],
                               xint = xint, yint = yint, gene_id = gene_id[i])
    if(!is.null(names(data)))
      p[[i]] <- p[[i]] + labs(subtitle = names(data)[i]) +
        theme(plot.subtitle = element_text(hjust = 0.5))
  }
  p <- patchwork::wrap_plots(p) +
    patchwork::plot_layout(guides = "collect")
  return(p)
})


