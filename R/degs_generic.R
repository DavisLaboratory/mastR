#' @include DE_functions.R plot.R
NULL

#' @title Get differentially expressed genes by comparing specified groups
#'
#' @description This function uses edgeR and limma to get 'UP' and 'DOWN' DEG
#'   lists, for multiple comparisons, DEGs can be obtained from intersection of
#'   all comparsion DEGs or by using product of p value ranks for multiple
#'   comparisons. Filter out low expressed genes and extract DE genes by using
#'   limma::voom and limma::treat, and also create an object `proc_data` to
#'   store processed data.
#'
#' @inheritParams de_analysis
#' @param data expression object
#' @param group_col vector or character, specify the group factor or column name of
#'           coldata for DE comparisons
#' @param slot character, specify which slot to use only for sce or seurat
#'             object, optional, default 'counts'
#' @param batch vector of column name(s) or dataframe, specify the batch effect
#'              factor(s), default NULL
#' @param ... params for [de_analysis()] and [DEGs_RP()] or [DEGs_Group()]
#'
#' @return A list of 'UP' and 'DOWN' data frame of all differentially expressed
#'         genes. Both are ordered by rank product or 'Rank' variable if
#'         keep.top is NULL
#'
#' @examples
#' data("im_data_6")
#' DEGs <- get_degs(im_data_6, group_col = "celltype:ch1",
#'                  target_group = "NK", gene_id = "ENSEMBL")
#'
#'@export
setGeneric("get_degs",
           function(data,
                    group_col,
                    target_group,
                    normalize = TRUE,
                    feature_selection = c('auto', "rankproduct", "none"),
                    slot = "counts",
                    batch = NULL,
                    ...)
           standardGeneric("get_degs"))

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'DGEList',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         normalize = TRUE,
         feature_selection = c('auto', "rankproduct", "none"),
         batch = NULL,
         ...) {

  feature_selection <- match.arg(feature_selection)
  stopifnot("Please provide column names as batch!" =
              is.null(batch) | is.vector(batch))

  DGE <- data
  DGE$samples$group <- DGE$samples[[group_col]]
  rm(data)

  ## standard DE analysis with edgeR and limma::voom pipeline
  voom_res <- de_analysis(
    dge = DGE,
    group_col = group_col,
    target_group = target_group,
    normalize = normalize,
    feature_selection = feature_selection,
    batch = batch,
    ...
  )

  ## assemble DEGs from comparisons by Rank Product or simply intersect/union
  if(feature_selection == "auto") {
    DEGs <- DEGs_RP(tfit = voom_res$tfit, ...)
    if(all(lengths(DEGs) < 5))
      DEGs <- DEGs_Group(tfit = voom_res$tfit, ...)
  }else if(feature_selection == "rankproduct") {
    DEGs <- DEGs_RP(tfit = voom_res$tfit, ...)
  }else {
    DEGs <- DEGs_Group(tfit = voom_res$tfit, ...)
  }

  ## save the processed DGEList
  DEGs$proc_data <- voom_res$proc_data

  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'matrix',
  group_col = 'vector',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         normalize = TRUE,
         feature_selection = c('auto', "rankproduct", "none"),
         batch = NULL,
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = group_col)
  if(!is.null(batch)) {
    DGE$samples <- data.frame(group = group_col, batch)
    batch <- colnames(DGE$samples)[-1]
  }
  group_col <- "group"
  rm(data)

  DEGs <- get_degs(data = DGE, group_col = group_col,
                   target_group = target_group,
                   normalize = normalize,
                   feature_selection = feature_selection,
                   batch = batch,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'Matrix',
  group_col = 'vector',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         normalize = TRUE,
         feature_selection = c('auto', "rankproduct", "none"),
         batch = NULL,
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = group_col)
  if(!is.null(batch)) {
    DGE$samples <- data.frame(group = group_col, batch)
    batch <- colnames(DGE$samples)[-1]
  }
  group_col <- "group"
  rm(data)

  DEGs <- get_degs(data = DGE, group_col = group_col,
                   target_group = target_group,
                   normalize = normalize,
                   feature_selection = feature_selection,
                   batch = batch,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'ExpressionSet',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         normalize = TRUE,
         feature_selection = c('auto', "rankproduct", "none"),
         batch = NULL,
         ...) {

  expr <- Biobase::exprs(data)
  coldata <- Biobase::pData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[group_col]])
  if(!is.null(batch)) batch <- make.names(batch)
  group_col <- make.names(group_col)
  rm(data, expr, coldata)

  DEGs <- get_degs(data = DGE, group_col = group_col,
                   target_group = target_group,
                   normalize = normalize,
                   feature_selection = feature_selection,
                   batch = batch,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'SummarizedExperiment',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         normalize = TRUE,
         feature_selection = c('auto', "rankproduct", "none"),
         slot = "counts",
         batch = NULL,
         ...) {

  expr <- SummarizedExperiment::assay(data, slot)
  coldata <- SummarizedExperiment::colData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[group_col]])
  if(!is.null(batch)) batch <- make.names(batch)
  group_col <- make.names(group_col)
  rm(data, expr, coldata)

  DEGs <- get_degs(data = DGE, group_col = group_col,
                   target_group = target_group,
                   normalize = normalize,
                   feature_selection = feature_selection,
                   batch = batch,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'Seurat',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         normalize = TRUE,
         feature_selection = c('auto', "rankproduct", "none"),
         slot = "counts",
         batch = NULL,
         ...) {

  expr <- Seurat::GetAssayData(data, slot = slot)
  coldata <- data@meta.data

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[group_col]])
  if(!is.null(batch)) batch <- make.names(batch)
  group_col <- make.names(group_col)
  rm(data, expr, coldata)

  DEGs <- get_degs(data = DGE, group_col = group_col,
                   target_group = target_group,
                   normalize = normalize,
                   feature_selection = feature_selection,
                   batch = batch,
                   ...)
  return(DEGs)
})


#' @title DE analysis pipeline
#'
#' @description Standard DE analysis by using edgeR and limma::voom pipeline
#'
#' @param dge DGEList object for DE analysis, including expr and samples info
#' @param group_col character, column name of coldata to specify the DE comparisons
#' @param target_group pattern, specify the group of interest, e.g. NK
#' @param normalize logical, if the expr in data is raw counts needs to be normalized
#' @param feature_selection one of 'auto', "rankproduct" or "none", choose if to
#'                          use rank product or not to select DEGs from multiple
#'                          comparisons of DE analysis, default 'auto' uses
#'                          'rankproduct' but change to 'none' if final genes < 5
#'                          for both UP and DOWN.
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
#' @param batch vector of character, column name(s) of coldata to be treated as
#'              batch effect factor, default NULL
#' @param ... omitted
#'
#' @return MArrayLM object generated by [limma::treat()]
#'
#' @export
de_analysis <- function(dge,
                        group_col,
                        target_group,
                        normalize = TRUE,
                        feature_selection = c('auto', "rankproduct", "none"),
                        group = FALSE,
                        filter = c(10, 10),
                        plot = FALSE,
                        lfc = 0,
                        p = 0.05,
                        markers = NULL,
                        gene_id = "SYMBOL",
                        batch = NULL,
                        ...) {

  stopifnot(is.logical(normalize), is.character(feature_selection),
            is.logical(group), is.numeric(filter),
            is.logical(plot), is.numeric(lfc),
            is.numeric(p), is.character(gene_id))

  feature_selection <- match.arg(feature_selection)

  ## filter low counts genes
  keep <- filterGenes(dge = dge,
                      group_col = group_col,
                      filter = filter,
                      normalize = normalize,
                      markers = markers,
                      gene_id = gene_id)

  ## plot distribution of unfiltered and filtered data
  if(plot) {
    plot_density(dge = dge, group_col = group_col, keep = keep,
                 normalize = normalize, filter = filter[1])
    plot_rle(dge = dge, group_col = group_col,
             keep = keep, normalize = normalize)
    plot_MDS(dge = dge, group_col = group_col,
             keep = keep, normalize = normalize)
  }

  ## DE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  ## normalization for raw counts data
  if(normalize)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

  ## voom linear fit by limma
  ## return a list of tfit and a processed DGEList
  res <- voom_lm_fit(dge = dge, group_col = group_col,
                     target_group = target_group,
                     feature_selection = feature_selection,
                     group = group, plot = plot, normalize = normalize,
                     lfc = lfc, p = p, batch = batch)

  return(res)
}
