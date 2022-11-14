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
#' @param ID vector or character, specify the group factor or column name of
#'           coldata for DE comparisons
#' @param slot character, specify which slot to use only for sce or seurat
#'             object, optional, default 'counts'
#' @param ... params for [de_analysis()] and [DEGs_RP()] or [DEGs_Group()]
#'
#' @return A list of 'UP' and 'DOWN' data frame of all differentially expressed
#'         genes. Both are ordered by rank product or 'Rank' variable if
#'         keep.top is NULL
#'
#' @examples
#' DEGs <- get_degs(im_data_6, ID = "celltype:ch1",
#'                  type = "NK", gene_id = "ENSEMBL")
#'
#'@export
setGeneric("get_degs",
           function(data,
                    ID,
                    type,
                    counts = TRUE,
                    method = c("RP", "Group"),
                    slot = "counts",
                    ...)
           standardGeneric("get_degs"))

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'DGEList',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = c("RP", "Group"),
         ...) {

  method <- match.arg(method)

  DGE <- data
  DGE$samples$group <- DGE$samples[[ID]]
  rm(data)

  ## standard DE analysis with edgeR and limma::voom pipeline
  tfit <- de_analysis(dge = DGE,
                      ID = ID,
                      type = type,
                      counts = counts,
                      method = method,
                      ...)

  ## assemble DEGs from comparisons by Rank Product or simple groups
  if(method == "RP") {
    DEGs <- DEGs_RP(tfit = tfit, ...)
  }else if(method == "Group") {
    DEGs <- DEGs_Group(tfit = tfit, ...)
  }else stop("Please provide a valid method, either 'RP' or 'Group'!")

  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'matrix',
  ID = 'vector',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = c("RP", "Group"),
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = ID)
  ID <- "group"
  rm(data)

  DEGs <- get_degs(data = DGE, ID = ID, type = type,
                   counts = counts, method = method,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'Matrix',
  ID = 'vector',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = c("RP", "Group"),
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = ID)
  ID <- "group"
  rm(data)

  DEGs <- get_degs(data = DGE, ID = ID, type = type,
                   counts = counts, method = method,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'ExpressionSet',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = c("RP", "Group"),
         ...) {

  expr <- Biobase::exprs(data)
  coldata <- Biobase::pData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DEGs <- get_degs(data = DGE, ID = ID, type = type,
                   counts = counts, method = method,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'SummarizedExperiment',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = c("RP", "Group"),
         slot = "counts",
         ...) {

  expr <- SummarizedExperiment::assay(data, slot)
  coldata <- SummarizedExperiment::colData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DEGs <- get_degs(data = DGE, ID = ID, type = type,
                   counts = counts, method = method,
                   ...)
  return(DEGs)
})

#' @rdname get_degs
setMethod("get_degs", signature(
  data = 'Seurat',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = c("RP", "Group"),
         slot = "counts",
         ...) {

  expr <- Seurat::GetAssayData(data, slot = slot)
  coldata <- data@meta.data

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DEGs <- get_degs(data = DGE, ID = ID, type = type,
                   counts = counts, method = method,
                   ...)
  return(DEGs)
})


#' @title DE analysis pipeline
#'
#' @description Standard DE analysis by using edgeR and limma::voom pipeline
#'
#' @param dge DGEList object for DE analysis, including expr and samples info
#' @param ID character, column name of coldata to specify the DE comparisons
#' @param type pattern, specify the group of interest, e.g. NK
#' @param counts logical, if the expr in data is raw counts data
#' @param method either 'RP' or 'Group', choose whether to use rank product or
#'               group other subsets for multiple comparisons for DE analysis,
#'               default 'RP'
#' @param group logical, TRUE to separate samples into only 2 groups:
#'              `type`` and 'Others'; FALSE to set each level as a group
#' @param filter a vector of 2 numbers, filter condition to remove low expression
#'               genes, the 1st for min.counts (if counts = TRUE) or CPM/TPM
#'               (if counts = F), the 2nd for samples size 'large.n'
#' @param plot logical, if to make plots to show QC before and after filtration
#' @param lfc num, cutoff of logFC for DE analysis
#' @param p num, cutoff of p value for DE analysis and permutation test if
#'          method = "RP"
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                kept anyway after filtration. Default 'NULL' means no special
#'                genes need to be kept.
#' @param gene_id character, specify the gene ID type of rownames of expression data
#'                when markers is not NULL, could be one of 'ENSEMBL', 'SYMBOL',
#'                'ENTREZ'..., default 'SYMBOL'
#' @param ... omitted
#'
#' @return MArrayLM object generated by [limma::treat()]
#'
#' @export
de_analysis <- function(dge,
                        ID,
                        type,
                        counts = TRUE,
                        method = c("RP", "Group"),
                        group = FALSE,
                        filter = c(10, 10),
                        plot = FALSE,
                        lfc = 0,
                        p = 0.05,
                        markers = NULL,
                        gene_id = "SYMBOL",
                        ...) {

  stopifnot(is.logical(counts), is.character(method), is.logical(group),
            is.numeric(filter), is.logical(plot), is.numeric(lfc),
            is.numeric(p), is.character(gene_id))

  method <- match.arg(method)

  ## filter low counts genes
  keep <- filterGenes(dge = dge, ID = ID, filter = filter, counts = counts,
                      markers = markers, gene_id = gene_id)

  ## plot distribution of unfiltered and filtered data
  if(plot) {
    plot_density(dge = dge, ID = ID, keep = keep,
                 counts = counts, filter = filter[1])
    plot_rle(dge = dge, ID = ID, keep = keep, counts = counts)
    plot_MDS(dge = dge, ID = ID, keep = keep, counts = counts)
  }

  ## DE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  ## normalization for raw counts data
  if(counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

  ## voom linear fit by limma
  tfit <- voom_lm_fit(dge = dge, ID = ID, type = type, method = method,
                      group = group, plot = plot, counts = counts,
                      lfc = lfc, p = p)

  return(tfit)
}
