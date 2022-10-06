#' @include DE_functions.R plot.R
NULL

#' @title Get differentially expressed genes by comparing specified groups
#'
#' @description This function uses edgeR and limma to get 'UP' and 'DOWN' DEG
#'   lists, for multiple comparisons, DEGs can be obtained from intersection of
#'   all comparsion DEGs or by using product of p value ranks for multiple
#'   comparisons. Filter out low expressed genes and extract DE genes by using
#'   limma::voom and limma::treat, and also create an object `proc_data` to store
#'   processed data.
#'
#' @param data expression object
#' @param ID vector or chr, to specify the group of DE comparisons
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
#' @param nperm num, permutation runs of simulating the distribution
#' @param rand NULL or num, default 123, if NULL the permutation would be random
#' @param assemble 'intersect' or 'union', whether to select intersected or
#'                  union genes of different comparisons, default 'intersect'
#' @param Rank chr, the variable for ranking DEGs, can be 'logFC', 'adj.P.Val'...
#'             default 'adj.P.Val'
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                kept anyway after filtration. Default 'NULL' means no special
#'                genes need to be kept.
#' @param gene_id chr, specify the gene ID type of rownames of expression data
#'                when markers is not NULL, could be one of 'ENSEMBL', 'SYMBOL',
#'                'ENTREZ'..., default 'SYMBOL'
#' @param keep.top NULL or num, whether to keep top n DEGs of specific comparison
#' @param keep.group NULL or pattern, specify the top DEGs of which comparison
#'                   or group to be kept
#' @param slot chr, specify which slot to use only for sce or seurat object,
#'             optional, default 'counts'
#'
#' @return A list of 'UP' and 'DOWN' data frame of all differentially expressed
#'         genes. Both are ordered by rank product or 'Rank' variable if
#'         keep.top is NULL
#'
#' @examples
#' DEGs <- get_DEGs(im_data_6, ID = "celltype:ch1",
#'                  type = "NK", gene_id = "ENSEMBL")
#'
#'@export
setGeneric("get_DEGs",
           function(data,
                    ID,
                    type,
                    counts = TRUE,
                    method = "RP",
                    group = FALSE,
                    filter = c(10, 10),
                    plot = FALSE,
                    lfc = 0,
                    p = 0.05,
                    nperm = 1e5,
                    rand = 123,
                    assemble = "intersect",
                    Rank = "adj.P.Val",
                    markers = NULL,
                    gene_id = "SYMBOL",
                    keep.top = NULL,
                    keep.group = NULL,
                    slot = "counts")
           standardGeneric("get_DEGs"))

#' @rdname get_DEGs
setMethod("get_DEGs", signature(
  data = 'DGEList',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = "RP",
         group = FALSE,
         filter = c(10, 10),
         plot = FALSE,
         lfc = 0,
         p = 0.05,
         nperm = 1e5,
         rand = 123,
         assemble = "intersect",
         Rank = "adj.P.Val",
         markers = NULL,
         gene_id = "SYMBOL",
         keep.top = NULL,
         keep.group = NULL) {

  stopifnot(is.numeric(filter), is.numeric(lfc), is.numeric(p),
            is.numeric(nperm), is.character(gene_id))

  DGE <- data
  DGE$samples$group <- DGE$samples[[ID]]
  rm(data)

  ## standard DE analysis with edgeR and limma::voom pipeline
  tfit <- DE_analysis(dge = DGE,
                      ID = ID,
                      type = type,
                      counts = counts,
                      method = method,
                      group = group,
                      filter = filter,
                      plot = plot,
                      lfc = lfc,
                      p = p,
                      markers = markers,
                      gene_id = gene_id)

  ## assemble DEGs from comparisons by Rank Product or simple groups
  if(method == "RP"){
    DEGs <- DEGs_RP(tfit = tfit, lfc = lfc, p = p, assemble = assemble,
                    Rank = Rank, rand = rand, nperm = nperm,
                    keep.top = keep.top, keep.group = keep.group)
  }else if(method == "Group") {
    DEGs <- DEGs_Group(tfit = tfit, lfc = lfc, p = p,
                       assemble = assemble, Rank = Rank,
                       keep.top = keep.top, keep.group = keep.group)
  }else stop("Please provide a valid method, either 'RP' or 'Group'!")

  return(DEGs)
})

#' @rdname get_DEGs
setMethod("get_DEGs", signature(
  data = 'matrix',
  ID = 'vector',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = "RP",
         group = FALSE,
         filter = c(10, 10),
         plot = FALSE,
         lfc = 0,
         p = 0.05,
         nperm = 1e5,
         rand = 123,
         assemble = "intersect",
         Rank = "adj.P.Val",
         markers = NULL,
         gene_id = "SYMBOL",
         keep.top = NULL,
         keep.group = NULL) {

  stopifnot(is.numeric(filter), is.numeric(lfc), is.numeric(p),
            is.numeric(nperm), is.character(gene_id))

  DGE <- edgeR::DGEList(counts = data, group = ID)
  ID <- "group"
  rm(data)

  DEGs <- get_DEGs(data = DGE, ID = ID, type = type, counts = counts,
                   method = method, group = group, filter = filter,
                   plot = plot, lfc = lfc, p = p, nperm = nperm,
                   rand = rand, assemble = assemble, Rank = Rank,
                   markers = markers, gene_id = gene_id,
                   keep.top = keep.top, keep.group = keep.group)
  return(DEGs)
})

#' @rdname get_DEGs
setMethod("get_DEGs", signature(
  data = 'Matrix',
  ID = 'vector',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = "RP",
         group = FALSE,
         filter = c(10, 10),
         plot = FALSE,
         lfc = 0,
         p = 0.05,
         nperm = 1e5,
         rand = 123,
         assemble = "intersect",
         Rank = "adj.P.Val",
         markers = NULL,
         gene_id = "SYMBOL",
         keep.top = NULL,
         keep.group = NULL) {

  stopifnot(is.numeric(filter), is.numeric(lfc), is.numeric(p),
            is.numeric(nperm), is.character(gene_id))

  DGE <- edgeR::DGEList(counts = data, group = ID)
  ID <- "group"
  rm(data)

  DEGs <- get_DEGs(data = DGE, ID = ID, type = type, counts = counts,
                   method = method, group = group, filter = filter,
                   plot = plot, lfc = lfc, p = p, nperm = nperm,
                   rand = rand, assemble = assemble, Rank = Rank,
                   markers = markers, gene_id = gene_id,
                   keep.top = keep.top, keep.group = keep.group)
  return(DEGs)
})

#' @rdname get_DEGs
setMethod("get_DEGs", signature(
  data = 'ExpressionSet',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = "RP",
         group = FALSE,
         filter = c(10, 10),
         plot = FALSE,
         lfc = 0,
         p = 0.05,
         nperm = 1e5,
         rand = 123,
         assemble = "intersect",
         Rank = "adj.P.Val",
         markers = NULL,
         gene_id = "SYMBOL",
         keep.top = NULL,
         keep.group = NULL) {

  stopifnot(is.numeric(filter), is.numeric(lfc), is.numeric(p),
            is.numeric(nperm), is.character(gene_id))
  expr <- Biobase::exprs(data)
  coldata <- Biobase::pData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DEGs <- get_DEGs(data = DGE, ID = ID, type = type, counts = counts,
                   method = method, group = group, filter = filter,
                   plot = plot, lfc = lfc, p = p, nperm = nperm,
                   rand = rand, assemble = assemble, Rank = Rank,
                   markers = markers, gene_id = gene_id,
                   keep.top = keep.top, keep.group = keep.group)
  return(DEGs)
})

#' @rdname get_DEGs
setMethod("get_DEGs", signature(
  data = 'SummarizedExperiment',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = "RP",
         group = FALSE,
         filter = c(10, 10),
         plot = FALSE,
         lfc = 0,
         p = 0.05,
         nperm = 1e5,
         rand = 123,
         assemble = "intersect",
         Rank = "adj.P.Val",
         markers = NULL,
         gene_id = "SYMBOL",
         keep.top = NULL,
         keep.group = NULL,
         slot = "counts") {

  stopifnot(is.numeric(filter), is.numeric(lfc), is.numeric(p),
            is.numeric(nperm), is.character(gene_id), is.character(slot))
  expr <- SummarizedExperiment::assay(data, slot)
  coldata <- SummarizedExperiment::colData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DEGs <- get_DEGs(data = DGE, ID = ID, type = type, counts = counts,
                   method = method, group = group, filter = filter,
                   plot = plot, lfc = lfc, p = p, nperm = nperm,
                   rand = rand, assemble = assemble, Rank = Rank,
                   markers = markers, gene_id = gene_id,
                   keep.top = keep.top, keep.group = keep.group)
  return(DEGs)
})

#' @rdname get_DEGs
setMethod("get_DEGs", signature(
  data = 'Seurat',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         counts = TRUE,
         method = "RP",
         group = FALSE,
         filter = c(10, 10),
         plot = FALSE,
         lfc = 0,
         p = 0.05,
         nperm = 1e5,
         rand = 123,
         assemble = "intersect",
         Rank = "adj.P.Val",
         markers = NULL,
         gene_id = "SYMBOL",
         keep.top = NULL,
         keep.group = NULL,
         slot = "counts") {

  stopifnot(is.numeric(filter), is.numeric(lfc), is.numeric(p),
            is.numeric(nperm), is.character(gene_id), is.character(slot))
  expr <- Seurat::GetAssayData(data, slot = slot)
  coldata <- data@meta.data

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DEGs <- get_DEGs(data = DGE, ID = ID, type = type, counts = counts,
                   method = method, group = group, filter = filter,
                   plot = plot, lfc = lfc, p = p, nperm = nperm,
                   rand = rand, assemble = assemble, Rank = Rank,
                   markers = markers, gene_id = gene_id,
                   keep.top = keep.top, keep.group = keep.group)
  return(DEGs)
})


#' @title DE analysis pipeline
#'
#' @description Standard DE analysis by using edgeR and limma::voom pipeline
#'
#' @param dge DGEList object for DE analysis, including expr and samples info
#' @param ID chr, column name of coldata to specify the DE comparisons
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
#' @param gene_id chr, specify the gene ID type of rownames of expression data
#'                when markers is not NULL, could be one of 'ENSEMBL', 'SYMBOL',
#'                'ENTREZ'..., default 'SYMBOL'
#'
#' @return MArrayLM object generated by [limma::treat()]
#'
#' @export
DE_analysis <- function(dge,
                        ID,
                        type,
                        counts = TRUE,
                        method = "RP",
                        group = FALSE,
                        filter = c(10, 10),
                        plot = FALSE,
                        lfc = 0,
                        p = 0.05,
                        markers = NULL,
                        gene_id = "SYMBOL") {

  stopifnot(is.logical(counts), is.character(method), is.logical(group),
            is.numeric(filter), is.logical(plot), is.numeric(lfc),
            is.numeric(p), is.character(gene_id))

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
  dge <- dge[keep, ,keep.lib.sizes = FALSE]
  ## normalization for raw counts data
  if(counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

  ## voom linear fit by limma
  tfit <- voom_lm_fit(dge = dge, ID = ID, type = type, method = method,
                      group = group, plot = plot, counts = counts,
                      lfc = lfc, p = p)

  return(tfit)
}
