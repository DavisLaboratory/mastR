#' @include DE_functions.R plot.R
NULL

#' @title Get DE analysis result table with statistics
#'
#' @description This function uses edgeR and limma to get DE analysis results
#'   lists for multiple comparisons. Filter out low expressed genes and obtain
#'   DE statistics by using limma::voom and limma::treat, and also create an
#'   object `proc_data` to store processed data.
#'
#' @param data expression object
#' @param ID chr, column name of coldata to specify the DE comparisons
#' @param type pattern, specify the group of interest, e.g. NK
#' @param slot chr, specify which slot to use only for sce or seurat object,
#'             optional, default 'counts'
#' @param ... params for function [DE_analysis()]
#'
#' @return A list of DE result table of all comparisons.
#'
#' @examples
#' DE_tables <- get_DE_table(im_data_6, ID = "celltype:ch1", type = "NK")
#'
#'@export
setGeneric("get_DE_table",
           function(data,
                    ID,
                    type,
                    slot = "counts",
                    ...)
           standardGeneric("get_DE_table"))

#' @rdname get_DE_table
setMethod("get_DE_table", signature(
  data = 'DGEList',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         ...) {

  DGE <- data
  DGE$samples$group <- DGE$samples[[ID]]
  rm(data)

  ## standard DE analysis with edgeR and limma::voom pipeline
  tfit <- DE_analysis(dge = DGE,
                      ID = ID,
                      type = type,
                      ...)

  ## save DE result tables into list
  DE_table <- list()
  for (i in 1:ncol(tfit)) {
    ## use limma::topTreat() to get statistics of DEA
    DE_table[[i]] <- limma::topTreat(tfit, coef = i, number = Inf) |> na.omit()
  }
  names(DE_table) <- colnames(tfit)

  return(DE_table)
})

#' @rdname get_DE_table
setMethod("get_DE_table", signature(
  data = 'matrix',
  ID = 'vector',
  type = 'character'
),
function(data,
         ID,
         type,
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = ID)
  ID <- "group"
  rm(data)

  DE_table <- get_DE_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_DE_table
setMethod("get_DE_table", signature(
  data = 'Matrix',
  ID = 'vector',
  type = 'character'
),
function(data,
         ID,
         type,
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = ID)
  ID <- "group"
  rm(data)

  DE_table <- get_DE_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_DE_table
setMethod("get_DE_table", signature(
  data = 'ExpressionSet',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         ...) {

  expr <- Biobase::exprs(data)
  coldata <- Biobase::pData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DE_table <- get_DE_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_DE_table
setMethod("get_DE_table", signature(
  data = 'SummarizedExperiment',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         slot = "counts",
         ...) {

  expr <- SummarizedExperiment::assay(data, slot)
  coldata <- SummarizedExperiment::colData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DE_table <- get_DE_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_DE_table
setMethod("get_DE_table", signature(
  data = 'Seurat',
  ID = 'character',
  type = 'character'
),
function(data,
         ID,
         type,
         slot = "counts",
         ...) {

  expr <- Seurat::GetAssayData(data, slot = slot)
  coldata <- data@meta.data

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[ID]])
  ID <- make.names(ID)
  rm(data, expr, coldata)

  DE_table <- get_DE_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})
