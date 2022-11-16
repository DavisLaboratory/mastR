#' @include DE_functions.R plot.R
NULL

#' @title Get DE analysis result table(s) with statistics
#'
#' @description This function uses edgeR and limma to get DE analysis results
#'   lists for multiple comparisons. Filter out low expressed genes and obtain
#'   DE statistics by using limma::voom and limma::treat, and also create an
#'   object `proc_data` to store processed data.
#'
#' @inheritParams get_degs
#' @param ... params for function [de_analysis()]
#'
#' @return A list of DE result table of all comparisons.
#'
#' @examples
#' DE_tables <- get_de_table(im_data_6, ID = "celltype:ch1", type = "NK")
#'
#'@export
setGeneric("get_de_table",
           function(data,
                    ID,
                    type,
                    slot = "counts",
                    ...)
           standardGeneric("get_de_table"))

#' @rdname get_de_table
setMethod("get_de_table", signature(
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
  voom_res <- de_analysis(
    dge = DGE,
    ID = ID,
    type = type,
    ...
  )

  ## save DE result tables into list
  DE_table <- list()
  for (i in seq_len(ncol(voom_res$tfit))) {
    ## use limma::topTreat() to get statistics of DEA
    DE_table[[i]] <- limma::topTreat(voom_res$tfit, coef = i, number = Inf) |> na.omit()
  }
  names(DE_table) <- colnames(voom_res$tfit)

  ## save the processed DGEList
  DE_table$proc_data <- voom_res$proc_data

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
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

  DE_table <- get_de_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
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

  DE_table <- get_de_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
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

  DE_table <- get_de_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
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

  DE_table <- get_de_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
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

  DE_table <- get_de_table(data = DGE, ID = ID, type = type, ...)

  return(DE_table)
})
