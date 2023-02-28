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
#' data("im_data_6")
#' DE_tables <- get_de_table(im_data_6, group_col = "celltype:ch1", target_group = "NK")
#'
#'@export
setGeneric("get_de_table",
           function(data,
                    group_col,
                    target_group,
                    slot = "counts",
                    ...)
           standardGeneric("get_de_table"))

#' @rdname get_de_table
setMethod("get_de_table", signature(
  data = 'DGEList',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         ...) {

  DGE <- data
  DGE$samples$group <- DGE$samples[[group_col]]
  rm(data)

  ## standard DE analysis with edgeR and limma::voom pipeline
  voom_res <- de_analysis(
    dge = DGE,
    group_col = group_col,
    target_group = target_group,
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
  group_col = 'vector',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = group_col)
  group_col <- "group"
  rm(data)

  DE_table <- get_de_table(data = DGE, group_col = group_col, target_group = target_group, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
  data = 'Matrix',
  group_col = 'vector',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         ...) {

  DGE <- edgeR::DGEList(counts = data, group = group_col)
  group_col <- "group"
  rm(data)

  DE_table <- get_de_table(data = DGE, group_col = group_col, target_group = target_group, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
  data = 'ExpressionSet',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         ...) {

  expr <- Biobase::exprs(data)
  coldata <- Biobase::pData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[group_col]])
  group_col <- make.names(group_col)
  rm(data, expr, coldata)

  DE_table <- get_de_table(data = DGE, group_col = group_col, target_group = target_group, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
  data = 'SummarizedExperiment',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         slot = "counts",
         ...) {

  expr <- SummarizedExperiment::assay(data, slot)
  coldata <- SummarizedExperiment::colData(data)

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[group_col]])
  group_col <- make.names(group_col)
  rm(data, expr, coldata)

  DE_table <- get_de_table(data = DGE, group_col = group_col, target_group = target_group, ...)

  return(DE_table)
})

#' @rdname get_de_table
setMethod("get_de_table", signature(
  data = 'Seurat',
  group_col = 'character',
  target_group = 'character'
),
function(data,
         group_col,
         target_group,
         slot = "counts",
         ...) {

  expr <- Seurat::GetAssayData(data, slot = slot)
  coldata <- data@meta.data

  DGE <- edgeR::DGEList(counts = expr,
                        samples = coldata,
                        group = coldata[[group_col]])
  group_col <- make.names(group_col)
  rm(data, expr, coldata)

  DE_table <- get_de_table(data = DGE, group_col = group_col, target_group = target_group, ...)

  return(DE_table)
})
