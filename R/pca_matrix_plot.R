#' @include plot.R
NULL

#' @title Make a matrix plot of PCA with top PCs
#'
#' @inheritParams pca_matrix_plot_init
#' @param data expression data, can be matrix, eSet, seurat...
#' @param slot character, specify the slot name of seurat or sce object,
#'             optional
#'
#' @return matrix plot of PCA
#'
#' @examples
#' pca_matrix_plot(data = im_data_6)
#'
#' @export
setGeneric("pca_matrix_plot",
           function(data,
                    features = "all",
                    slot = "counts",
                    counts = TRUE,
                    group_by = NULL,
                    scale = TRUE,
                    n = 4,
                    loading = FALSE,
                    n_loadings = 10,
                    gene_id = "SYMBOL")
             standardGeneric("pca_matrix_plot"))

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = "matrix"
),
function(data,
         features = "all",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  p <- pca_matrix_plot_init(data = data, features = features,
                            counts = counts, group_by = group_by,
                            scale = scale, n = n, loading = loading,
                            n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = "Matrix"
),
function(data,
         features = "all",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  p <- pca_matrix_plot_init(data = data, features = features,
                            counts = counts, group_by = group_by,
                            scale = scale, n = n, loading = loading,
                            n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = "data.frame"
),
function(data,
         features = "all",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  p <- pca_matrix_plot_init(data = data |> as.matrix(), features = features,
                            counts = counts, group_by = group_by,
                            scale = scale, n = n, loading = loading,
                            n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = 'ExpressionSet'
),
function(data,
         features = "all",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  if(is.character(group_by))
    group_by <- data[[group_by]]
  data <- Biobase::exprs(data)

  p <- pca_matrix_plot(data = data, features = features,
                       counts = counts, group_by = group_by,
                       scale = scale, n = n, loading = loading,
                       n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = 'DGEList'
),
function(data,
         features = "all",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  if(is.character(group_by))
    group_by <- data$samples[[group_by]]
  data <- data$counts

  p <- pca_matrix_plot(data = data, features = features,
                       counts = counts, group_by = group_by,
                       scale = scale, n = n, loading = loading,
                       n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = 'SummarizedExperiment'
),
function(data,
         features = "all",
         slot = "counts",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  if(is.character(group_by))
    group_by <- data@colData[[group_by]]
  data <- SummarizedExperiment::assay(data, slot)

  p <- pca_matrix_plot(data = data, features = features,
                       counts = counts, group_by = group_by,
                       scale = scale, n = n, loading = loading,
                       n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})

#' @rdname pca_matrix_plot
setMethod("pca_matrix_plot", signature(
  data = 'Seurat'
),
function(data,
         features = "all",
         slot = "counts",
         counts = TRUE,
         group_by = NULL,
         scale = TRUE,
         n = 4,
         loading = FALSE,
         n_loadings = 10,
         gene_id = "SYMBOL") {

  if(is.character(group_by))
    group_by <- data@meta.data[[group_by]]
  data <- Seurat::GetAssayData(data, slot = slot)

  p <- pca_matrix_plot(data = data, features = features,
                       counts = counts, group_by = group_by,
                       scale = scale, n = n, loading = loading,
                       n_loadings = n_loadings, gene_id = gene_id)
  return(p)
})
