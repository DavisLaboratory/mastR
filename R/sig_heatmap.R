#' @include plot.R
#' @import ggplot2 patchwork
NULL

#' Heatmap original markers and screened signature
#'
#' Compare the heatmap before and after screening.
#'
#' @inheritParams sig_gseaplot
#' @inheritParams sig_boxplot
#' @param markers a vector of gene names, listed the gene symbols of original
#'                markers pool
#' @param scale could be one of 'none' (default), 'row' or 'column'
#' @param ranks_plot logical, if to use ranks instead of expression of genes to
#'                   draw heatmap
#' @param ... params for [ComplexHeatmap::Heatmap()]
#'
#' @return patchwork object of pheatmaps
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' sig_heatmap(data = im_data_6, sigs = NK_markers$HGNC_Symbol[1:10],
#'             group_col = "celltype:ch1",
#'             gene_id = "ENSEMBL")
#'
#' @export
setGeneric("sig_heatmap",
           function(data,
                    sigs,
                    group_col,
                    markers,
                    normalize = FALSE,
                    scale = c("none", "row", "column"),
                    gene_id = "SYMBOL",
                    ranks_plot = FALSE,
                    slot = "counts",
                    ...)
             standardGeneric("sig_heatmap"))

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'matrix',
  sigs = 'character',
  group_col = 'vector',
  markers = 'missing'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         ...) {

  scale <- match.arg(scale)
  stopifnot(is.logical(normalize),
            is.character(gene_id),
            is.logical(ranks_plot))

  p <- heatmap_init(
    expr = data,
    sigs = list(Sig = sigs),
    by = factor(group_col),
    normalize = normalize,
    scale = scale,
    gene_id = gene_id,
    ranks_plot = ranks_plot,
    ...
  )

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'matrix',
  sigs = 'character',
  group_col = 'vector',
  markers = 'vector'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         ...) {
  scale <- match.arg(scale)
  stopifnot(is.logical(normalize),
            is.character(gene_id),
            is.logical(ranks_plot))

  ## build sigs
  sigs <- list(Original = setdiff(markers, sigs), Screen_Sig = sigs)

  p1 <- heatmap_init(
    expr = data,
    sigs = sigs,
    by = factor(group_col),
    normalize = normalize,
    scale = scale,
    gene_id = gene_id,
    ranks_plot = ranks_plot,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    name = "expr_orig",
    column_title = "Expression of Original Markers",
    ...
  )
  p2 <- heatmap_init(
    expr = data,
    sigs = sigs["Screen_Sig"],
    by = factor(group_col),
    normalize = normalize,
    scale = scale,
    gene_id = gene_id,
    ranks_plot = ranks_plot,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    name = "expr_scr",
    column_title = "Expression of Screened Signature",
    ...
  )
  p2@top_annotation <- p1@top_annotation
  ## convert to ggplot
  p1 <- patchwork::wrap_elements(grid::grid.grabExpr(ComplexHeatmap::draw(p1)))
  p2 <- patchwork::wrap_elements(grid::grid.grabExpr(ComplexHeatmap::draw(p2)))

  p <- patchwork::wrap_plots(p1, p2)
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'matrix',
  sigs = 'list',
  group_col = 'vector',
  markers = 'missing'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         ...) {

  scale <- match.arg(scale)
  stopifnot(is.logical(normalize),
            is.character(gene_id),
            is.logical(ranks_plot))

  p <- heatmap_init(
    expr = data,
    sigs = sigs,
    by = factor(group_col),
    normalize = normalize,
    scale = scale,
    gene_id = gene_id,
    ranks_plot = ranks_plot,
    ...
  )

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'Matrix',
  group_col = 'vector'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         ...) {

  stopifnot(is.logical(normalize), is.character(scale),
            is.character(gene_id), is.logical(ranks_plot))

  if(missing(markers)) {
    p <- sig_heatmap(data = as.matrix(data), sigs = sigs,
                     group_col = group_col,
                     normalize = normalize, scale = scale,
                     gene_id = gene_id, ranks_plot = ranks_plot, ...)
  } else {
    p <- sig_heatmap(data = as.matrix(data), sigs = sigs,
                     group_col = group_col, markers = markers,
                     normalize = normalize, scale = scale,
                     gene_id = gene_id, ranks_plot = ranks_plot, ...)
  }

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'data.frame',
  group_col = 'vector'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         ...) {

  if(missing(markers)) {
    p <- sig_heatmap(data = as.matrix(data), sigs = sigs,
                     group_col = group_col,
                     normalize = normalize, scale = scale,
                     gene_id = gene_id, ranks_plot = ranks_plot, ...)
  } else {
    p <- sig_heatmap(data = as.matrix(data), sigs = sigs,
                     group_col = group_col, markers = markers,
                     normalize = normalize, scale = scale,
                     gene_id = gene_id, ranks_plot = ranks_plot, ...)
  }

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'DGEList',
  group_col = 'character'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         ...) {

  if(missing(markers)) {
    p <- sig_heatmap(data = data[[slot]], sigs = sigs,
                     group_col = data$samples[[group_col]],
                     normalize = normalize,
                     scale = scale, gene_id = gene_id,
                     ranks_plot = ranks_plot, ...)
  } else {
    p <- sig_heatmap(data = data[[slot]], sigs = sigs,
                     group_col = data$samples[[group_col]],
                     markers = markers, normalize = normalize,
                     scale = scale, gene_id = gene_id,
                     ranks_plot = ranks_plot, ...)
  }

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'ExpressionSet',
  group_col = 'character'
),
function(data,
         sigs,
         group_col,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         ...) {

  if(missing(markers)) {
    p <- sig_heatmap(data = Biobase::exprs(data), sigs = sigs,
                     group_col = data[[group_col]],
                     normalize = normalize, scale = scale,
                     gene_id = gene_id, ranks_plot = ranks_plot, ...)
  } else {
    p <- sig_heatmap(data = Biobase::exprs(data), sigs = sigs,
                     group_col = data[[group_col]], markers = markers,
                     normalize = normalize, scale = scale,
                     gene_id = gene_id, ranks_plot = ranks_plot, ...)
  }

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'Seurat',
  group_col = 'character'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         ...) {

  if(missing(markers)) {
    p <- sig_heatmap(data = Seurat::GetAssayData(data, slot = slot),
                     sigs = sigs, group_col = data@meta.data[[group_col]],
                     normalize = normalize,
                     scale = scale, gene_id = gene_id,
                     ranks_plot = ranks_plot, ...)
  } else {
    p <- sig_heatmap(data = Seurat::GetAssayData(data, slot = slot),
                     sigs = sigs, group_col = data@meta.data[[group_col]],
                     markers = markers, normalize = normalize,
                     scale = scale, gene_id = gene_id,
                     ranks_plot = ranks_plot, ...)
  }

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'SummarizedExperiment',
  group_col = 'character'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         ...) {

  if(missing(markers)) {
    p <- sig_heatmap(data = SummarizedExperiment::assay(data, slot),
                     sigs = sigs, group_col = data@colData[[group_col]],
                     normalize = normalize,
                     scale = scale, gene_id = gene_id,
                     ranks_plot = ranks_plot, ...)
  } else {
    p <- sig_heatmap(data = SummarizedExperiment::assay(data, slot),
                     sigs = sigs, group_col = data@colData[[group_col]],
                     markers = markers, normalize = normalize,
                     scale = scale, gene_id = gene_id,
                     ranks_plot = ranks_plot, ...)
  }

  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'list',
  group_col = 'character'
),
function(data,
         sigs,
         group_col,
         markers,
         normalize = FALSE,
         scale = "none",
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         ...) {

  stopifnot(is.logical(normalize), is.character(scale),
            is.character(gene_id), is.logical(ranks_plot))

  if(length(normalize) == 1)
    normalize <- rep(normalize, length(data))
  if(length(group_col) == 1)
    group_col <- rep(group_col, length(data))
  if(length(slot) == 1)
    slot <- rep(slot, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))

  p <- list()
  for (i in seq_along(data)) {
    if(missing(markers)) {
      p[[i]] <- sig_heatmap(data = data[[i]],
                            sigs = sigs, group_col = group_col[[i]],
                            normalize = normalize[i],
                            scale = scale, gene_id = gene_id[i],
                            slot = slot[i], ranks_plot = ranks_plot, ...)
    } else {
      p[[i]] <- sig_heatmap(data = data[[i]],
                            sigs = sigs, group_col = group_col[[i]],
                            markers = markers, normalize = normalize[i],
                            scale = scale, gene_id = gene_id[i],
                            slot = slot[i], ranks_plot = ranks_plot, ...)
    }

    if(!is.null(names(data))) {
      if(is.ggplot(p[[i]])) {
        p[[i]] <- p[[i]] + ggtitle(names(data)[i])
      } else {
        p[[i]] <- ComplexHeatmap::draw(p[[i]], column_title = names(data)[i]) |>
          grid::grid.grabExpr() |>
          patchwork::wrap_elements()
      }
    } else {
      if(is.ggplot(p[[i]])) {
        p[[i]] <- p[[i]] + ggtitle(i)
      } else {
        p[[i]] <- ComplexHeatmap::draw(p[[i]], column_title = i) |>
          grid::grid.grabExpr() |>
          patchwork::wrap_elements()
      }
    }
  }

  p <- patchwork::wrap_plots(p)
  return(p)
})
