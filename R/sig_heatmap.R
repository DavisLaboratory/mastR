#' @include plot.R
#' @import ggplot2 patchwork
NULL

#' Heatmap original markers and screened signature
#'
#' Compare the heatmap before and after screening.
#'
#' @inheritParams sig_boxplot
#' @param markers a vector of gene names, listed the gene symbols of original
#'                markers pool
#' @param scale could be 'none', 'row' or 'column'
#' @param min_max logical, if to use min_max to normalize data by rows
#' @param ranks_plot logical, if to use ranks instead of expression of genes to
#'                   draw heatmap
#' @param col vactor of color to draw heatmap
#'
#' @return patchwork object of pheatmaps
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' sig_heatmap(data = im_data_6, sigs = NK_markers$HGNC_Symbol[1:10],
#'             ID = "celltype:ch1", markers = NK_markers$HGNC_Symbol,
#'             gene_id = "ENSEMBL")
#'
#' @export
setGeneric("sig_heatmap",
           function(data,
                    sigs,
                    ID,
                    markers,
                    counts = TRUE,
                    scale = "none",
                    min_max = FALSE,
                    gene_id = "SYMBOL",
                    ranks_plot = FALSE,
                    slot = "counts",
                    col = colorRampPalette(c("#76B7B2", "#E15759"))(256))
             standardGeneric("sig_heatmap"))

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'matrix',
  sigs = 'vector',
  ID = 'vector',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  stopifnot(is.logical(counts), is.character(scale), is.logical(min_max),
            is.character(gene_id), is.logical(ranks_plot))

  p <- heatmap_init(expr = data, sigs = sigs, by = ID, markers = markers,
                    counts = counts, scale = scale, min_max = min_max,
                    gene_id = gene_id, ranks_plot = ranks_plot, col = col)

  p <- ((p[[1]]) + p[[2]]) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Heatmaps of logCPM",
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'Matrix',
  sigs = 'vector',
  ID = 'vector',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  stopifnot(is.logical(counts), is.character(scale), is.logical(min_max),
            is.character(gene_id), is.logical(ranks_plot))

  p <- heatmap_init(expr = data, sigs = sigs, by = ID, markers = markers,
                    counts = counts, scale = scale, min_max = min_max,
                    gene_id = gene_id, ranks_plot = ranks_plot, col = col)

  p <- (p[[1]] + p[[2]]) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Heatmaps of logCPM",
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'data.frame',
  sigs = 'vector',
  ID = 'vector',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  p <- sig_heatmap(data = as.matrix(data), sigs = sigs, ID = ID,
                   markers = markers, counts = counts, scale = scale,
                   min_max = min_max, gene_id = gene_id,
                   ranks_plot = ranks_plot, col = col)
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'DGEList',
  sigs = 'vector',
  ID = 'character',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  p <- sig_heatmap(data = data$counts, sigs = sigs, ID = data$samples[[ID]],
                   markers = markers, counts = counts, scale = scale,
                   min_max = min_max, gene_id = gene_id,
                   ranks_plot = ranks_plot, col = col)
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'ExpressionSet',
  sigs = 'vector',
  ID = 'character',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  p <- sig_heatmap(data = Biobase::exprs(data), sigs = sigs,
                   ID = data[[ID]], markers = markers, counts = counts,
                   scale = scale, min_max = min_max, gene_id = gene_id,
                   ranks_plot = ranks_plot, col = col)
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'Seurat',
  sigs = 'vector',
  ID = 'character',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  p <- sig_heatmap(data = Seurat::GetAssayData(data, slot = slot),
                   sigs = sigs, ID = data@meta.data[[ID]],
                   markers = markers, counts = counts,
                   scale = scale, min_max = min_max, gene_id = gene_id,
                   ranks_plot = ranks_plot, col = col)
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'SummarizedExperiment',
  sigs = 'vector',
  ID = 'character',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  p <- sig_heatmap(data = SummarizedExperiment::assay(data, slot),
                   sigs = sigs, ID = data@colData[[ID]],
                   markers = markers, counts = counts,
                   scale = scale, min_max = min_max, gene_id = gene_id,
                   ranks_plot = ranks_plot, col = col)
  return(p)
})

#' @rdname sig_heatmap
setMethod("sig_heatmap", signature(
  data = 'list',
  sigs = 'vector',
  ID = 'character',
  markers = 'vector'
),
function(data,
         sigs,
         ID,
         markers,
         counts = TRUE,
         scale = "none",
         min_max = FALSE,
         gene_id = "SYMBOL",
         ranks_plot = FALSE,
         slot = "counts",
         col = colorRampPalette(c("#76B7B2", "#E15759"))(256)) {

  stopifnot(is.logical(counts), is.character(scale), is.logical(min_max),
            is.character(gene_id), is.logical(ranks_plot))

  if(length(counts) == 1)
    counts <- rep(counts, length(data))
  if(length(ID) == 1)
    ID <- rep(ID, length(data))
  if(length(slot) == 1)
    slot <- rep(slot, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))

  p <- list()
  for (i in seq_along(data)) {
    p[[i]] <- sig_heatmap(data = data[[i]],
                          sigs = sigs, ID = ID[[i]],
                          markers = markers, counts = counts[i],
                          scale = scale, min_max = min_max,
                          gene_id = gene_id[i], slot = slot[i],
                          ranks_plot = ranks_plot, col = col)
    if(!is.null(names(data)))
      p[[i]] <- p[[i]] +
        patchwork::plot_annotation(
          subtitle = names(data)[i],
          theme = theme(plot.subtitle = element_text(hjust = 0.5))
        )
    p[[i]] <- patchwork::patchworkGrob(p[[i]]) |> ggpubr::as_ggplot()
  }

  p <- patchwork::wrap_plots(p) +
    patchwork::plot_layout(guides = "collect")
  return(p)
})
