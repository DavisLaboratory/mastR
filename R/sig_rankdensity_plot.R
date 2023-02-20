#' Plot rank density
#'
#' Show the rank density of given signature in the given comparison.
#'
#' @inheritParams sig_boxplot
#' @param aggregate logical, if to aggregate expression according to `ID`,
#'                  default FALSE
#'
#' @return ggplot or patchwork
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' sig_rankdensity_plot(
#'   data = im_data_6, sigs = NK_markers$HGNC_Symbol[1:10],
#'   ID = "celltype:ch1", gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric("sig_rankdensity_plot",
           function(data,
                    sigs,
                    ID,
                    counts = TRUE,
                    aggregate = FALSE,
                    slot = "counts",
                    gene_id = "SYMBOL")
             standardGeneric("sig_rankdensity_plot"))

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'matrix',
  sigs = 'vector',
  ID = 'vector'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- rankdensity_init(expr = data, sigs = sigs, by = ID, counts = counts,
                        aggregate = aggregate, gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'Matrix',
  sigs = 'vector',
  ID = 'vector'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- rankdensity_init(expr = data, sigs = sigs, by = ID, counts = counts,
                        aggregate = aggregate, gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'data.frame',
  sigs = 'vector',
  ID = 'vector'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- rankdensity_init(expr = as.matrix(data), sigs = sigs, by = ID,
                        counts = counts, aggregate = aggregate,
                        gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'DGEList',
  sigs = 'vector',
  ID = 'character'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- sig_rankdensity_plot(data = data$counts, sigs = sigs,
                            ID = data$samples[[ID]], counts = counts,
                            aggregate = aggregate, gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'ExpressionSet',
  sigs = 'vector',
  ID = 'character'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- sig_rankdensity_plot(data = Biobase::exprs(data), sigs = sigs,
                            ID = data[[ID]], counts = counts,
                            aggregate = aggregate, gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'Seurat',
  sigs = 'vector',
  ID = 'character'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         slot = "counts",
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- sig_rankdensity_plot(data = Seurat::GetAssayData(data, slot = slot),
                            sigs = sigs, ID = data@meta.data[[ID]],
                            counts = counts, aggregate = aggregate,
                            gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'SummarizedExperiment',
  sigs = 'vector',
  ID = 'character'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         slot = "counts",
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  p <- sig_rankdensity_plot(data = SummarizedExperiment::assay(data, slot),
                            sigs = sigs, ID = data@colData[[ID]],
                            counts = counts, aggregate = aggregate,
                            gene_id = gene_id)

  return(p)
})

#' @rdname sig_rankdensity_plot
setMethod("sig_rankdensity_plot", signature(
  data = 'list',
  sigs = 'vector',
  ID = 'character'
),
function(data,
         sigs,
         ID,
         counts = TRUE,
         aggregate = FALSE,
         slot = "counts",
         gene_id = "SYMBOL") {

  stopifnot(is.character(gene_id), is.logical(counts), is.logical(aggregate))

  if(length(counts) == 1)
    counts <- rep(counts, length(data))
  if(length(ID) == 1)
    ID <- rep(ID, length(data))
  if(length(slot) == 1)
    slot <- rep(slot, length(data))
  if(length(aggregate) == 1)
    aggregate <- rep(aggregate, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))

  p <- list()
  for (i in seq_along(data)) {
    p[[i]] <- sig_rankdensity_plot(data = data[[i]], sigs = sigs,
                                   ID = ID[[i]], counts = counts[i],
                                   slot = slot[i],
                                   aggregate = aggregate[i],
                                   gene_id = gene_id[i])

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
