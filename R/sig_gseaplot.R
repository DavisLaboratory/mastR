#' @include plot.R
NULL

#' Visualize GSEA result with input list of gene symbols.
#'
#' Visualize GSEA result with multiple lists of genes by using `clusterProfiler`.
#'
#' @inheritParams sig_boxplot
#' @param digits num, specify the number of significant digits of pvalue table
#' @param method one of "gseaplot" and "dotplot", how to plot GSEA result
#' @param size column name of [clusterProfiler::GSEA()] result, used for dot
#'             size when method = "dotplot"
#' @param ... params for function [get_de_table()]
#'
#' @return patchwork object for all comparisons
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' sig_gseaplot(sigs = list(A = NK_markers$HGNC_Symbol[1:15],
#'                          B = NK_markers$HGNC_Symbol[20:40],
#'                          C = NK_markers$HGNC_Symbol[60:75]),
#'              data = im_data_6, group_col = "celltype:ch1",
#'              target_group = "NK", gene_id = "ENSEMBL")
#'
#' @export
setGeneric("sig_gseaplot",
           function(data,
                    sigs,
                    group_col,
                    target_group,
                    gene_id = "SYMBOL",
                    digits = 2,
                    slot = "counts",
                    method = c("dotplot", "gseaplot"),
                    size = "enrichmentScore",
                    ...)
             standardGeneric("sig_gseaplot"))

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'ANY',
  sigs = 'vector',
  group_col = 'ANY',
  target_group = 'ANY'
),
function(data,
         sigs,
         group_col,
         target_group,
         gene_id = "SYMBOL",
         digits = 2,
         slot = "counts",
         method = c("dotplot", "gseaplot"),
         size = "enrichmentScore",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))
  method <- match.arg(method)

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, group_col = group_col,
                       target_group = target_group,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       slot = slot,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables
  gsets <- data.frame(SYMBOL = sigs, set = "Signature")
  ids <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  ## gsea
  gse <- gsea_analysis(tDEG = tDEG, gsets = gsets,
                       gene_id = gene_id, digits = digits)

  if(method == "gseaplot") {
    p <- gsea_plot_init(gse)
  } else p <- gsea_dotplot_init(gse, size = size)

  p <-  patchwork::wrap_plots(p) + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'ANY',
  sigs = 'list',
  group_col = 'ANY',
  target_group = 'ANY'
),
function(data,
         sigs,
         group_col,
         target_group,
         gene_id = "SYMBOL",
         digits = 2,
         slot = "counts",
         method = c("dotplot", "gseaplot"),
         size = "enrichmentScore",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))
  method <- match.arg(method)

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, group_col = group_col,
                       target_group = target_group,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       slot = slot,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables

  if(is.null(names(sigs)))
    names(sigs) <- seq_along(sigs)  ## set gene list names
  gsets <- utils::stack(sigs)
  colnames(gsets) <- c("SYMBOL", "set")
  ids <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  ## gsea
  gse <- gsea_analysis(tDEG = tDEG, gsets = gsets,
                       gene_id = gene_id, digits = digits)

  if(method == "gseaplot") {
    p <- gsea_plot_init(gse)
  } else p <- gsea_dotplot_init(gse, size = size)

  p <- patchwork::wrap_plots(p) + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'DGEList',
  sigs = 'vector',
  group_col = 'ANY',
  target_group = 'ANY'
),
function(data,
         sigs,
         group_col,
         target_group,
         gene_id = "SYMBOL",
         digits = 2,
         method = c("dotplot", "gseaplot"),
         size = "enrichmentScore",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))
  method <- match.arg(method)

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, group_col = group_col,
                       target_group = target_group,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables
  gsets <- data.frame(SYMBOL = sigs, set = "Signature")
  ids <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  ## gsea
  gse <- gsea_analysis(tDEG = tDEG, gsets = gsets,
                       gene_id = gene_id, digits = digits)

  if(method == "gseaplot") {
    p <- gsea_plot_init(gse)
  } else p <- gsea_dotplot_init(gse, size = size)

  p <- patchwork::wrap_plots(p) + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'DGEList',
  sigs = 'list',
  group_col = 'ANY',
  target_group = 'ANY'
),
function(data,
         sigs,
         group_col,
         target_group,
         gene_id = "SYMBOL",
         digits = 2,
         method = c("dotplot", "gseaplot"),
         size = "enrichmentScore",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))
  method <- match.arg(method)

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, group_col = group_col,
                       target_group = target_group,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables

  if(is.null(names(sigs)))
    names(sigs) <- seq_along(sigs)  ## set gene list names
  gsets <- utils::stack(sigs)
  colnames(gsets) <- c("SYMBOL", "set")
  ids <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  ## gsea
  gse <- gsea_analysis(tDEG = tDEG, gsets = gsets,
                       gene_id = gene_id, digits = digits)

  if(method == "gseaplot") {
    p <- gsea_plot_init(gse)
  } else p <- gsea_dotplot_init(gse, size = size)

  p <- patchwork::wrap_plots(p) + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'list',
  sigs = 'ANY',
  group_col = 'ANY',
  target_group = 'ANY'
),
function(data,
         sigs,
         group_col,
         target_group,
         gene_id = "SYMBOL",
         digits = 2,
         slot = "counts",
         method = c("dotplot", "gseaplot"),
         size = "enrichmentScore",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))
  method <- match.arg(method)

  if(length(target_group) == 1)
    target_group <- rep(target_group, length(data))
  if(length(group_col) == 1)
    group_col <- rep(group_col, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))
  if(length(slot) == 1)
    slot <- rep(slot, length(data))

  p <- list()
  for (i in seq_along(data)) {
    p[[i]] <- sig_gseaplot(data = data[[i]],
                           sigs = sigs,
                           group_col = group_col[[i]],
                           target_group = target_group[i],
                           gene_id = gene_id[i],
                           digits = digits,
                           slot = slot[i],
                           method = method,
                           size = size,
                           ...)
    p[[i]] <- patchwork::patchworkGrob(p[[i]]) |> ggpubr::as_ggplot()
    if(!is.null(names(data)))
      p[[i]] <- p[[i]] + ggplot2::labs(subtitle = names(data)[i]) +
        ggplot2::theme(plot.subtitle = element_text(hjust = 0.5))
  }
  p <- patchwork::wrap_plots(p) +
    patchwork::plot_layout(guides = "collect")

  return(p)
})
