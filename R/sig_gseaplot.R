#' @include plot.R
NULL

#' Visualize GSEA result with input list of gene symbols.
#'
#' Visualize GSEA result with multiple lists of genes by using `clusterProfiler`.
#'
#' @inheritParams sig_boxplot
#' @param digits num, specify the number of significant digits of pvalue table
#' @param ... params for function [get_de_table()]
#'
#' @return patchwork object for all comparisons
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' sig_gseaplot(sigs = list(A = NK_markers$HGNC_Symbol[1:15],
#'                          B = NK_markers$HGNC_Symbol[20:40],
#'                          C = NK_markers$HGNC_Symbol[60:75]),
#'              data = im_data_6, ID = "celltype:ch1",
#'              type = "NK", gene_id = "ENSEMBL")
#'
#' @export
setGeneric("sig_gseaplot",
           function(data,
                    sigs,
                    ID,
                    type,
                    gene_id = "SYMBOL",
                    digits = 2,
                    slot = "counts",
                    ...)
             standardGeneric("sig_gseaplot"))

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'ANY',
  sigs = 'vector',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         sigs,
         ID,
         type,
         gene_id = "SYMBOL",
         digits = 2,
         slot = "counts",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, ID = ID, type = type,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       slot = slot,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables
  gsets <- data.frame(SYMBOL = sigs, set = "Signature")
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  p <- gsea_plot_init(tDEG = tDEG, gsets = gsets,
                      gene_id = gene_id, digits = digits) |>
    patchwork::wrap_plots() + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'ANY',
  sigs = 'list',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         sigs,
         ID,
         type,
         gene_id = "SYMBOL",
         digits = 2,
         slot = "counts",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, ID = ID, type = type,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       slot = slot,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables

  if(is.null(names(sigs)))
    names(sigs) <- seq_along(sigs)  ## set gene list names
  gsets <- utils::stack(sigs)
  colnames(gsets) <- c("SYMBOL", "set")
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  p <- gsea_plot_init(tDEG = tDEG, gsets = gsets,
                      gene_id = gene_id, digits = digits) |>
    patchwork::wrap_plots() + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'DGEList',
  sigs = 'vector',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         sigs,
         ID,
         type,
         gene_id = "SYMBOL",
         digits = 2,
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, ID = ID, type = type,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables
  gsets <- data.frame(SYMBOL = sigs, set = "Signature")
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  p <- gsea_plot_init(tDEG = tDEG, gsets = gsets,
                      gene_id = gene_id, digits = digits) |>
    patchwork::wrap_plots() + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'DGEList',
  sigs = 'list',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         sigs,
         ID,
         type,
         gene_id = "SYMBOL",
         digits = 2,
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, ID = ID, type = type,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       ...)

  tDEG <- tDEG[which(names(tDEG) != "proc_data")]  ## only keep DEG tables

  if(is.null(names(sigs)))
    names(sigs) <- seq_along(sigs)  ## set gene list names
  gsets <- utils::stack(sigs)
  colnames(gsets) <- c("SYMBOL", "set")
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = TRUE, by.x = "SYMBOL", by.y = "SYMBOL")

  p <- gsea_plot_init(tDEG = tDEG, gsets = gsets,
                      gene_id = gene_id, digits = digits) |>
    patchwork::wrap_plots() + patchwork::plot_layout(guides = "collect")
  return(p)
})

#' @rdname sig_gseaplot
setMethod("sig_gseaplot", signature(
  data = 'list',
  sigs = 'ANY',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         sigs,
         ID,
         type,
         gene_id = "SYMBOL",
         digits = 2,
         slot = "counts",
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  if(length(type) == 1)
    type <- rep(type, length(data))
  if(length(ID) == 1)
    ID <- rep(ID, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))
  if(length(slot) == 1)
    slot <- rep(slot, length(data))

  p <- list()
  for (i in seq_along(data)) {
    p[[i]] <- sig_gseaplot(data = data[[i]],
                           sigs = sigs,
                           ID = ID[[i]],
                           type = type[i],
                           gene_id = gene_id[i],
                           digits = digits,
                           slot = slot[i],
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

utils::globalVariables(c("org.Hs.eg.db"))
