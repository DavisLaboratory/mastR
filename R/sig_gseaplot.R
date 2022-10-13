#' @include plot.R
#' @import ggplot2 clusterProfiler enrichplot grid ggplotify
NULL

#' Visualize GSEA result with input list of gene symbols.
#'
#' Visualize GSEA result with multiple lists of genes by using `clusterProfiler`.
#'
#' @param data expression data objects, can be DGEList, eSet, Seurat, sce...
#' @param sigs list(s) of gene symbols for GSEA
#' @param ID chr or vector, specify the column name in coldata or factor to compare
#' @param type pattern, specify the group of interest to be compared
#' @param gene_id chr, indicate the ID type of expression data's rowname,
#'                could be one of 'ENSEMBL', 'SYMBOL', ... default 'SYMBOL'
#' @param digits num, specify the number of significant digits of pvalue table
#' @param slot chr, indicate which slot used for seurat/sce object
#' @param ... params for function [get_de_table()]
#'
#' @return patchwork object for all comparisons
#'
#' @examples
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
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, ID = ID, type = type,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       ...)
  gsets <- reshape2::melt(sigs)
  colnames(gsets) <- "SYMBOL"
  gsets$set <- "Signature"
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = T, by.x = "SYMBOL", by.y = "SYMBOL")

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
         ...) {

  stopifnot(is.character(gene_id), is.numeric(digits))

  ## get DEGs tables list with statistics
  tDEG <- get_de_table(data = data, ID = ID, type = type,
                       markers = Reduce(union, sigs),
                       gene_id = gene_id,
                       ...)
  gsets <- reshape2::melt(sigs)
  colnames(gsets) <- c("SYMBOL","set")
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = T, by.x = "SYMBOL", by.y = "SYMBOL")

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
  gsets <- reshape2::melt(sigs)
  colnames(gsets) <- "SYMBOL"
  gsets$set <- "Signature"
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = T, by.x = "SYMBOL", by.y = "SYMBOL")

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
  gsets <- reshape2::melt(sigs)
  colnames(gsets) <- c("SYMBOL","set")
  ids <- AnnotationDbi::select(org.Hs.eg.db, gsets[,1],
                               columns = gene_id,
                               keytype = "SYMBOL")
  gsets <- merge(gsets, ids, all = T, by.x = "SYMBOL", by.y = "SYMBOL")

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
  for (i in 1:length(data)) {
    p[[i]] <- sig_gseaplot(data = data[[i]],
                           sigs = sigs,
                           ID = ID[[i]],
                           type = type[i],
                           gene_id = gene_id[i],
                           digits = digits,
                           slot = slot[i],
                           ...)
    p[[i]] <- ggplotify::as.ggplot(p[[i]])
    if(!is.null(names(data)))
      p[[i]] <- p[[i]] + ggplot2::labs(subtitle = names(data)[i]) +
        ggplot2::theme(plot.subtitle = element_text(hjust = 0.5))
  }
  p <- patchwork::wrap_plots(p) +
    patchwork::plot_layout(guides = "collect")

  return(p)
})

