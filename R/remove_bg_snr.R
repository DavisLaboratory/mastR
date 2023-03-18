#' @title Remove markers with background noise.
#'
#' @description Specify signatures against specific tissues or cell lines,
#'   and generate CCLE_tpm_new after retrieving specific tissues.
#'
#' @param bg_data 'CCLE' or expression object as background data
#' @param sig_data expression object, can be matrix or DGEList, as signal data
#' @param b_group_col vector or character, to specify the group of background
#'                     target_groups, or column name of [depmap::depmap_metadata()],
#'                     e.g. 'primary_disease'
#' @param b_target_group pattern, specify the target_group of interest in bg_data,
#'                       e.g. 'colorectal'
#' @param s_group_col vector or character, to specify the group of signal
#'                    target_groups, or column name of group
#' @param s_target_group pattern, specify the target group of interest in sig_data
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                filtered. Must be gene SYMBOLs.
#' @param snr num, the cutoff of SNR to screen markers which are not or lowly
#'            expressed in bg_data
#' @param ... params for [grep()] to find matched cell lines in bg_data
#' @param filter NULL or a vector of 2 num, filter condition to remove low
#'               expression genes in bg_data, the 1st for logcounts,
#'               the 2nd for samples size
#' @param gene_id character, specify the gene ID type of rownames of expression
#'                data, could be one of 'ENSEMBL', 'SYMBOL', 'ENTREZ'...,
#'                default 'SYMBOL'
#' @param b_slot character, specify which slot to use of DGEList, sce or seurat
#'               object for bg_data, optional, default 'counts'
#' @param s_slot character, specify which slot to use of DGEList, sce or seurat
#'               object for sig_data, optional, default 'counts'
#' @param ccle_tpm ccle_tpm data from [depmap::depmap_TPM()], only used when
#'                 data = 'CCLE', default NULL
#' @param ccle_meta ccle_meta data from [depmap::depmap_metadata()], only used
#'                  when data = 'CCLE', default NULL
#'
#' @return a vector of gene symbols
#'
#' @examples
#' data("NK_markers", "ccle_crc_5")
#' remove_bg_snr(ccle_crc_5, Biobase::exprs(im_data_6),
#'               "cancer", "CRC", ## for bg_data
#'               im_data_6$`celltype:ch1`, "NK", ## for sig_data
#'               markers = NK_markers$HGNC_Symbol[30:40],
#'               gene_id = c("SYMBOL", "ENSEMBL"))
#'
#' @export
setGeneric("remove_bg_snr",
           function(bg_data = 'CCLE',
                    sig_data,
                    b_group_col,
                    b_target_group,
                    s_group_col,
                    s_target_group,
                    markers,
                    snr = 1,
                    ...,
                    filter = NULL,
                    gene_id = "SYMBOL",
                    b_slot = "counts",
                    s_slot = "counts",
                    ccle_tpm = NULL,
                    ccle_meta = NULL)
           standardGeneric("remove_bg_snr"))

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'matrix',
  sig_data = 'matrix',
  b_group_col = 'vector',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL") {

  stopifnot("snr must be numeric!" = is.numeric(snr),
            "markers must be vector!" = is.vector(markers))

  ## convert markers to the same ID type
  stopifnot("Please provide a vector of gene symbols!" = is.vector(markers))
  markers <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                   markers,
                                   gene_id, "SYMBOL")

  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, 2)

  ## update markers, only keep markers in sig_data
  m_in <- markers$SYMBOL[markers[[gene_id[2]]] %in% rownames(sig_data)]
  m_out <- setdiff(markers$SYMBOL, m_in)
  if(length(m_out) > 0)
    warning(paste("Gene", m_out, "is not in the sig_data! So remove it!\n"))
  markers <- markers[markers$SYMBOL %in% m_in,]

  ## select target samples
  idx <- grep(b_target_group, b_group_col, ...)
  stopifnot("No matched samples in bg_data for b_target_group!" =
              length(idx) > 0)
  bg_data <- bg_data[,idx]

  idx <- grep(s_target_group, s_group_col, ...)
  stopifnot("No matched samples in sig_data for s_target_group!" =
              length(idx) > 0)
  sig_data <- sig_data[,idx]

  ## filter low expression genes for bg_data
  if(!is.null(filter)) {
    stopifnot("Must provide 2 numbers for filter!" =
                (length(filter) == 2) & is.numeric(filter))
    idx <- (Matrix::rowSums(bg_data > filter[1]) > filter[2])
    bg_data <- bg_data[idx,]
  }

  ## scale data by column
  bg_data <- scale_0_1(bg_data)
  sig_data <- scale_0_1(sig_data)

  ## select markers
  ## bg_data common genes
  idx <- na.omit(match(rownames(bg_data), markers[[gene_id[1]]]))
  bg_data <- bg_data[markers[[gene_id[1]]][idx],]
  m1 <- markers[idx, "SYMBOL"]
  rownames(bg_data) <- m1
  ## get markers without valid expression in bg_data
  m2 <- setdiff(markers$SYMBOL, m1)

  if(length(m1) == 0) {
    message("No gene in markers is in bg_data! Return all markers!")
    return(m2)
  }

  ## sig_data common genes
  idx <- na.omit(match(rownames(sig_data), markers[[gene_id[2]]]))
  sig_data <- sig_data[markers[[gene_id[2]]][idx],]
  rownames(sig_data) <- markers[idx, "SYMBOL"]

  ## compute SNR
  sig_med <- Biobase::rowMedians(sig_data[m1,], na.rm = TRUE)
  bg_med <- Biobase::rowMedians(bg_data[m1,], na.rm = TRUE)
  bg_sd <- sapply(m1, \(m) sd(bg_data[m,], na.rm = TRUE))
  snrs <- (sig_med - bg_med) / bg_sd

  m1 <- m1[snrs > snr]

  final_m <- union(m1, m2)

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'matrix',
  sig_data = 'DGEList',
  b_group_col = 'vector',
  b_target_group = 'character',
  s_group_col = 'character',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         s_slot = "counts") {

  final_m <- remove_bg_snr(bg_data = bg_data,
                           sig_data = sig_data[[s_slot]],
                           b_group_col = b_group_col,
                           b_target_group = b_target_group,
                           s_group_col = sig_data$samples[[s_group_col]],
                           s_target_group = s_target_group,
                           markers = markers,
                           snr = snr,
                           ...,
                           filter = filter,
                           gene_id = gene_id)

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'DGEList',
  sig_data = 'ANY',
  b_group_col = 'character',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         b_slot = "counts",
         s_slot = "counts") {

  stopifnot("sig_data must be matrix or DGEList!" =
              is.matrix(sig_data) | is(sig_data, "DGEList"))

  final_m <- remove_bg_snr(bg_data = bg_data[[b_slot]],
                           sig_data = sig_data,
                           b_group_col = bg_data$samples[[b_group_col]],
                           b_target_group = b_target_group,
                           s_group_col = s_group_col,
                           s_target_group = s_target_group,
                           markers = markers,
                           snr = snr,
                           ...,
                           filter = filter,
                           gene_id = gene_id,
                           s_slot = s_slot)

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'ExpressionSet',
  sig_data = 'ANY',
  b_group_col = 'character',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         s_slot = "counts") {

  stopifnot("sig_data must be matrix or DGEList!" =
              is.matrix(sig_data) | is(sig_data, "DGEList"))

  final_m <- remove_bg_snr(bg_data = Biobase::exprs(bg_data),
                           sig_data = sig_data,
                           b_group_col = Biobase::pData(bg_data)[[b_group_col]],
                           b_target_group = b_target_group,
                           s_group_col = s_group_col,
                           s_target_group = s_target_group,
                           markers = markers,
                           snr = snr,
                           ...,
                           filter = filter,
                           gene_id = gene_id,
                           s_slot = s_slot)

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'SummarizedExperiment',
  sig_data = 'ANY',
  b_group_col = 'character',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         b_slot = "counts",
         s_slot = "counts") {

  stopifnot("sig_data must be matrix or DGEList!" =
              is.matrix(sig_data) | is(sig_data, "DGEList"))

  final_m <- remove_bg_snr(
    bg_data = SummarizedExperiment::assay(bg_data, b_slot),
    sig_data = sig_data,
    b_group_col = SummarizedExperiment::colData(bg_data)[[b_group_col]],
    b_target_group = b_target_group,
    s_group_col = s_group_col,
    s_target_group = s_target_group,
    markers = markers,
    snr = snr,
    ...,
    filter = filter,
    gene_id = gene_id,
    s_slot = s_slot
  )

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'Seurat',
  sig_data = 'ANY',
  b_group_col = 'character',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         b_slot = "counts",
         s_slot = "counts") {

  stopifnot("sig_data must be matrix or DGEList!" =
              is.matrix(sig_data) | is(sig_data, "DGEList"))

  final_m <- remove_bg_snr(
    bg_data = Seurat::GetAssayData(bg_data, slot = b_slot),
    sig_data = sig_data,
    b_group_col = bg_data@meta.data[[b_group_col]],
    b_target_group = b_target_group,
    s_group_col = s_group_col,
    s_target_group = s_target_group,
    markers = markers,
    snr = snr,
    ...,
    filter = filter,
    gene_id = gene_id,
    s_slot = s_slot
  )

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'character',
  sig_data = 'ANY',
  b_group_col = 'vector',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY',
  filter = "ANY"
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         s_slot = "counts",
         ccle_tpm = NULL,
         ccle_meta = NULL) {

  stopifnot("bg_data can't be character other than 'CCLE'!" = bg_data == 'CCLE',
            "filter must be numeric!" = is.numeric(filter))

  ## load latest version of CCLE TPM from depmap
  if(is.null(ccle_tpm))
    ccle_tpm <- depmap::depmap_TPM()
  if(is.null(ccle_meta))
    ccle_meta <- depmap::depmap_metadata()
  ## grep target cell lines
  idx <- grep(b_target_group, ccle_meta[[b_group_col]], ...)
  ccle_meta <- ccle_meta[idx,]
  ccle_tpm <- dplyr::inner_join(ccle_tpm, ccle_meta, by = "depmap_id")

  ccle <- tidyr::pivot_wider(
    ccle_tpm[, c("gene_name", "rna_expression", "depmap_id")],
    names_from = "depmap_id",
    values_from = "rna_expression"
  )
  g_n <- ccle$gene_name
  ccle$gene_name <- NULL
  ccle <- as.matrix(ccle)
  rownames(ccle) <- g_n

  final_m <- remove_bg_snr(bg_data = ccle,
                           sig_data = sig_data,
                           b_group_col = rep(b_target_group, ncol(ccle)),
                           b_target_group = b_target_group,
                           s_group_col = s_group_col,
                           s_target_group = s_target_group,
                           markers = markers,
                           snr = snr,
                           ...,
                           filter = filter,
                           gene_id = gene_id,
                           s_slot = s_slot)

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'missing',
  sig_data = 'ANY',
  b_group_col = 'vector',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY',
  filter = "ANY"
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL",
         s_slot = "counts",
         ccle_tpm = NULL,
         ccle_meta = NULL) {

  final_m <- remove_bg_snr(bg_data = 'CCLE',
                           sig_data = sig_data,
                           b_group_col = b_group_col,
                           b_target_group = b_target_group,
                           s_group_col = s_group_col,
                           s_target_group = s_target_group,
                           markers = markers,
                           snr = snr,
                           ...,
                           filter = filter,
                           gene_id = gene_id,
                           s_slot = s_slot,
                           ccle_tpm = ccle_tpm,
                           ccle_meta = ccle_meta)

  return(final_m)
})

#' @rdname remove_bg_snr
setMethod("remove_bg_snr", signature(
  bg_data = 'ANY',
  sig_data = 'ANY',
  b_group_col = 'vector',
  b_target_group = 'character',
  s_group_col = 'vector',
  s_target_group = 'character',
  markers = 'ANY'
),
function(bg_data = 'CCLE',
         sig_data,
         b_group_col,
         b_target_group,
         s_group_col,
         s_target_group,
         markers,
         snr = 1,
         ...,
         filter = NULL,
         gene_id = "SYMBOL") {

  final_m <- remove_bg_snr(bg_data = as.matrix(bg_data),
                           sig_data = as.matrix(sig_data),
                           b_group_col = b_group_col,
                           b_target_group = b_target_group,
                           s_group_col = s_group_col,
                           s_target_group = s_target_group,
                           markers = markers,
                           snr = snr,
                           ...,
                           filter = filter,
                           gene_id = gene_id)

  return(final_m)
})

#helper: scale data by column to range(0,1)
scale_0_1 <- function(mat) {
  mins <- apply(mat, 2, min)
  maxs <- apply(mat, 2, max)
  mat <- scale(mat, center = mins, scale = maxs - mins)
  return(mat)
}


utils::globalVariables(c("gene_name", "rna_expression"))
