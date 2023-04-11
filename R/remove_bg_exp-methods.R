#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "matrix",
    bg_data = "matrix",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL") {
    stopifnot(
      "s_group_col must be a vector!" =
        is.vector(s_group_col) | is.null(s_group_col),
      "s_target_group must be a one-length regex pattern!" =
        (is.character(s_target_group) & length(s_target_group) == 1) |
          is.null(s_target_group),
      "b_group_col must be a vector!" =
        is.vector(b_group_col) | is.null(b_group_col),
      "b_target_group must be a one-length regex pattern!" =
        (is.character(b_target_group) & length(b_target_group) == 1) |
          is.null(b_target_group),
      "snr must be numeric!" = is.numeric(snr)
    )

    ## select target samples in sig_data
    if (!is.null(s_group_col) & !is.null(s_target_group)) {
      idx <- grep(s_target_group, s_group_col, ...)
      stopifnot(
        "No matched samples in sig_data for s_target_group!" =
          length(idx) > 0
      )
      sig_data <- sig_data[, idx]
    } else {
      message("All samples in sig_data are used!")
    }

    ## select target samples in bg_data
    if (!is.null(b_group_col) & !is.null(b_target_group)) {
      idx <- grep(b_target_group, b_group_col, ...)
      stopifnot(
        "No matched samples in bg_data for b_target_group!" =
          length(idx) > 0
      )
      bg_data <- bg_data[, idx]
    } else {
      message("All samples in bg_data are used!")
    }

    ## filter low expression genes for bg_data
    if (!is.null(filter)) {
      stopifnot(
        "Must provide 2 numbers for filter!" =
          (length(filter) == 2) & is.numeric(filter)
      )
      idx <- (Matrix::rowSums(bg_data > filter[1], na.rm = TRUE) > filter[2])
      stopifnot("No gene kept after filtration!" = sum(idx, na.rm = TRUE) > 0)
      bg_data <- bg_data[idx, ]
    }

    final_m <- remove_bg_exp_mat(
      sig_mat = sig_data,
      bg_mat = bg_data,
      markers = markers,
      snr = snr,
      gene_id = gene_id
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "DGEList",
    bg_data = "matrix",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts") {
    stopifnot(
      "s_group_col must be a column name!" =
        is.character(s_group_col) | is.null(s_group_col)
    )

    if (!is.null(s_group_col)) {
      s_group_col <- sig_data$samples[[s_group_col]]
    }

    final_m <- remove_bg_exp(
      bg_data = bg_data,
      sig_data = sig_data[[s_slot]],
      b_group_col = b_group_col,
      b_target_group = b_target_group,
      s_group_col = s_group_col,
      s_target_group = s_target_group,
      markers = markers,
      snr = snr,
      ...,
      filter = filter,
      gene_id = gene_id
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "DGEList",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts",
           b_slot = "counts") {
    stopifnot(
      "b_group_col must be a column name!" =
        is.character(b_group_col) | is.null(b_group_col),
      "sig_data must be matrix or DGEList!" =
        is.matrix(sig_data) | is(sig_data, "DGEList")
    )

    if (!is.null(b_group_col)) {
      b_group_col <- bg_data$samples[[b_group_col]]
    }

    final_m <- remove_bg_exp(
      bg_data = bg_data[[b_slot]],
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
      s_slot = s_slot
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "ExpressionSet",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts") {
    stopifnot(
      "b_group_col must be a column name!" =
        is.character(b_group_col) | is.null(b_group_col),
      "sig_data must be matrix or DGEList!" =
        is.matrix(sig_data) | is(sig_data, "DGEList")
    )

    if (!is.null(b_group_col)) {
      b_group_col <- Biobase::pData(bg_data)[[b_group_col]]
    }

    final_m <- remove_bg_exp(
      bg_data = Biobase::exprs(bg_data),
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
      s_slot = s_slot
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "SummarizedExperiment",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts",
           b_slot = "counts") {
    stopifnot(
      "b_group_col must be a column name!" =
        is.character(b_group_col) | is.null(b_group_col),
      "sig_data must be matrix or DGEList!" =
        is.matrix(sig_data) | is(sig_data, "DGEList")
    )

    if (!is.null(b_group_col)) {
      b_group_col <- SummarizedExperiment::colData(bg_data)[[b_group_col]]
    }

    final_m <- remove_bg_exp(
      bg_data = SummarizedExperiment::assay(bg_data, b_slot),
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
      s_slot = s_slot
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "Seurat",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts",
           b_slot = "counts") {
    stopifnot(
      "b_group_col must be a column name!" =
        is.character(b_group_col) | is.null(b_group_col),
      "sig_data must be matrix or DGEList!" =
        is.matrix(sig_data) | is(sig_data, "DGEList")
    )

    if (!is.null(b_group_col)) {
      b_group_col <- slot(bg_data, "meta.data")[[b_group_col]]
    }

    final_m <- remove_bg_exp(
      bg_data = SeuratObject::GetAssayData(bg_data, slot = b_slot),
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
      s_slot = s_slot
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "character",
    markers = "vector",
    filter = "ANY"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts",
           ccle_tpm = NULL,
           ccle_meta = NULL) {
    stopifnot(
      "bg_data can't be character other than 'CCLE'!" = bg_data == "CCLE",
      "b_group_col must be a column name!" =
        is.character(b_group_col) | is.null(b_group_col),
      "Must provide 2 numbers for filter!" =
        (length(filter) == 2) & is.numeric(filter)
    )

    ## load latest version of CCLE TPM from depmap
    if (is.null(ccle_tpm)) {
      ccle_tpm <- depmap::depmap_TPM()
    }
    if (is.null(ccle_meta)) {
      ccle_meta <- depmap::depmap_metadata()
    }
    ## grep target cell lines
    if (!is.null(b_group_col) & !is.null(b_target_group)) {
      idx <- grep(b_target_group, ccle_meta[[b_group_col]], ...)
      ccle_meta <- ccle_meta[idx, ]
    } else {
      message("All cell lines in CCLE are used!")
    }
    ccle_tpm <- dplyr::inner_join(ccle_tpm, ccle_meta, by = "depmap_id")

    ## convert long data into wide matrix
    ccle_tpm <- ccle_2_wide(ccle = ccle_tpm)

    final_m <- remove_bg_exp(
      bg_data = ccle_tpm,
      sig_data = sig_data,
      b_group_col = rep(b_target_group, ncol(ccle_tpm)),
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
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "missing",
    markers = "vector",
    filter = "ANY"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts",
           ccle_tpm = NULL,
           ccle_meta = NULL) {
    final_m <- remove_bg_exp(
      bg_data = "CCLE",
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
      ccle_meta = ccle_meta
    )

    return(final_m)
  }
)

#' @rdname remove_bg_exp
setMethod(
  "remove_bg_exp", signature(
    sig_data = "ANY",
    bg_data = "ANY",
    markers = "vector"
  ),
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL") {
    final_m <- remove_bg_exp(
      bg_data = as.matrix(bg_data),
      sig_data = as.matrix(sig_data),
      b_group_col = b_group_col,
      b_target_group = b_target_group,
      s_group_col = s_group_col,
      s_target_group = s_target_group,
      markers = markers,
      snr = snr,
      ...,
      filter = filter,
      gene_id = gene_id
    )

    return(final_m)
  }
)

#' Remove genes show high signal in the background expression data from markers.
#'
#' @param sig_mat expression matrix of interested signal data
#' @param bg_mat expression matrix of interested background data
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                filtered. Must be gene SYMBOLs.
#' @param snr num, the cutoff of SNR to screen markers which are not or lowly
#'            expressed in bg_data
#' @param gene_id character, specify the gene ID types of row names of sig_mat
#'                and bg_mat data, could be one of 'ENSEMBL', 'SYMBOL', 'ENTREZ'...,
#'                default 'SYMBOL'
#'
#' @return a vector of genes after filtration
#' @export
#'
#' @examples
#' data("im_data_6", "nk_markers", "ccle_crc_5")
#' remove_bg_exp_mat(
#'   sig_mat = Biobase::exprs(im_data_6),
#'   bg_mat = ccle_crc_5$counts,
#'   markers = nk_markers$HGNC_Symbol[30:40],
#'   gene_id = c("ENSEMBL", "SYMBOL")
#' )
remove_bg_exp_mat <- function(
    sig_mat,
    bg_mat,
    markers,
    snr = 1,
    gene_id = "SYMBOL") {
  stopifnot(
    "sig_mat must be a matrix!" = is.matrix(sig_mat),
    "bg_mat must be a matrix!" = is.matrix(bg_mat),
    "markers must be a vector of genes!" = is.vector(markers),
    "snr must be a number!" = is.numeric(snr)
  )

  ## convert markers to the same ID type
  if (length(union(gene_id, "SYMBOL")) == 1) {
    markers <- data.frame(SYMBOL = markers)
  } else {
    markers <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      markers,
      gene_id, "SYMBOL"
    )
  }

  if (length(gene_id) == 1) {
    gene_id <- rep(gene_id, 2)
  }

  ## update markers, only keep markers in sig_mat
  m_in <- markers$SYMBOL[markers[[gene_id[1]]] %in% rownames(sig_mat)]
  stopifnot("No gene in markers is in sig_mat!" = length(m_in) > 0)
  m_out <- setdiff(markers$SYMBOL, m_in)
  if (length(m_out) > 0) {
    warning(paste("Gene", m_out, "is not in the sig_mat! So remove it!\n"))
  }
  markers <- subset(markers, SYMBOL %in% m_in)

  ## scale data by column
  bg_mat <- scale_0_1(bg_mat)
  sig_mat <- scale_0_1(sig_mat)

  ## select markers
  ## bg_mat common genes
  idx <- na.omit(match(rownames(bg_mat), markers[[gene_id[2]]]))
  if (length(idx) == 0) {
    warning("No gene in markers is in bg_mat! Return all markers!")
    return(unique(markers$SYMBOL))
  }
  bg_mat <- bg_mat[markers[[gene_id[2]]][idx], ]
  if (length(idx) == 1) {
    bg_mat <- matrix(bg_mat, nrow = length(idx))
  }
  m1 <- markers[idx, "SYMBOL"]
  rownames(bg_mat) <- m1
  ## get markers without valid expression in bg_mat
  m2 <- setdiff(markers$SYMBOL, m1)

  ## sig_mat common genes
  idx <- na.omit(match(rownames(sig_mat), markers[[gene_id[1]]]))
  sig_mat <- sig_mat[markers[[gene_id[1]]][idx], ]
  if (length(idx) == 1) {
    sig_mat <- matrix(sig_mat, nrow = length(idx))
  }
  rownames(sig_mat) <- markers[idx, "SYMBOL"]

  ## compute SNR
  if (length(m1) == 1) {
    sig_med <- median(sig_mat[m1, ], na.rm = TRUE)
    bg_med <- median(bg_mat[m1, ], na.rm = TRUE)
    bg_sd <- sd(bg_mat[m1, ], na.rm = TRUE)
    snrs <- (sig_med - bg_med) / bg_sd
  } else {
    sig_med <- Biobase::rowMedians(sig_mat[m1, ], na.rm = TRUE)
    bg_med <- Biobase::rowMedians(bg_mat[m1, ], na.rm = TRUE)
    bg_sd <- apply(bg_mat[m1, ], 1, sd, na.rm = TRUE)
    snrs <- (sig_med - bg_med) / bg_sd
  }

  m1 <- m1[snrs > snr]

  final_m <- union(m1, m2)

  return(final_m)
}

#' @title Convert CCLE data from long data to wide data.
#'
#' @description Convert CCLE data downloaded by [depmap::depmap_TPM()] from long
#'     data into wide matrix, with row names are gene names and column names are
#'     depmap IDs.
#' @param ccle CCLE data downloaded by [depmap::depmap_TPM()]
#'
#' @return a matrix
#' @export
#'
#' @examples
#' data("ccle_crc_5")
#' ccle <- data.frame(
#'   gene_name = rownames(ccle_crc_5),
#'   ccle_crc_5$counts
#' ) |>
#'   tidyr::pivot_longer(
#'     -gene_name,
#'     names_to = "depmap_id",
#'     values_to = "rna_expression"
#'   )
#' ccle_wide <- ccle_2_wide(ccle)
ccle_2_wide <- function(ccle) {
  ccle <- tidyr::pivot_wider(
    ccle[, c("gene_name", "rna_expression", "depmap_id")],
    names_from = "depmap_id",
    values_from = "rna_expression"
  )
  g_n <- ccle$gene_name
  ccle$gene_name <- NULL
  ccle <- as.matrix(ccle)
  rownames(ccle) <- g_n
  return(ccle)
}

# helper: scale data by column to range(0,1)
scale_0_1 <- function(mat) {
  mins <- apply(mat, 2, min)
  maxs <- apply(mat, 2, max)
  mat <- scale(mat, center = mins, scale = maxs - mins)
  return(mat)
}


utils::globalVariables(c("gene_name", "rna_expression", "SYMBOL"))
