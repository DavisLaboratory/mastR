#' @title Specify cell target_group signature against adherent cancer cell lines.
#'
#' @description Specify signatures against specific tissues or cell lines,
#'   and generate CCLE_tpm_new after retrieving specific tissues.
#'
#' @param data 'CCLE' or expression object
#' @param group_col vector or character, to specify the group of tumor/tissue target_groups
#' @param target_group pattern, specify the target_group of interest, e.g. 'colorectal', match
#'             to primary_disease of CCLE if data is 'CCLE'
#' @param ... params for [grep()] to find matched cell lines in CCLE
#' @param log logical, if to do log transformation on expression
#' @param q num, quantile cutoff to screen markers which are not or lowly
#'          expressed in relevant samples
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                filtered. Default 'NULL' means all genes
#' @param ignore.case logical, if to ignore the case of tissue pattern
#' @param slot character, specify which slot to use only for sce or seurat
#'             object, optional, default 'counts'
#'
#' @return a vector of gene symbols
#'
#' @examples
#' data("NK_markers", "ccle_crc_5")
#' filter_non_tissue(ccle_crc_5, group_col = "cancer", target_group = "CRC",
#'                   markers = NK_markers$HGNC_Symbol)
#'
#' @export
setGeneric("filter_non_tissue",
           function(data = 'CCLE',
                    group_col,
                    target_group,
                    ...,
                    log = FALSE,
                    q = 0.25,
                    markers = NULL,
                    ignore.case = TRUE,
                    slot = "counts")
           standardGeneric("filter_non_tissue"))

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'matrix',
  group_col = 'vector',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  stopifnot(is.logical(log), is.numeric(q), is.logical(ignore.case))

  ## select interested samples
  idx <- grep(target_group, group_col, ignore.case = ignore.case, ...)
  data <- data[,idx]

  ## log-transformation on expr
  if (log == TRUE) {
    data <- log1p(data)
  }

  ## set threshold for expression of non-0 expressed genes to
  ## filter out genes highly expressed by tumor or tissue cells per se
  threshod_25th <- quantile(data[data != 0], q)

  non_tissue_genes <- apply(data, 1, \(x) {
    median(x, na.rm = TRUE) < threshod_25th
  })
  non_tissue_genes <- rownames(data)[non_tissue_genes]

  ## only keep genes which are alo in markers when markers provided
  if(!is.null(markers))
    non_tissue_genes <- intersect(non_tissue_genes, markers)
  non_tissue_genes <- union(setdiff(markers, rownames(data)),
                            non_tissue_genes)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'Matrix',
  group_col = 'vector',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  stopifnot(is.logical(log), is.numeric(q), is.logical(ignore.case))

  ## select interested samples
  idx <- grep(target_group, group_col, ignore.case = ignore.case, ...)
  data <- data[,idx]

  ## log-transformation on expr
  if (log == TRUE) {
    data <- log1p(data)
  }

  ## set threshold for expression of non-0 expressed genes to
  ## filter out genes highly expressed by tumor or tissue cells per se
  threshod_25th <- quantile(data[data != 0], q)

  non_tissue_genes <- apply(data, 1, \(x) {
    median(x, na.rm = TRUE) < threshod_25th
  })
  non_tissue_genes <- rownames(data)[non_tissue_genes]

  ## only keep genes which are alo in markers when markers provided
  if(!is.null(markers))
    non_tissue_genes <- intersect(non_tissue_genes, markers)
  non_tissue_genes <- union(setdiff(markers, rownames(data)),
                            non_tissue_genes)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'data.frame',
  group_col = 'vector',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = as.matrix(data),
                                        group_col = group_col, target_group = target_group,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'character',
  group_col = 'missing',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  stopifnot("data can't be character other than 'CCLE'!" = data == 'CCLE')

  ## load latest version of CCLE TPM from depmap
  CCLE <- depmap::depmap_TPM()
  CCLE_meta <- depmap::depmap_metadata()

  ## convert CCLE from long to wide data format
  CCLE <- CCLE[,c("gene_name", "cell_line", "rna_expression")] |>
    tidyr::pivot_wider(
      names_from = "cell_line",
      values_from = "rna_expression"
      ) |>
    data.frame()
  rownames(CCLE) <- CCLE$gene_name
  CCLE$gene_name <- NULL

  ## match depmap ID of CCLE TPM to CCLE meta to get disease names
  idx <- match(colnames(CCLE), CCLE_meta$cell_line)
  cell_lines <- CCLE_meta$primary_disease[idx]

  non_tissue_genes <- filter_non_tissue(data = CCLE,
                                        group_col = cell_lines,
                                        target_group = target_group,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'missing',
  group_col = 'missing',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = 'CCLE',
                                        target_group = target_group,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'DGEList',
  group_col = 'character',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = data$counts,
                                        group_col = data$samples[[group_col]],
                                        target_group = target_group,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'ExpressionSet',
  group_col = 'character',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = data@assayData$exprs,
                                        group_col = data[[group_col]],
                                        target_group = target_group,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'SummarizedExperiment',
  group_col = 'character',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE,
         slot = "counts") {

  stopifnot(is.character(slot))

  non_tissue_genes <- filter_non_tissue(data = SummarizedExperiment::assay(data, slot),
                                        group_col = SummarizedExperiment::colData(data)[[group_col]],
                                        target_group = target_group,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'Seurat',
  group_col = 'character',
  target_group = 'character'
),
function(data = 'CCLE',
         group_col,
         target_group,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE,
         slot = "counts") {

  stopifnot(is.character(slot))

  non_tissue_genes <- filter_non_tissue(
    data = Seurat::GetAssayData(data, slot = slot),
    group_col = data@meta.data[[group_col]],
    target_group = target_group,
    log = log, q = q,
    markers = markers,
    ignore.case = ignore.case,
    ...
  )

  return(non_tissue_genes)
})
