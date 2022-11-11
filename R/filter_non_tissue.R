#' @title Specify cell type signature against adherent cancer cell lines.
#'
#' @description Specify signatures against specific tissues or cell lines,
#'   and generate CCLE_tpm_new after retrieving specific tissues.
#'
#' @param data 'CCLE' or expression object
#' @param ID vector or character, to specify the group of tumor/tissue types
#' @param type pattern, specify the type of interest, e.g. 'colorectal', match
#'             to primary_disease of CCLE if data is 'CCLE'
#' @param ... params for [grep()] to find matched cell lines in CCLE
#' @param log logical, if to do log transformation on expression
#' @param q num, quantile cutoff to screen markers which are not or lowly
#'          expressed in relevant samples
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                filtered. Default 'NULL' means all genes in data to be filtered
#' @param ignore.case logical, if to ignore the case of tissue pattern
#' @param slot character, specify which slot to use only for sce or seurat object,
#'             optional, default 'counts'
#'
#' @return a vector of gene symbols
#'
#' @examples
#' filter_non_tissue(type = "colorectal", markers = NK_markers$HGNC_Symbol)
#'
#' @export
setGeneric("filter_non_tissue",
           function(data = 'CCLE',
                    ID,
                    type,
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
  ID = 'vector',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  stopifnot(is.logical(log), is.numeric(q), is.logical(ignore.case))

  ## select interested samples
  idx <- grep(type, ID, ignore.case = ignore.case, ...)
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
  ID = 'vector',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  stopifnot(is.logical(log), is.numeric(q), is.logical(ignore.case))

  ## select interested samples
  idx <- grep(type, ID, ignore.case = ignore.case, ...)
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
  ID = 'vector',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = as.matrix(data),
                                        ID = ID, type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'character',
  ID = 'missing',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
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
  CCLE <- tidyr::pivot_wider(CCLE[,c("gene_name", "cell_line", "rna_expression")],
                             names_from = "cell_line",
                             values_from = "rna_expression") |> data.frame()
  rownames(CCLE) <- CCLE$gene_name
  CCLE$gene_name <- NULL

  ## match depmap ID of CCLE TPM to CCLE meta to get disease names
  idx <- match(colnames(CCLE), CCLE_meta$cell_line)
  cell_lines <- CCLE_meta$primary_disease[idx]

  non_tissue_genes <- filter_non_tissue(data = CCLE,
                                        ID = cell_lines,
                                        type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'missing',
  ID = 'missing',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = 'CCLE',
                                        type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'DGEList',
  ID = 'character',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = data$counts,
                                        ID = data$samples[[ID]],
                                        type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'ExpressionSet',
  ID = 'character',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE) {

  non_tissue_genes <- filter_non_tissue(data = data@assayData$exprs,
                                        ID = data[[ID]],
                                        type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'SummarizedExperiment',
  ID = 'character',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE,
         slot = "counts") {

  stopifnot(is.character(slot))

  non_tissue_genes <- filter_non_tissue(data = SummarizedExperiment::assay(data, slot),
                                        ID = SummarizedExperiment::colData(data)[[ID]],
                                        type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

#' @rdname filter_non_tissue
setMethod("filter_non_tissue", signature(
  data = 'Seurat',
  ID = 'character',
  type = 'character'
),
function(data = 'CCLE',
         ID,
         type,
         ...,
         log = FALSE,
         q = 0.25,
         markers = NULL,
         ignore.case = TRUE,
         slot = "counts") {

  stopifnot(is.character(slot))

  non_tissue_genes <- filter_non_tissue(data = Seurat::GetAssayData(data, slot = slot),
                                        ID = data@meta.data[[ID]],
                                        type = type,
                                        log = log, q = q,
                                        markers = markers,
                                        ignore.case = ignore.case,
                                        ...)

  return(non_tissue_genes)
})

