#' Filter specific cell type signature genes against other subsets.
#'
#' Specify the signature of the subset matched 'type' against other subsets,
#' and generate proc_data after removing low expression genes and limma::voom
#' transformation.
#'
#' @inheritParams de_analysis
#' @param data A list of expression data objects
#' @param ID vec or chr, specify the group factor or column name of coldata for
#'           DE comparisons
#' @param dir chr, could be 'UP' or 'DOWN' to use only up- or down-expressed
#'            genes
#' @param comb 'RRA' or Fun, keep all passing genes or only intersected genes,
#'              could be `union` or `intersect` or `setdiff` or customed Fun,
#'              or could be 'RRA' to use Robust Rank Aggregation method for
#'              integrating multi-lists if DEGs, default 'union'
#' @param filter (list of) vector of 2 numbers, filter condition to remove low
#'               expression genes, the 1st for min.counts (if counts = TRUE) or
#'               CPM/TPM (if counts = F), the 2nd for samples size 'large.n'
#' @param s_thres num, threshold of score if comb = 'RRA'
#' @param ... other params for [get_degs()]
#'
#' @return a vector of gene symbols
#'
#' @examples
#' sigs <- filter_subset_sig(im_data_6, "celltype:ch1", "NK",
#'                           markers = NK_markers$HGNC_Symbol,
#'                           gene_id = "ENSEMBL")
#'
#' @export
setGeneric("filter_subset_sig",
           function(data,
                    markers,
                    ID,
                    type,
                    counts = TRUE,
                    dir = "UP",
                    gene_id = "SYMBOL",
                    method = "RP",
                    comb = union,
                    filter = c(10,10),
                    s_thres = 0.05,
                    ...)
           standardGeneric("filter_subset_sig"))

#' @rdname filter_subset_sig
setMethod("filter_subset_sig", signature(
  data = 'list',
  markers = 'vector',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         markers,
         ID,
         type,
         counts = TRUE,
         dir = "UP",
         gene_id = "SYMBOL",
         method = "RP",
         comb = union,
         filter = c(10,10),
         s_thres = 0.05,
         ...) {

  stopifnot(is.vector(markers), is.character(dir),
            is.character(method), is.numeric(s_thres),
            is.character(gene_id))

  markers <- AnnotationDbi::select(org.Hs.eg.db, markers, gene_id, "SYMBOL")
  # tryCatch(
  #   {
  #     markers$ENSEMBL[match(c("KLRA1P", "TRBC1", "TRDC"), markers$SYMBOL)] <- c(
  #       "ENSG00000256667",
  #       "ENSG00000211751",
  #       "ENSG00000211829"
  #     )
  #   },
  #   error = function(e) {
  #     "not containing KLRA1P,TRBC1,TRDC"
  #   }
  # )

  if(length(ID) == 1)
    ID <- rep(ID, length(data))
  if(length(type) == 1)
    type <- rep(type, length(data))
  if(length(counts) == 1)
    counts <- rep(counts, length(data))
  if(length(gene_id) == 1)
    gene_id <- rep(gene_id, length(data))
  if(!is.list(filter) & length(filter) == 2)
    filter <- lapply(seq_along(data), \(x) filter)

  NK_against_subsets <- list()
  for (i in 1:length(data)) {
    DEGs <- get_degs(data = data[[i]], ID = ID[[i]], type = type[[i]],
                     counts = counts[[i]], method = method,
                     markers = markers$SYMBOL |> unique(),
                     filter = filter[[i]],
                     gene_id = gene_id[[i]], ...)[[dir]]

    NK_against_subsets[[i]] <- markers$SYMBOL[markers[[gene_id[[i]]]] %in% DEGs]
  }
  if((is.character(comb)) && (comb == "RRA")) {
    if(length(NK_against_subsets) == 1) return(unlist(NK_against_subsets))
    ag <- RobustRankAggreg::aggregateRanks(NK_against_subsets)
    return(ag$Name[ag$Score < s_thres])
  }else {
    return(Reduce(comb, NK_against_subsets))
  }
})

#' @rdname filter_subset_sig
setMethod("filter_subset_sig", signature(
  data = 'DGEList',
  markers = 'vector',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         markers,
         ID,
         type,
         counts = TRUE,
         dir = "UP",
         gene_id = "SYMBOL",
         method = "RP",
         comb = union,
         filter = c(10,10),
         s_thres = 0.05,
         ...) {

  stopifnot(is.vector(markers), is.character(dir),
            is.character(method), is.numeric(s_thres),
            is.character(gene_id))

  markers <- AnnotationDbi::select(org.Hs.eg.db, markers, gene_id, "SYMBOL")
  # tryCatch(
  #   {
  #     markers$ENSEMBL[match(c("KLRA1P", "TRBC1", "TRDC"), markers$SYMBOL)] <- c(
  #       "ENSG00000256667",
  #       "ENSG00000211751",
  #       "ENSG00000211829"
  #     )
  #   },
  #   error = function(e) {
  #     "not containing KLRA1P,TRBC1,TRDC"
  #   }
  # )

  DEGs <- get_degs(data = data, ID = ID, type = type,
                   counts = counts, method = method,
                   markers = markers$SYMBOL |> unique(),
                   filter = filter,
                   gene_id = gene_id, ...)[[dir]]

  NK_against_subsets <- markers$SYMBOL[markers[[gene_id]] %in% DEGs]

  return(NK_against_subsets)
})

#' @rdname filter_subset_sig
setMethod("filter_subset_sig", signature(
  data = 'ANY',
  markers = 'vector',
  ID = 'ANY',
  type = 'ANY'
),
function(data,
         markers,
         ID,
         type,
         counts = TRUE,
         dir = "UP",
         gene_id = "SYMBOL",
         method = "RP",
         comb = union,
         filter = c(10,10),
         s_thres = 0.05,
         ...) {

  stopifnot(is.vector(markers), is.character(dir),
            is.character(method), is.numeric(s_thres),
            is.character(gene_id))

  markers <- AnnotationDbi::select(org.Hs.eg.db, markers, gene_id, "SYMBOL")
  # tryCatch(
  #   {
  #     markers$ENSEMBL[match(c("KLRA1P", "TRBC1", "TRDC"), markers$SYMBOL)] <- c(
  #       "ENSG00000256667",
  #       "ENSG00000211751",
  #       "ENSG00000211829"
  #     )
  #   },
  #   error = function(e) {
  #     "not containing KLRA1P,TRBC1,TRDC"
  #   }
  # )

  DEGs <- get_degs(data = data, ID = ID, type = type,
                   counts = counts, method = method,
                   markers = markers$SYMBOL |> unique(),
                   filter = filter,
                   gene_id = gene_id, ...)[[dir]]

  NK_against_subsets <- markers$SYMBOL[markers[[gene_id]] %in% DEGs]

  return(NK_against_subsets)
})

