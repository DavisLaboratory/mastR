#' Filter specific cell type signature genes against other subsets.
#'
#' Specify the signature of the subset matched 'target_group' against other subsets,
#' either "union", "intersect" or "RRA" can be specified when input is a list
#' of datasets to integrate the signatures into one.
#'
#' @inheritParams process_data
#' @inheritParams select_sig
#' @param data An expression data or a list of expression data objects
#' @param group_col vector or character, specify the group factor or column name of
#'           coldata for DE comparisons
#' @param dir character, could be 'UP' or 'DOWN' to use only up- or
#'            down-expressed genes
#' @param comb 'RRA' or Fun for combining sigs from multiple datasets, keep all
#'             passing genes or only intersected genes, could be `union` or
#'             `intersect` or `setdiff` or customized Fun, or could be 'RRA' to
#'             use Robust Rank Aggregation method for integrating multi-lists
#'             of sigs, default 'union'
#' @param filter (list of) vector of 2 numbers, filter condition to remove low
#'               expression genes, the 1st for min.counts (if normalize = TRUE)
#'               or CPM/TPM (if normalize = FALSE), the 2nd for samples size 'large.n'
#' @param s_thres num, threshold of score if comb = 'RRA'
#' @param ... other params for [get_degs()]
#'
#' @return a vector of gene symbols
#'
#' @examples
#' data("im_data_6", "NK_markers")
#' sigs <- filter_subset_sig(im_data_6, "celltype:ch1", "NK",
#'   markers = NK_markers$HGNC_Symbol,
#'   gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric(
  "filter_subset_sig",
  function(data,
           group_col,
           target_group,
           markers = NULL,
           normalize = TRUE,
           dir = "UP",
           gene_id = "SYMBOL",
           feature_selection = c("auto", "rankproduct", "none"),
           comb = union,
           filter = c(10, 10),
           s_thres = 0.05,
           ...) {
    standardGeneric("filter_subset_sig")
  }
)

#' @rdname filter_subset_sig
setMethod(
  "filter_subset_sig", signature(
    data = "list",
    group_col = "ANY",
    target_group = "ANY"
  ),
  function(data,
           group_col,
           target_group,
           markers = NULL,
           normalize = TRUE,
           dir = "UP",
           gene_id = "SYMBOL",
           feature_selection = c("auto", "rankproduct", "none"),
           comb = union,
           filter = c(10, 10),
           s_thres = 0.05,
           slot = "counts",
           batch = NULL,
           ...) {
    stopifnot(
      is.character(dir), is.character(feature_selection),
      is.numeric(s_thres), is.character(gene_id)
    )

    if (!is.null(markers)) {
      stopifnot("Please provide a vector of gene symbols!" = is.vector(markers))
      markers <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        markers,
        gene_id, "SYMBOL"
      )
    }
    # tryCatch(
    #   {
    #     markers$ENSEMBL[match(c("KLRA1P", "TRBC1", "TRDC"),
    #                           markers$SYMBOL)] <- c(
    #       "ENSG00000256667",
    #       "ENSG00000211751",
    #       "ENSG00000211829"
    #     )
    #   },
    #   error = function(e) {
    #     "not containing KLRA1P,TRBC1,TRDC"
    #   }
    # )

    if (length(group_col) == 1) {
      group_col <- rep(group_col, length(data))
    }
    if (length(target_group) == 1) {
      target_group <- rep(target_group, length(data))
    }
    if (length(normalize) == 1) {
      normalize <- rep(normalize, length(data))
    }
    if (length(slot) == 1) {
      slot <- rep(slot, length(data))
    }
    if (length(batch) == 1) {
      batch <- rep(batch, length(data))
    }
    if (length(gene_id) == 1) {
      gene_id <- rep(gene_id, length(data))
    }
    if (!is.list(filter) & length(filter) == 2) {
      filter <- lapply(seq_along(data), \(x) filter)
    }

    NK_against_subsets <- list()
    for (i in seq_along(data)) {
      DEGs <- get_degs(
        data = data[[i]], group_col = group_col[[i]],
        target_group = target_group[[i]],
        normalize = normalize[[i]],
        slot = slot[[i]],
        batch = batch[[i]],
        feature_selection = feature_selection,
        markers = unique(markers$SYMBOL),
        filter = filter[[i]],
        gene_id = gene_id[[i]], ...
      )
      DEGs <- GSEABase::geneIds(DEGs$DEGs[[dir]])

      if (is.null(markers)) {
        NK_against_subsets[[i]] <- DEGs
      } else {
        NK_against_subsets[[i]] <- markers$SYMBOL[markers[[gene_id[[i]]]] %in% DEGs]
      }
    }
    if ((is.character(comb)) && (comb == "RRA")) {
      if (length(NK_against_subsets) == 1) {
        return(unlist(NK_against_subsets))
      }
      ag <- RobustRankAggreg::aggregateRanks(NK_against_subsets)
      return(ag$Name[ag$Score < s_thres])
    } else {
      return(Reduce(comb, NK_against_subsets))
    }
  }
)

#' @rdname filter_subset_sig
setMethod(
  "filter_subset_sig", signature(
    data = "DGEList",
    group_col = "ANY",
    target_group = "ANY"
  ),
  function(data,
           group_col,
           target_group,
           markers = NULL,
           normalize = TRUE,
           dir = "UP",
           gene_id = "SYMBOL",
           feature_selection = c("auto", "rankproduct", "none"),
           comb = union,
           filter = c(10, 10),
           s_thres = 0.05,
           ...) {
    stopifnot(
      is.character(dir), is.character(feature_selection),
      is.numeric(s_thres), is.character(gene_id)
    )

    if (!is.null(markers)) {
      stopifnot("Please provide a vector of gene symbols!" = is.vector(markers))
      markers <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        markers,
        gene_id, "SYMBOL"
      )
    }
    # tryCatch(
    #   {
    #     markers$ENSEMBL[match(c("KLRA1P", "TRBC1", "TRDC"),
    #                           markers$SYMBOL)] <- c(
    #       "ENSG00000256667",
    #       "ENSG00000211751",
    #       "ENSG00000211829"
    #     )
    #   },
    #   error = function(e) {
    #     "not containing KLRA1P,TRBC1,TRDC"
    #   }
    # )

    DEGs <- get_degs(
      data = data, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      markers = unique(markers$SYMBOL),
      filter = filter,
      gene_id = gene_id, ...
    )
    DEGs <- GSEABase::geneIds(DEGs$DEGs[[dir]])

    if (is.null(markers)) {
      NK_against_subsets <- DEGs
    } else {
      NK_against_subsets <- markers$SYMBOL[markers[[gene_id]] %in% DEGs]
    }

    return(NK_against_subsets)
  }
)

#' @rdname filter_subset_sig
setMethod(
  "filter_subset_sig", signature(
    data = "ANY",
    group_col = "ANY",
    target_group = "ANY"
  ),
  function(data,
           group_col,
           target_group,
           markers = NULL,
           normalize = TRUE,
           dir = "UP",
           gene_id = "SYMBOL",
           feature_selection = c("auto", "rankproduct", "none"),
           comb = union,
           filter = c(10, 10),
           s_thres = 0.05,
           ...) {
    stopifnot(
      is.character(dir), is.character(feature_selection),
      is.numeric(s_thres), is.character(gene_id)
    )

    if (!is.null(markers)) {
      stopifnot("Please provide a vector of gene symbols!" = is.vector(markers))
      markers <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        markers,
        gene_id, "SYMBOL"
      )
    }
    # tryCatch(
    #   {
    #     markers$ENSEMBL[match(c("KLRA1P", "TRBC1", "TRDC"),
    #                           markers$SYMBOL)] <- c(
    #       "ENSG00000256667",
    #       "ENSG00000211751",
    #       "ENSG00000211829"
    #     )
    #   },
    #   error = function(e) {
    #     "not containing KLRA1P,TRBC1,TRDC"
    #   }
    # )

    DEGs <- get_degs(
      data = data, group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      markers = unique(markers$SYMBOL),
      filter = filter,
      gene_id = gene_id, ...
    )
    DEGs <- GSEABase::geneIds(DEGs$DEGs[[dir]])

    if (is.null(markers)) {
      NK_against_subsets <- DEGs
    } else {
      NK_against_subsets <- markers$SYMBOL[markers[[gene_id]] %in% DEGs]
    }

    return(NK_against_subsets)
  }
)
