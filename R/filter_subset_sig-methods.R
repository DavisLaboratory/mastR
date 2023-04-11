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

    NK_against_subsets <- lapply(seq_along(data), \(i) {
      NK_against_subsets <- filter_subset_sig(
        data = data[[i]],
        group_col = group_col[[i]],
        target_group = target_group[[i]],
        markers = markers,
        normalize = normalize[[i]],
        dir = dir,
        gene_id = gene_id[[i]],
        feature_selection = feature_selection,
        comb = comb,
        filter = filter[[i]],
        s_thres = s_thres,
        batch = batch[[i]],
        slot = slot[[i]],
        ...
      )
    })

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

    DEGs <- get_degs(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      markers = markers,
      filter = filter,
      gene_id = gene_id,
      ...
    )
    DEGs <- GSEABase::geneIds(DEGs$DEGs[[dir]])

    if (is.null(markers)) {
      NK_against_subsets <- DEGs
    } else {
      stopifnot("Please provide a vector of gene symbols!" = is.vector(markers))
      markers <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        markers,
        gene_id, "SYMBOL"
      )
      idx <- na.omit(match(DEGs, markers[[gene_id]]))
      NK_against_subsets <- markers$SYMBOL[idx]
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

    DEGs <- get_degs(
      data = data,
      group_col = group_col,
      target_group = target_group,
      normalize = normalize,
      feature_selection = feature_selection,
      markers = markers,
      filter = filter,
      gene_id = gene_id, ...
    )
    DEGs <- GSEABase::geneIds(DEGs$DEGs[[dir]])

    if (is.null(markers)) {
      NK_against_subsets <- DEGs
    } else {
      stopifnot("Please provide a vector of gene symbols!" = is.vector(markers))
      markers <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        markers,
        gene_id, "SYMBOL"
      )
      idx <- na.omit(match(DEGs, markers[[gene_id]]))
      NK_against_subsets <- markers$SYMBOL[idx]
    }

    return(NK_against_subsets)
  }
)
