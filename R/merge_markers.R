#' @importFrom jsonlite toJSON fromJSON
NULL

#' Merge markers list into one.
#'
#' Merge markers collected from different DB into one 'GeneSet' object, saved a
#' data.frame in json format under `longDescription` with 'TRUE' and '-' to
#' indicate which DB each gene is from, this can be shown via [jsonlite::fromJSON()]
#'
#' @param markers_list list of markers to be merged, can be a list of
#'                     'GeneSet'/vector or a 'GeneSetCollection' object
#' @param plot logical, if to make UpSetR plot for given list
#'
#' @return A GeneSet class of union genes in the given list
#' @export
#'
#' @examples
#' Panglao <- get_panglao_sig(type = "NK cells")
#' merge_markers(markers = list(NK_markers = NK_markers$HGNC_Symbol,
#'                              PanglaoDB = GSEABase::geneIds(Panglao)))
merge_markers <- function(markers_list, plot = FALSE) {

  stopifnot("markers_list must be a list!" = is.list(markers_list),
            "plot is not a logical value!" = is.logical(plot),
            "provide more than 1 value for plot!" = length(plot) == 1)

  if(is.null(names(markers_list)))
    stop("Please provide named list!")

  ## convert list into 'GeneSetCollection'
  if(all(sapply(markers_list, \(x) is(x, "GeneSet")))) {
    markers_list <- GSEABase::GeneSetCollection(markers_list)
  }else if(all(sapply(markers_list, is.vector))) {
    markers_list <- lapply(names(markers_list),
                           \(n) GSEABase::GeneSet(markers_list[[n]],
                                                  setName = n)) |>
      GSEABase::GeneSetCollection()
  }else {
    stop("Elements in the list should all be the same class!")
  }

  ## merge and get union of all genesets
  markers <- BiocGenerics::Reduce('|', markers_list)
  GSEABase::setName(markers) <- "merged_markers_pool"

  markers_list <- GSEABase::geneIds(markers_list)

  ## plot upset diagram for overlapping genes among gene-sets in the list
  if(plot == TRUE) {
    UpSetR::upset(UpSetR::fromList(markers_list),
                  nsets = length(markers_list)) |> print()
  }

  ## create tibble to show merged genes' origin
  markers_list <- lapply(names(markers_list), \(x) {
    d <- cbind(markers_list[[x]], TRUE)
    colnames(d) <- c("Gene", x)
    d
    })
  markers_list <- Reduce(function(x, y) merge(x, y, all = TRUE),
                         markers_list) |>
    tidyr::as_tibble()

  markers_list$y <- NULL
  markers_list[is.na(markers_list)] <- "-"

  ## save the tbl into GeneSet
  GSEABase::longDescription(markers) <- jsonlite::toJSON(markers_list) |>
    as.character()  ## convert tbl into json character format

  return(markers)
}
