#' Merge markers list into one.
#'
#' Merge markers collected from different DB into one 'GeneSet' object, saved a
#' data.frame in json format under `longDescription` with 'TRUE' and '-' to
#' indicate which DB each gene is from, this can be shown via
#' [jsonlite::fromJSON()].
#'
#' @param ... GeneSet or GeneSetCollection object to be merged
#'
#' @return A GeneSet class of union genes in the given list
#' @export
#'
#' @examples
#' Markers <- merge_markers(mastR::msigdb_gobp_nk[1:3])
#' jsonlite::fromJSON(GSEABase::longDescription(Markers))
merge_markers <- function(...) {

  ## input must be GeneSet or GeneSetCollection
  if(!all(sapply(list(...), class) %in% c("GeneSet", "GeneSetCollection")))
    stop("Only accept GeneSet or GeneSetCollection as input")

  ## convert into 'GeneSetCollection'
  gsc <- GSEABase::GeneSetCollection(c(...))

  ## merge and get union of all genesets
  markers <- Reduce('|', gsc)
  GSEABase::setName(markers) <- "merged_markers_pool"

  markers_list <- GSEABase::geneIds(gsc)

  ## create tibble to show merged genes' origin
  markers_list <- lapply(names(markers_list), \(x) {
    d <- cbind(markers_list[[x]], TRUE)
    colnames(d) <- c("Gene", x)
    d
    })
  markers_list <- as.data.frame(Reduce(function(x, y) merge(x, y, all = TRUE),
                                       markers_list))

  markers_list$y <- NULL
  markers_list[is.na(markers_list)] <- "-"

  ## save the tbl into GeneSet
  ## convert tbl into json character format
  GSEABase::longDescription(markers) <- as.character(jsonlite::toJSON(markers_list))

  return(markers)
}
