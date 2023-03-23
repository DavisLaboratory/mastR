#' Convert gene-set list into GeneSetCollection
#'
#' @param ... vector of genes or list of genes
#'
#' @return GeneSetCollection
#' @export
#'
#' @examples
#' data("msigdb_gobp_nk")
#' gls2gsc(GSEABase::geneIds(msigdb_gobp_nk[1:3]))
setGeneric(
  "gls2gsc",
  function(...) {
    standardGeneric("gls2gsc")
  }
)

#' @rdname gls2gsc
setMethod(
  "gls2gsc", signature(
    "list"
  ),
  function(...) {
    gsc <- c(...)
    ## name gene_list
    if (is.null(names(gsc))) {
      names(gsc) <- seq_along(gsc)
    }
    idx <- which(names(gsc) == "")
    names(gsc)[idx] <- idx
    ## convert list into gsc
    gsc <- lapply(seq_along(gsc), function(i) {
      GSEABase::GeneSet(gsc[[i]], setName = names(gsc)[i])
    })
    gsc <- GSEABase::GeneSetCollection(gsc)

    return(gsc)
  }
)

#' @rdname gls2gsc
setMethod(
  "gls2gsc", signature(
    "vector"
  ),
  function(...) {
    gsc <- make_named_list(...)
    gsc <- gls2gsc(gsc)

    return(gsc)
  }
)

# helper: make named list
make_named_list <- function(...) {
  structure(list(...), names = as.list(substitute(list(...)))[-1L])
}
