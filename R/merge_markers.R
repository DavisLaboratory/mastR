#' Merge markers tbl into one.
#'
#' Merge markers collected from different DB into one tbl with 'TRUE' and '-' to
#' indicate which DB each gene is from.
#'
#' @param markers_list list of markers to be merged, can be list of vectors or
#'                     GeneSetCollection object
#' @param plot logical, if to make UpSetR plot for given list
#'
#' @return A data frame with genes symbols and columns indicate where the genes
#'         come from
#' @export
#'
#' @examples
#' Msig <- get_msigdb_sig(
#'   pattern = "natural_killer_cell_mediated",
#'   ignore.case = TRUE
#' )
#' Panglao <- get_panglao_sig(type = "NK cells")
#' merge_markers(markers = list(NK_markers = NK_markers$HGNC_Symbol,
#'                              MSig = GSEABase::geneIds(Msig),
#'                              PanglaoDB = GSEABase::geneIds(Panglao)))
merge_markers <- function(markers_list, plot = FALSE) {

  if(all(sapply(markers_list, \(x) is(x, "GeneSet")))) {
    markers_list <- GSEABase::GeneSetCollection(markers_list)
  }else if(any(sapply(markers_list, \(x) is(x, "GeneSet")))) {
    stop("Elements in the list should all be the same class!")
  }

  if(is(markers_list, "GeneSetCollection"))
    markers_list <- GSEABase::geneIds(markers_list)

  if(is.null(names(markers_list)))
    stop("Please provide named list!")

  markers_list <- lapply(names(markers_list), \(x) {
    d <- cbind(markers_list[[x]], TRUE)
    colnames(d) <- c("Gene", x)
    d
    })

  markers <- Reduce(function(x, y) merge(x, y, all = TRUE),
                    markers_list) |>
    tidyr::as_tibble()

  markers$y <- NULL
  markers[is.na(markers)] <- "-"

  if (plot == TRUE) UpSetR::upset(ifelse(markers[, -1] == TRUE, 1, 0) |>
                                    as.data.frame(),
                                  nsets = dim(markers)[2] - 1) |> print()

  return(markers)
}
