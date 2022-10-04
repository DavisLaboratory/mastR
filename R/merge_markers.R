#' Merge markers tbl into one.
#'
#' Merge markers collected from different DB into one tbl with 'TRUE' and '-' to
#' indicate which DB each gene is from.
#'
#' @param markers_list list of tibble, each tibble containing at least one column
#'                     'HGNC_Symbol'
#' @param plot logical, if to make UpSetR plot
#'
#' @return A data frame with genes symbols and columns indicate where the genes
#'         come from
#' @export
#'
#' @examples
#' Msig <- get_MSig(
#'   species = "Homo sapiens", cat = "C5", subcat = "GO:BP",
#'   pattern = "natural_killer_cell_mediated"
#' )
#' Panglao <- get_Panglao(species = "Hs", type = "NK cells")
#' merge_markers(markers = list(NK_markers, Msig, Panglao))
merge_markers <- function(markers_list, plot = FALSE) {
  markers <- Reduce(function(x, y) merge(x, y, all = T),
                    markers_list) |>
    tidyr::as_tibble()

  markers$y <- NULL
  markers[is.na(markers)] <- "-"

  if (plot == TRUE) UpSetR::upset(ifelse(markers[, -1] == TRUE, 1, 0) |>
                                    as.data.frame(),
                                  nsets = dim(markers)[2] - 1) |> print()
  return(markers)
}
