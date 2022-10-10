#' Merge markers tbl into one.
#'
#' Merge markers collected from different DB into one tbl with 'TRUE' and '-' to
#' indicate which DB each gene is from.
#'
#' @param markers_list list of markers to be merged
#' @param plot logical, if to make UpSetR plot for given list
#'
#' @return A data frame with genes symbols and columns indicate where the genes
#'         come from
#' @export
#'
#' @examples
#' Msig <- get_msigdb_sig(
#'   species = "Homo sapiens", cat = "C5", subcat = "GO:BP",
#'   pattern = "natural_killer_cell_mediated"
#' )
#' Panglao <- get_panglao_sig(species = "Hs", type = "NK cells")
#' merge_markers(markers = list(NK_markers = NK_markers$HGNC_Symbol,
#'                              MSig = Msig,
#'                              PanglaoDB = Panglao))
merge_markers <- function(markers_list, plot = FALSE) {

  if(is.null(names(markers_list)))
    stop("Please provide named list!")

  markers_list <- lapply(names(markers_list), \(x) {
    d <- cbind(markers_list[[x]], TRUE)
    colnames(d) <- c("Gene", x)
    d
    })

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
