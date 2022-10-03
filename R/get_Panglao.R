#' Extract immune subset markers from PanglaoDB website.
#'
#' Extract specific immune subset markers for 'Hs' or 'Mm', the markers are
#' retrived from up-to-date PanglaoDB website.
#'
#' @param species chr, default 'Hs', could be 'Hs', 'Mm' or 'Mm Hs',
#'                specify the species of interest
#' @param type character vector, cell type name(s) of interest,
#'             available subsets could be listed by [list_Panglao_types()]
#'
#' @return a tibble, first column is gene symbol, second column is markers origin
#' @export
#'
#' @examples
#' get_Panglao(species = "Hs", type = "NK cells")
get_Panglao <- function(species = "Hs", type) {

  web <- rvest::read_html("https://panglaodb.se/markers.html?cell_type='all_cells'")
  table <- rvest::html_table(web)[[1]]

  stopifnot("Please provide valid species, either 'Hs' or 'Mm' or 'Mm Hs'!" = grepl("Hs|Mm", species),
            "Please provide valid cell types, which can be listed by list_Panglao_types()!" = type %in% table$`Cell type`)

  markers <- subset(table, grepl(species, Species) & `Cell type` %in% type)
  markers <- tibble::tibble(HGNC_Symbol = markers$`Official gene symbol`,
                            PanglaoDB = T)
  markers$PanglaoDB <- as.character(markers$PanglaoDB)
  markers <- markers[!duplicated(markers),]  ## remove duplicated or overlapped genes
  return(markers)
}
