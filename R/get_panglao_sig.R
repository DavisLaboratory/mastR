#' Extract immune subset markers from PanglaoDB website.
#'
#' Extract specific immune subset markers for 'Hs' or 'Mm', the markers are
#' retrived from up-to-date PanglaoDB website.
#'
#' @param species chr, default 'Hs', could be 'Hs', 'Mm' or 'Mm Hs',
#'                specify the species of interest
#' @param type character vector, cell type name(s) of interest,
#'             available subsets could be listed by [list_panglao_types()]
#'
#' @return a vector of markers
#' @export
#'
#' @examples
#' get_panglao_sig(species = "Hs", type = "NK cells")
get_panglao_sig <- function(species = "Hs", type) {

  web <- rvest::read_html("https://panglaodb.se/markers.html?cell_type='all_cells'")
  table <- rvest::html_table(web)[[1]]

  stopifnot("Please provide valid species, either 'Hs' or 'Mm' or 'Mm Hs'!" = grepl("Hs|Mm", species),
            "Please provide valid cell types, which can be listed by list_Panglao_types()!" = type %in% table$`Cell type`)

  markers <- subset(table, grepl(species, Species) & `Cell type` %in% type)
  markers <- markers$`Official gene symbol` |> unique()
  return(markers)
}
