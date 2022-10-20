#' Extract immune subset markers from PanglaoDB website.
#'
#' Extract specific immune subset markers for 'Hs' or 'Mm', the markers are
#' retrived from up-to-date PanglaoDB website.
#'
#' @param type character vector, cell type name(s) of interest,
#'             available subsets could be listed by [list_panglao_types()]
#' @param species character, default 'Hs', could be 'Hs', 'Mm' or 'Mm Hs',
#'                specify the species of interest
#'
#' @return a 'GeneSet' class object containing genes of given type(s)
#' @export
#'
#' @examples
#' get_panglao_sig(type = "NK cells")
#' get_panglao_sig(type = c("NK cells", "T cells"))
get_panglao_sig <- function(type, species = c('Hs', 'Mm', 'Mm Hs')) {

  stopifnot(is.character(type))

  species <- match.arg(species)

  web <- rvest::read_html("https://panglaodb.se/markers.html?cell_type='all_cells'")
  table <- rvest::html_table(web)[[1]]

  stopifnot("Please provide valid species, either 'Hs' or 'Mm' or 'Mm Hs'!" = grepl("Hs|Mm", species),
            "Please provide valid cell types, which can be listed by list_Panglao_types()!" = type %in% table$`Cell type`)

  gs <- subset(table, grepl(species, Species) & `Cell type` %in% type)
  gs <- gs$`Official gene symbol` |> unique()
  gs <- GSEABase::GeneSet(gs, setName = "PanglaoDB",
                          geneIdType = GSEABase::SymbolIdentifier())
  return(gs)
}
