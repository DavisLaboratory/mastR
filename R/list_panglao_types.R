#' Show the summary info of available cell types in PanglaoDB.
#'
#' Show the name and number of each cell type in PanglaoDB. Help users
#' know which subset(s) marker list(s) could be retrieved by PanglaoDB.
#'
#' @param organ character, specify the tissue or organ label to list cell types
#'
#' @return a vector of available cell types of the organ in PanglaoDB
#' @export
#'
#' @examples
#' list_panglao_types(organ = "Immune system")
list_panglao_types <- function(organ) {
  organ <- match.arg(organ, choices = list_panglao_organs())

  av_organs <- list_panglao_organs()
  stopifnot("Please provide a valid organ!" = organ %in% av_organs)

  web <- rvest::read_html("https://panglaodb.se/markers.html?cell_type=%27choose%27")

  av_cells <- web |>
    rvest::html_elements(paste0("optgroup[label='", organ, "'] option")) |>
    rvest::html_text()
  return(av_cells)
}
