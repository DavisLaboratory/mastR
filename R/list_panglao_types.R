#' Show the summary info of available cell types in PanglaoDB.
#'
#' Show the name and number of each cell type in PanglaoDB. Help users
#' know which subset(s) marker list(s) could be retrived by PanglaoDB.
#'
#' @param type chr, specify if to list all available organ types or cell types
#'             for specific organ
#' @param organ chr, specify the tissue or organ label to list cell types,
#'              optional when type = 'cell'
#'
#' @return a vector of available organ types or cell types in PanglaoDB
#' @export
#'
#' @examples
#' list_panglao_types(organ = "Immune system")
list_panglao_types <- function(type = "organ", organ) {

  stopifnot("Please provide a valid type, either 'cell' or 'organ'!" = type %in% c('cell', "organ"))

  web <- rvest::read_html("https://panglaodb.se/markers.html?cell_type=%27choose%27")

  av_organs <- rvest::html_elements(web, "optgroup") |>
    rvest::html_attr("label")
  if(type == "organ") {
    return(av_organs)
  }else{
    stopifnot("Please provide a valid organ!" = organ %in% av_organs)

    av_cells <- web |>
    rvest::html_elements(paste0("optgroup[label='",organ,"'] option")) |>
    rvest::html_text()
    return(av_cells)
  }
}
