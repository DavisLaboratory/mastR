#' Show the summary info of available organs in PanglaoDB.
#'
#' Show the name of organs available in PanglaoDB. Help users know which organs
#' could be retrived by PanglaoDB.
#'
#' @return a vector of available organ types or cell types in PanglaoDB
#' @export
#'
#' @examples
#' list_panglao_organs()
list_panglao_organs <- function() {
  web <- rvest::read_html("https://panglaodb.se/markers.html?cell_type=%27choose%27")

  av_organs <- rvest::html_attr(rvest::html_elements(web, "optgroup"), "label")

  return(av_organs)
}
