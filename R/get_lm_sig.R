#' Extract specific subset markers from LM7 or/and LM22
#'
#' @param lm7.pattern character string containing a regular expression
#'                    (or character string for fixed = TRUE), to be matched in
#'                    the given subsets in LM7
#' @param lm22.pattern character string containing a regular expression
#'                    (or character string for fixed = TRUE), to be matched in
#'                    the given subsets in LM22
#' @param ... params for function [grep()]
#'
#' @return A vector of markers
#' @export
#'
#' @examples
#' get_lm_sig(lm7.pattern = "NK", lm22.pattern = "NK cells")
get_lm_sig <- function(lm7.pattern = FALSE, lm22.pattern = FALSE, ...){


  if(lm7.pattern == FALSE & lm22.pattern == FALSE)
    stop("Must provide at least one of lm7.pattern and lm22.pattern param!")

  if(lm7.pattern != FALSE){
    markers_7 <- subset(LM7, grepl(lm7.pattern, Subset, ...))$Gene
    markers <- markers_7
  }
  if(lm22.pattern != FALSE){
    idx <- grep(lm22.pattern, colnames(LM22)[-1], ...)
    idx <- Reduce(function(x,y){x == 1 | y == 1}, LM22[,-1][,idx]) |>
      as.logical()
    markers_22 <- LM22$Gene[idx]
    markers <- markers_22
  }
  if(lm7.pattern != FALSE & lm22.pattern != FALSE){
    markers <- union(markers_7, markers_22)
  }
  return(markers)
}
