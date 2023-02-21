#' Extract specific subset markers from LM7 or/and LM22
#'
#' Extract markers for subsets matched to the given pattern from LM7/LM22, and
#' save the matched genes in 'GeneSet' class object, if both pattern are
#' provided, the output would be a 'GeneSetCollection' class object with
#' setName: LM7, LM22.
#'
#' @param lm7.pattern character string containing a regular expression,
#'                    to be matched in the given subsets in LM7
#' @param lm22.pattern character string containing a regular expression,
#'                     to be matched in the given subsets in LM22
#' @param ... params for function [grep()]
#'
#' @return A GeneSet or GeneSetCollection for matched subsets in LM7 and/or LM22
#' @export
#'
#' @examples
#' data("LM7", "LM22")
#' get_lm_sig(lm7.pattern = "NK", lm22.pattern = "NK cells")
get_lm_sig <- function(lm7.pattern, lm22.pattern) {


  if(missing(lm7.pattern) && missing(lm22.pattern))
    stop("Must provide at least one of lm7.pattern and lm22.pattern param!")

  ## retrieve markers from LM7 signature matrix if lm7.pattern is given
  if(!missing(lm7.pattern)) {
    stopifnot(is.character(lm7.pattern))
    LM7 <- mastR::LM7
    keep <- grepl(lm7.pattern, LM7$Subset, ...)
    gs_7 <- LM7$Gene[keep]
    gs_7 <- GSEABase::GeneSet(gs_7, setName = "LM7",
                              geneIdType = GSEABase::SymbolIdentifier())
    gs <- gs_7
  }

  ## retrieve markers from LM22 signature matrix if lm22.pattern is given
  if(!missing(lm22.pattern)) {
    stopifnot(is.character(lm22.pattern))

    idx <- grep(lm22.pattern, colnames(LM22)[-1], ...)
    idx <- Reduce(function(x, y) {x == 1 | y == 1}, LM22[,-1][,idx]) |>
      as.logical()
    gs_22 <- LM22$Gene[idx]
    gs_22 <- GSEABase::GeneSet(gs_22, setName = "LM22",
                              geneIdType = GSEABase::SymbolIdentifier())
    gs <- gs_22
  }
  if(!missing(lm7.pattern) && !missing(lm22.pattern)) {
    gs <- GSEABase::GeneSetCollection(gs_7, gs_22)
  }
  return(gs)
}

utils::globalVariables(c("LM7", "LM22", "Subset"))
