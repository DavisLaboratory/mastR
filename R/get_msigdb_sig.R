#' Collect genes from MSigDB.
#'
#' Collect genes of relevant MSigDB genesets.
#'
#' @param species chr, species of interest, available species can be listed by
#'                [msigdbr::msigdbr_species()]
#' @param cat chr, category of interest, available collections can be listed by
#'            [msigdbr::msigdbr_collections()]
#' @param subcat chr, subcategory of interest, available collections can be
#'               listed by [msigdbr::msigdbr_collections()]
#' @param pattern pattern, to be matched the MsigDB gs_name of interest,
#'                e.g. 'natural_killer_cell_mediated'
#' @param plot logical, if to plot UpSetR diagram
#' @param ... params for [grep()], used to match pattern to gs_name
#'
#' @return A GeneSet object
#' @export
#'
#' @examples
#' get_msigdb_sig(
#'   species = "Homo sapiens", cat = "C5", subcat = "GO:BP",
#'   pattern = "natural_killer_cell_mediated", ignore.case = TRUE
#' )
get_msigdb_sig <- function(species = "Homo sapiens", cat = "C5",
                           subcat = "GO:BP", pattern,
                           plot = FALSE, ...) {
  msigdb <- msigdbr::msigdbr(species = species, category = cat,
                             subcategory = subcat)
  msig_terms <- grep(pattern, msigdb$gs_name, value = T, ...) |> unique()
  if(length(msig_terms) == 0){
    message("No relevant term was found!")
    return()
  }
  msig_set <- msigdb$gene_symbol[msigdb$gs_name %in% msig_terms] |> unique()

  tmp <- msigdb[msigdb$gs_name %in% msig_terms, c("gs_name", "gene_symbol")]
  if (plot) {
    ls <- split(tmp$gene_symbol, tmp$gs_name)
    if(length(ls) > 1) {
      UpSetR::upset(UpSetR::fromList(ls),
                    nsets = length(msig_terms)) |> print()
    }else message("Only one gene-set is matched!")
  }

  msig_set <- GSEABase::GeneSet(msig_set,
                                setName = paste("MSigDB", pattern, sep = "_"),
                                geneIdType = GSEABase::SymbolIdentifier())
  return(msig_set)
}
