#' Collect genes from MSigDB.
#'
#' Collect genes of relevant MSigDB genesets.
#'
#' @param pattern pattern, to be matched the MsigDB gs_name of interest,
#'                e.g. 'NATURAL_KILLER_CELL_MEDIATED'
#' @param species chr, species of interest, can be 'hs' or 'mm'
#' @param cat chr, stating the category(s) to be retrieved.
#'            The category(s) must be one from [msigdb::listCollections()]
#' @param subcat chr, stating the sub-category(s) to be retrieved.
#'               The sub-category(s) must be one from [msigdb::listSubCollections()]
#' @param id a character, representing the ID type to use ("SYM" for gene
#'           symbols and "EZID" for Entrez IDs)
#' @param version a character, stating the version of MSigDB to be retrieved
#'                (should be >= 7.2). See [msigdb::getMsigdbVersions()].
#' @param plot logical, if to plot UpSetR diagram
#' @param ... params for [grep()], used to match pattern to gs_name
#'
#' @return A GeneSet object
#' @export
#'
#' @examples
#' get_msigdb_sig(
#'   pattern = "natural_killer_cell_mediated",
#'   species = "hs", subcat = "GO:BP",
#'   version = '7.4',
#'   ignore.case = TRUE
#' )
get_msigdb_sig <- function(pattern,
                           species = c("hs", "mm"),
                           cat = NULL,
                           subcat = NULL,
                           id = c("SYM", "EZID"),
                           version = msigdb::getMsigdbVersions(),
                           plot = FALSE, ...) {
  ## load msigdb from ExperimentHub
  msigdb <- msigdb::getMsigdb(org = species, id = id, version = version)

  ## subset msigdb if given cat or subcat
  if(!(is.null(cat) & is.null(subcat))) {
     msigdb <- msigdb::subsetCollection(msigdb,
                                        collection = cat,
                                        subcollection = subcat)
  }

  msig_gs <- grep(pattern, names(msigdb), ...)
  if(length(msig_gs) == 0){
    message("No matched gene set is found in MSigDB!")
    return()
  }
  msigdb <- msigdb[msig_gs]

  if(plot) {
    if(length(msigdb) > 1) {
      UpSetR::upset(UpSetR::fromList(msigdb |> GSEABase::geneIds()),
                    nsets = length(msigdb)) |> print()
    }else message("Only one gene-set is matched!")
  }

  msigdb <- Reduce('|', msigdb)
  msigdb@setName <- paste("MSigDB", pattern, sep = "_")
  return(msigdb)
}
