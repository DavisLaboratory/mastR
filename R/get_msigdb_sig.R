#' Collect genes from MSigDB.
#'
#' Collect genes of relevant MSigDB gene-sets, of which the gene-set names are
#' matched to the given pattern by using [grep()] function. By setting cat and
#' subcat, matching can be constrained in the union of given categories and
#' subcategories.
#'
#' @param pattern pattern pass to [grep()], to match the MsigDB gene-set name of
#'                interest, e.g. 'NATURAL_KILLER_CELL_MEDIATED'
#' @param species character, species of interest, can be 'hs' or 'mm'
#' @param cat character, stating the category(s) to be retrieved.
#'            The category(s) must be one from [msigdb::listCollections()],
#'            see details in [msigdb::subsetCollection()]
#' @param subcat character, stating the sub-category(s) to be retrieved.
#'               The sub-category(s) must be one from [msigdb::listSubCollections()],
#'               see details in [msigdb::subsetCollection()]
#' @param id a character, representing the ID type to use ("SYM" for gene
#'           symbols and "EZID" for Entrez IDs)
#' @param version a character, stating the version of MSigDB to be retrieved
#'                (should be >= 7.2). See [msigdb::getMsigdbVersions()].
#' @param plot logical, if to plot UpSetR diagram
#' @param ... params for [grep()], used to match pattern to gene-set names
#'
#' @return A GeneSet object containing all matched gene-sets in MSigDB
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

  stopifnot(is.character(pattern), is.logical(plot))

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
    stop("No matched gene set is found in MSigDB!")
  }
  msigdb <- msigdb[msig_gs]

  if(plot == TRUE) {
    if(length(msigdb) > 1) {
      UpSetR::upset(UpSetR::fromList(msigdb |> GSEABase::geneIds()),
                    nsets = length(msigdb)) |> print()
    }else warning("Only one gene-set is matched!")
  }

  msigdb <- Reduce('|', msigdb)
  GSEABase::setName(msigdb) <- paste("MSigDB", pattern, sep = "_")[1]
  return(msigdb)
}
