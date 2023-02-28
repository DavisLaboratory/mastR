#' Collect genes from MSigDB or provided GeneSetCollection.
#'
#' Collect genesets from MSigDB or given GeneSetCollection, of which the gene-set
#' names are matched to the given regex pattern by using [grep()] function.
#' By setting cat and subcat, matching can be constrained in the union of given
#' categories and subcategories if gsc = 'msigdb'.
#'
#' @param gsc 'msigdb' or GeneSetCollection to be searched
#' @param pattern pattern pass to [grep()], to match the MsigDB gene-set name of
#'                interest, e.g. 'NATURAL_KILLER_CELL_MEDIATED'
#' @param cat character, stating the category(s) to be retrieved.
#'            The category(s) must be one from [msigdb::listCollections()],
#'            see details in [msigdb::subsetCollection()]
#' @param subcat character, stating the sub-category(s) to be retrieved.
#'               The sub-category(s) must be one from
#'               [msigdb::listSubCollections()], see details in
#'               [msigdb::subsetCollection()]
#' @param species character, species of interest, can be 'hs' or 'mm'
#' @param id a character, representing the ID type to use ("SYM" for gene
#'           symbols and "EZID" for Entrez IDs)
#' @param version a character, stating the version of MSigDB to be retrieved
#'                (should be >= 7.2). See [msigdb::getMsigdbVersions()].
#' @param ... params for [grep()], used to match pattern to gene-set names
#'
#' @return A GeneSet object containing all matched gene-sets in MSigDB
#' @export
#'
#' @examples
#' get_gsc_sig(
#'   gsc = msigdb_gobp_nk,
#'   pattern = "natural_killer_cell_mediated",
#'   subcat = "GO:BP",
#'   ignore.case = TRUE
#' )
setGeneric("get_gsc_sig",
           function(gsc = "msigdb",
                    pattern,
                    cat = NULL,
                    subcat = NULL,
                    species = c("hs", "mm"),
                    id = c("SYM", "EZID"),
                    version = msigdb::getMsigdbVersions(),
                    ...)
           standardGeneric("get_gsc_sig"))

#' @rdname get_gsc_sig
setMethod("get_gsc_sig", signature(
  gsc = 'GeneSetCollection',
  pattern = 'character'
),
function(gsc = "msigdb",
         pattern,
         cat = NULL,
         subcat = NULL,
         ...) {

  ## subset gsc if given cat or subcat
  if(!(is.null(cat) && is.null(subcat))) {
    gsc <- msigdb::subsetCollection(gsc,
                                    collection = cat,
                                    subcollection = subcat)
  }

  sub_gsc <- grep(pattern, names(gsc), ...)
  if(length(sub_gsc) == 0) {
    stop("No matched gene set is found in gsc!")
  }
  gsc <- gsc[sub_gsc]

  return(gsc)
})

#' @rdname get_gsc_sig
setMethod("get_gsc_sig", signature(
  gsc = 'character',
  pattern = 'character'
),
function(gsc = "msigdb",
         pattern,
         cat = NULL,
         subcat = NULL,
         species = c("hs", "mm"),
         id = c("SYM", "EZID"),
         version = msigdb::getMsigdbVersions(),
         ...) {

  stopifnot("gsc must be 'msigdb' or GeneSetCollection!" = gsc == "msigdb")

  ## load msigdb from ExperimentHub
  gsc <- msigdb::getMsigdb(org = species, id = id, version = version)

  gsc <- get_gsc_sig(gsc = gsc, pattern = pattern, cat = cat, subcat = subcat)
  return(gsc)
})
