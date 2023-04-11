#' @rdname get_gsc_sig
setMethod(
  "get_gsc_sig", signature(
    gsc = "GeneSetCollection",
    pattern = "character"
  ),
  function(gsc = "msigdb",
           pattern,
           cat = NULL,
           subcat = NULL,
           ...) {
    ## subset gsc if given cat or subcat
    if (!(is.null(cat) && is.null(subcat))) {
      gsc <- msigdb::subsetCollection(gsc,
        collection = cat,
        subcollection = subcat
      )
    }

    sub_gsc <- grep(pattern, names(gsc), ...)
    if (length(sub_gsc) == 0) {
      stop("No matched gene set is found in gsc!")
    }
    gsc <- gsc[sub_gsc]

    return(gsc)
  }
)

#' @rdname get_gsc_sig
setMethod(
  "get_gsc_sig", signature(
    gsc = "character",
    pattern = "character"
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
  }
)
