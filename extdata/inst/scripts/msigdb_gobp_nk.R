## code to prepare `msigdb_gobp_nk` dataset goes here

### load msigdb
msigdb <- msigdb::getMsigdb(org = "hs", id = "SYM", version = "7.4")
msigdb_subcollection <- msigdb::subsetCollection(msigdb,
  collection = "c5",
  subcollection = "GO:BP"
)
### grep NK relevant genesets
keep <- grep("NATURAL_KILLER", names(msigdb_subcollection))
msigdb_gobp_nk <- msigdb_subcollection[keep]

usethis::use_data(msigdb_gobp_nk, overwrite = TRUE)
