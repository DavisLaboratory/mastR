msigdb_collection <- msigdb::getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
msigdb_subcollection <- msigdb::subsetCollection(msigdb, collection = "c5", subcollection = "GO:BP")

usethis::use_data(msigdb_subcollection, overwrite = TRUE)
