test_that("gls2gsc works", {
  data("msigdb_gobp_nk")
  ## test geneset list
  gsc1 <- gls2gsc(GSEABase::geneIds(msigdb_gobp_nk[1:3]))
  expect_setequal(GSEABase::geneIds(gsc1),
                  GSEABase::geneIds(msigdb_gobp_nk[1:3]))
  ## test vector
  gsc2 <- gls2gsc(msigdb_gobp_nk[[1]]@geneIds,
                  msigdb_gobp_nk[[2]]@geneIds,
                  msigdb_gobp_nk[[3]]@geneIds)
  expect_setequal(GSEABase::geneIds(gsc1), GSEABase::geneIds(gsc2))
})
