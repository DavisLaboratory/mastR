test_that("merge_markers works", {

  data("msigdb_gobp_nk")

  ## test GeneSetCollection
  markers <- merge_markers(msigdb_gobp_nk)

  expect_s4_class(markers, "GeneSet")
  expect_setequal(markers@geneIds,
                  Reduce(f = union, GSEABase::geneIds(msigdb_gobp_nk)))

  ## test GeneSet
  markers <- merge_markers(msigdb_gobp_nk[[1]],
                           msigdb_gobp_nk[[2]])

  expect_s4_class(markers, "GeneSet")
  expect_setequal(markers@geneIds,
                  Reduce(f = union, GSEABase::geneIds(msigdb_gobp_nk[1:2])))

  ## test on mixed GeneSet and GeneSetCollection
  markers <- merge_markers(msigdb_gobp_nk[[1]], msigdb_gobp_nk[2:3])
  expect_s4_class(markers, "GeneSet")
  expect_setequal(markers@geneIds,
                  Reduce(f = union, GSEABase::geneIds(msigdb_gobp_nk[1:3])))

  ## test non- GeneSet or GeneSetCollection
  expect_error(merge_markers(msigdb_gobp_nk[[1]]@geneIds))
})
