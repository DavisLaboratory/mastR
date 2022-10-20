test_that("merge_markers works", {

  Msig <- get_msigdb_sig(
    pattern = "natural_killer_cell_mediated",
    ignore.case = TRUE
  )
  Panglao <- get_panglao_sig(species = "Hs", type = "NK cells")

  ## test GeneSetCollection
  gsc <- GSEABase::GeneSetCollection(list(NK_markers$HGNC_Symbol |>
                                            GSEABase::GeneSet(setName = "NK_curson",
                                                              geneIdType = GSEABase::SymbolIdentifier()),
                                          Msig,
                                          Panglao))
  markers <- merge_markers(markers = gsc)

  expect_s4_class(markers, "GeneSet")
  expect_setequal(markers@geneIds, list(NK_markers$HGNC_Symbol,
                                        Msig@geneIds,
                                        Panglao@geneIds) |>
                    Reduce(f = union))

  ## test markers vector list
  markers <- merge_markers(markers = list(NK = NK_markers$HGNC_Symbol,
                                          MSigDB = Msig@geneIds,
                                          PanglaoDB = Panglao@geneIds))

  expect_s4_class(markers, "GeneSet")
  expect_setequal(markers@geneIds, list(NK_markers$HGNC_Symbol,
                                        Msig@geneIds,
                                        Panglao@geneIds) |>
                    Reduce(f = union))
  ## test un-named list
  expect_error(merge_markers(markers = list(NK_markers, Msig@geneIds, Panglao@geneIds)))
})
