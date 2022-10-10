test_that("merge_markers works", {

  Msig <- get_msigdb_sig(
    species = "Homo sapiens", cat = "C5", subcat = "GO:BP",
    pattern = "natural_killer_cell_mediated"
  )
  Panglao <- get_panglao_sig(species = "Hs", type = "NK cells")
  markers <- merge_markers(markers = list(NK = NK_markers$HGNC_Symbol,
                                          MSigDB = Msig,
                                          PanglaoDB = Panglao))

  expect_s3_class(markers, "tbl")
  expect_setequal(markers$Gene, list(NK_markers$HGNC_Symbol, Msig, Panglao) |>
                    Reduce(f = union))
  ## test un-named list
  expect_error(merge_markers(markers = list(NK_markers, Msig, Panglao)))
})
