test_that("merge_markers works", {

  Msig <- get_MSig(
    species = "Homo sapiens", cat = "C5", subcat = "GO:BP",
    pattern = "natural_killer_cell_mediated"
  )
  Panglao <- get_Panglao(species = "Hs", type = "NK cells")
  markers <- merge_markers(markers = list(NK_markers, Msig, Panglao))

  expect_s3_class(markers, "tbl")
  expect_setequal(markers$HGNC_Symbol, list(NK_markers$HGNC_Symbol,
                                            Msig$HGNC_Symbol,
                                            Panglao$HGNC_Symbol) |> Reduce(f = union))
})
