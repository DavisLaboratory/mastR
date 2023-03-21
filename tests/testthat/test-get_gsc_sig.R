test_that("get_gsc_sig works", {
  MSig <- get_gsc_sig(
    gsc = mastR::msigdb_gobp_nk,
    pattern = "natural_killer_cell_mediated",
    ignore.case = TRUE
  )

  expect_s4_class(MSig, "GeneSetCollection")
  expect_error(get_gsc_sig(
    gsc = mastR::msigdb_gobp_nk,
    pattern = "natural_killer_cell_mediated"
  ))
  expect_warning(get_gsc_sig(
    gsc = mastR::msigdb_gobp_nk,
    pattern = c(
      "NATURAL_KILLER_CELL_MEDIATED",
      "T_CELL"
    )
  ))
})
