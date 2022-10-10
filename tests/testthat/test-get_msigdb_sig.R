test_that("get_msigdb_sig works", {
  MSig <- get_msigdb_sig(pattern = "natural_killer_cell_mediated")

  expect_type(MSig, "character")
})
