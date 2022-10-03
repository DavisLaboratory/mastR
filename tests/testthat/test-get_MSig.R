test_that("get_MSig works", {
  MSig <- get_MSig(pattern = "natural_killer_cell_mediated")

  expect_s3_class(MSig, "tbl")
})
