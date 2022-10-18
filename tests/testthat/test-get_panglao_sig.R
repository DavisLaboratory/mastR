test_that("get_panglao_sig works", {
  Panglao <- get_panglao_sig(type = "NK cells")

  expect_s4_class(Panglao, "GeneSet")
  expect_error(get_panglao_sig(type = "unknow"))
  expect_error(get_panglao_sig(species = "unknow"))
})
