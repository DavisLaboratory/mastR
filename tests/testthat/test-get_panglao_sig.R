test_that("get_panglao_sig works", {
  Panglao <- get_panglao_sig(type = "NK cells")

  expect_type(Panglao, "character")
  expect_error(get_panglao_sig(type = "unknow"))
  expect_error(get_panglao_sig(species = "unknow"))
})
