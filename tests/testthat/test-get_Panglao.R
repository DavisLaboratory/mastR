test_that("get_Panglao works", {
  Panglao <- get_Panglao(type = "NK cells")

  expect_s3_class(Panglao, "tbl")
  expect_error(get_Panglao(type = "unknow"))
  expect_error(get_Panglao(species = "unknow"))
})
