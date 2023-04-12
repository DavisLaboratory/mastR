test_that("list_panglao_types works", {
  av_types <- list_panglao_types(organ = "Immune system")
  expect_type(av_types, "character")
  av_types <- list_panglao_types(organ = "Brain")
  expect_type(av_types, "character")
  av_types <- list_panglao_types(organ = "Liver")
  expect_type(av_types, "character")

  expect_error(list_panglao_types(organ = "unknow"))
  expect_error(list_panglao_types(organ = c("Immune system", "Brain")))
})
