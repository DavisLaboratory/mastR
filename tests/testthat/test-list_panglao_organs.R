test_that("list_panglao_organs works", {
  av_organs <- list_panglao_organs()
  expect_type(av_organs, "character")
})
