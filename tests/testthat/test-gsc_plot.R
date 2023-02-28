test_that("gsc_plot works", {
  data("msigdb_gobp_nk")
  ## test on gsc
  p <- gsc_plot(msigdb_gobp_nk[1:3])
  expect_true(is(p, "upset"))

  ## test on gs
  p <- gsc_plot(msigdb_gobp_nk[[1]], msigdb_gobp_nk[[2]])
  expect_true(is(p, "upset"))

  ## test on mixed gs and gsc
  p <- gsc_plot(msigdb_gobp_nk[[1]], msigdb_gobp_nk[[2]], msigdb_gobp_nk[3:4])
  expect_true(is(p, "upset"))
})
