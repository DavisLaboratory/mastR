test_that("get_lm_sig works", {
  lm7 <- get_lm_sig(lm7.pattern = "\\w+")
  lm22 <- get_lm_sig(lm22.pattern = "\\w+")
  lm <- get_lm_sig(lm7.pattern = "\\w+", lm22.pattern = "\\w+")

  expect_setequal(lm7, LM7$Gene[!is.na(LM7$Subset)])
  expect_setequal(lm22, LM22$Gene)
  expect_setequal(union(lm7, lm22), lm)
})
