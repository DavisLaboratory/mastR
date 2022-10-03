test_that("get_LM works", {
  lm7 <- get_LM(lm7.pattern = "\\w+")
  lm22 <- get_LM(lm22.pattern = "\\w+")
  lm <- get_LM(lm7.pattern = "\\w+", lm22.pattern = "\\w+")

  expect_setequal(lm7$HGNC_Symbol, LM7$Gene[!is.na(LM7$Subset)])
  expect_setequal(lm22$HGNC_Symbol, LM22$Gene)
  expect_setequal(union(lm7$HGNC_Symbol, lm22$HGNC_Symbol), lm$HGNC_Symbol)
})
