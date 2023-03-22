test_that("get_lm_sig works", {
  data("LM7", "LM22")
  lm7 <- get_lm_sig(lm7.pattern = "\\w+")
  lm22 <- get_lm_sig(lm22.pattern = "\\w+")
  lm <- get_lm_sig(lm7.pattern = "\\w+", lm22.pattern = "\\w+")

  expect_setequal(GSEABase::geneIds(lm7), LM7$Gene[!is.na(LM7$Subset)])
  expect_setequal(GSEABase::geneIds(lm22), LM22$Gene)
  expect_setequal(GSEABase::geneIds(lm7 | lm22), Reduce(f = union, GSEABase::geneIds(lm)))
  expect_error(get_lm_sig())
})
