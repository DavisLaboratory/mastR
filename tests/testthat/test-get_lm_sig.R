test_that("get_lm_sig works", {
  data("lm7", "lm22")
  LM7 <- get_lm_sig(lm7.pattern = "\\w+")
  LM22 <- get_lm_sig(lm22.pattern = "\\w+")
  LM <- get_lm_sig(lm7.pattern = "\\w+", lm22.pattern = "\\w+")

  expect_setequal(GSEABase::geneIds(LM7), lm7$Gene[!is.na(lm7$Subset)])
  expect_setequal(GSEABase::geneIds(LM22), lm22$Gene)
  expect_setequal(GSEABase::geneIds(LM7 | LM22), Reduce(f = union, GSEABase::geneIds(LM)))
  expect_error(get_lm_sig())
})
