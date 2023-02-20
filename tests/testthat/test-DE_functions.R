test_that("filterGenes works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  expect_type(filterGenes(dge = dge, ID = "group"), 'logical')
  expect_type(filterGenes(dge = dge, ID = "group", counts = FALSE), 'logical')
})

test_that("voom_lm_fit works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  tfit_RP <- voom_lm_fit(dge = dge, ID = "group",
                      type = "NK", method = "RP")
  tfit_GP <- voom_lm_fit(dge = dge, ID = "group",
                         type = "NK", method = "Group",
                         group = TRUE)

  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6) + 1,
                        group = im_data_6$`celltype:ch1`)

  tfit_GP_c <- voom_lm_fit(dge = dge, ID = "group",
                           type = "NK", method = "Group",
                           group = TRUE, counts = FALSE)

  expect_true(is(tfit_RP$tfit, "MArrayLM"))
  expect_true(is(tfit_GP$tfit, "MArrayLM"))
  expect_true(is(tfit_GP_c$tfit, "MArrayLM"))
  expect_true(is(tfit_RP$proc_data, "DGEList"))
  expect_true(is(tfit_GP$proc_data, "DGEList"))
  expect_true(is(tfit_GP_c$proc_data, "DGEList"))
})

test_that("DEGs_RP works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)
  tfit <- voom_lm_fit(dge = dge, ID = "group",
                      type = "NK", method = "RP")$tfit

  DEGs_i <- DEGs_RP(tfit = tfit, assemble = "intersect")
  DEGs_u <- DEGs_RP(tfit = tfit, assemble = "union")

  expect_true(is(DEGs_i, "list"))
  expect_true(is(DEGs_u, "list"))
  expect_identical(names(DEGs_i), c("UP", "DOWN"))
  expect_identical(names(DEGs_u), c("UP", "DOWN"))
})

test_that("DEGs_Group works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)
  tfit <- voom_lm_fit(dge = dge, ID = "group",
                      type = "NK", method = "Group")$tfit

  DEGs_i <- DEGs_Group(tfit = tfit, assemble = "intersect")
  DEGs_u <- DEGs_Group(tfit = tfit, assemble = "union")

  expect_true(is(DEGs_i, "list"))
  expect_true(is(DEGs_u, "list"))
  expect_true(all(DEGs_i$UP %in% DEGs_u$UP))
  expect_true(all(DEGs_i$DOWN %in% DEGs_u$DOWN))
})
