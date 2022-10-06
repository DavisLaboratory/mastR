test_that("filterGenes works", {
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  expect_type(filterGenes(dge = dge, ID = "group"), 'logical')
  expect_type(filterGenes(dge = dge, ID = "group", counts = FALSE), 'logical')
})

test_that("voom_lm_fit works", {
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

  expect_true(is(tfit_RP, "MArrayLM"))
  expect_true(is(tfit_GP, "MArrayLM"))
  expect_true(is(tfit_GP_c, "MArrayLM"))
})

test_that("DEGs_RP works", {
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)
  tfit <- voom_lm_fit(dge = dge, ID = "group",
                      type = "NK", method = "RP")

  DEGs_i <- DEGs_RP(tfit = tfit, assemble = "intersect")
  DEGs_u <- DEGs_RP(tfit = tfit, assemble = "union")

  expect_true(is(DEGs_i, "list"))
  expect_true(is(DEGs_u, "list"))
  expect_identical(names(DEGs_i), c("UP", "DOWN"))
  expect_identical(names(DEGs_u), c("UP", "DOWN"))
})

test_that("DEGs_Group works", {
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)
  tfit <- voom_lm_fit(dge = dge, ID = "group",
                      type = "NK", method = "Group")

  DEGs_i <- DEGs_Group(tfit = tfit, assemble = "intersect")
  DEGs_u <- DEGs_Group(tfit = tfit, assemble = "union")

  expect_true(is(DEGs_i, "list"))
  expect_true(is(DEGs_u, "list"))
  expect_true(all(DEGs_i$UP %in% DEGs_u$UP))
  expect_true(all(DEGs_i$DOWN %in% DEGs_u$DOWN))
})
