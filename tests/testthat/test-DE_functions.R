test_that("filterGenes works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  expect_type(filterGenes(dge = dge, group_col = "group"), 'logical')
  expect_type(filterGenes(dge = dge, group_col = "group", normalize = FALSE),
              'logical')
})

test_that("voom_fit_treat works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  tfit_RP <- voom_fit_treat(dge = dge, group_col = "group",
                            target_group = "NK", feature_selection = "rankproduct")
  tfit_GP <- voom_fit_treat(dge = dge, group_col = "group",
                            target_group = "NK", feature_selection = "none",
                            group = TRUE)

  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6) + 1,
                        group = im_data_6$`celltype:ch1`)

  tfit_GP_c <- voom_fit_treat(dge = dge, group_col = "group",
                              target_group = "NK", feature_selection = "none",
                              group = TRUE, normalize = FALSE)

  expect_true(is(tfit_RP$tfit, "MArrayLM"))
  expect_true(is(tfit_GP$tfit, "MArrayLM"))
  expect_true(is(tfit_GP_c$tfit, "MArrayLM"))
  expect_true(is(tfit_RP, "DGEList"))
  expect_true(is(tfit_GP, "DGEList"))
  expect_true(is(tfit_GP_c, "DGEList"))
})

test_that("DEGs_RP works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)
  tfit <- voom_fit_treat(dge = dge, group_col = "group",
                         target_group = "NK",
                         feature_selection = "rankproduct")$tfit

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
  tfit <- voom_fit_treat(dge = dge, group_col = "group",
                         target_group = "NK",
                         feature_selection = "none")$tfit

  DEGs_i <- DEGs_Group(tfit = tfit, assemble = "intersect")
  DEGs_u <- DEGs_Group(tfit = tfit, assemble = "union")

  expect_true(is(DEGs_i, "list"))
  expect_true(is(DEGs_u, "list"))
  expect_true(all(DEGs_i$UP %in% DEGs_u$UP))
  expect_true(all(DEGs_i$DOWN %in% DEGs_u$DOWN))
})
