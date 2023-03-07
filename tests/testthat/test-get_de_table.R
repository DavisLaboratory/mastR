test_that("get_de_table works", {
  data("im_data_6")
  ## test Expression object
  DEG_table <- get_de_table(im_data_6, group_col = "celltype:ch1",
                            target_group = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test DGEList object
  data <- process_data(im_data_6, group_col = "celltype:ch1",
                       target_group = "NK")
  DEG_table <- get_de_table(data, group_col = "celltype.ch1",
                            target_group = "NK", normalize = FALSE)
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test matrix object
  DEG_table <- get_de_table(im_data_6@assayData$exprs,
                            group_col = im_data_6$`celltype:ch1`,
                            target_group = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test Matrix object
  DEG_table <- get_de_table(im_data_6@assayData$exprs |> Matrix::Matrix(),
                            group_col = im_data_6$`celltype:ch1`,
                            target_group = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))


  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = data$samples)
  DEG_table <- get_de_table(data_seurat,
                            group_col = "celltype.ch1",
                            target_group = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  DEG_table <- get_de_table(data_sce,
                            group_col = "celltype.ch1",
                            target_group = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))
})
