test_that("get_DE_table works", {
  ## test Expression object
  DEG_table <- get_DE_table(im_data_6, ID = "celltype:ch1", type = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test DGEList object
  data <- proc_data
  DEG_table <- get_DE_table(data, ID = "celltype.ch1",
                            type = "NK", counts = FALSE)
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test matrix object
  DEG_table <- get_DE_table(im_data_6@assayData$exprs,
                            ID = im_data_6$`celltype:ch1`,
                            type = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test Matrix object
  DEG_table <- get_DE_table(im_data_6@assayData$exprs |> Matrix::Matrix(),
                            ID = im_data_6$`celltype:ch1`,
                            type = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = data$samples)
  DEG_table <- get_DE_table(data_seurat,
                            ID = "celltype.ch1",
                            type = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  DEG_table <- get_DE_table(data_sce,
                            ID = "celltype.ch1",
                            type = "NK")
  expect_true(all(sapply(DEG_table, is.data.frame)))
})
