test_that("get_DEGs works", {
  ## test Expression object
  DEGs <- get_DEGs(im_data_6, ID = "celltype:ch1", type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN"))

  ## test DGEList object
  data <- proc_data
  DEGs <- get_DEGs(data, ID = "celltype.ch1",
                   type = "NK", counts = FALSE)
  expect_identical(names(DEGs), c("UP", "DOWN"))

  ## test matrix object
  DEGs <- get_DEGs(im_data_6@assayData$exprs,
                   ID = im_data_6$`celltype:ch1`,
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN"))

  ## test Matrix object
  DEGs <- get_DEGs(im_data_6@assayData$exprs |> Matrix::Matrix(),
                   ID = im_data_6$`celltype:ch1`,
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN"))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = data$samples)
  DEGs <- get_DEGs(data_seurat,
                   ID = "celltype.ch1",
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN"))

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  DEGs <- get_DEGs(data_sce,
                   ID = "celltype.ch1",
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN"))
})
