test_that("get_degs works", {
  ## test Expression object
  DEGs <- get_degs(im_data_6, ID = "celltype:ch1", type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN", "proc_data"))

  ## test DGEList object
  data <- DEGs$proc_data
  DEGs <- get_degs(data, ID = "celltype.ch1",
                   type = "NK", counts = FALSE)
  expect_identical(names(DEGs), c("UP", "DOWN", "proc_data"))

  ## test matrix object
  DEGs <- get_degs(im_data_6@assayData$exprs,
                   ID = im_data_6$`celltype:ch1`,
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN", "proc_data"))

  ## test Matrix object
  DEGs <- get_degs(im_data_6@assayData$exprs |> Matrix::Matrix(),
                   ID = im_data_6$`celltype:ch1`,
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN", "proc_data"))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = data$samples)
  DEGs <- get_degs(data_seurat,
                   ID = "celltype.ch1",
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN", "proc_data"))

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  DEGs <- get_degs(data_sce,
                   ID = "celltype.ch1",
                   type = "NK")
  expect_identical(names(DEGs), c("UP", "DOWN", "proc_data"))
})
