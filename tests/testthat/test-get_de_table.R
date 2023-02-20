test_that("get_de_table works", {
  data("im_data_6")
  ## test Expression object
  DEG_table <- get_de_table(im_data_6, ID = "celltype:ch1", type = "NK")
  expect_true(all(sapply(DEG_table[-length(DEG_table)], is.data.frame)))
  expect_true(is(DEG_table$proc_data, "DGEList"))

  ## test DGEList object
  data <- DEG_table$proc_data
  DEG_table <- get_de_table(data, ID = "celltype.ch1",
                            type = "NK", counts = FALSE)
  expect_true(all(sapply(DEG_table[-length(DEG_table)], is.data.frame)))
  expect_true(is(DEG_table$proc_data, "DGEList"))

  ## test matrix object
  DEG_table <- get_de_table(im_data_6@assayData$exprs,
                            ID = im_data_6$`celltype:ch1`,
                            type = "NK")
  expect_true(all(sapply(DEG_table[-length(DEG_table)], is.data.frame)))
  expect_true(is(DEG_table$proc_data, "DGEList"))

  ## test Matrix object
  DEG_table <- get_de_table(im_data_6@assayData$exprs |> Matrix::Matrix(),
                            ID = im_data_6$`celltype:ch1`,
                            type = "NK")
  expect_true(all(sapply(DEG_table[-length(DEG_table)], is.data.frame)))
  expect_true(is(DEG_table$proc_data, "DGEList"))


  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = data$samples)
  DEG_table <- get_de_table(data_seurat,
                            ID = "celltype.ch1",
                            type = "NK")
  expect_true(all(sapply(DEG_table[-length(DEG_table)], is.data.frame)))
  expect_true(is(DEG_table$proc_data, "DGEList"))

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  DEG_table <- get_de_table(data_sce,
                            ID = "celltype.ch1",
                            type = "NK")
  expect_true(all(sapply(DEG_table[-length(DEG_table)], is.data.frame)))
  expect_true(is(DEG_table$proc_data, "DGEList"))
})
