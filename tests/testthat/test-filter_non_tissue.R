test_that("filter_non_tissue works", {
  data("NK_markers", "ccle_crc_5")

  ## test default 'CCLE'
  non_tissue_genes <- filter_non_tissue(
    # data = "CCLE",
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_true(all(non_tissue_genes %in% NK_markers$HGNC_Symbol))

  ## test data.frame object
  non_tissue_genes_1 <- filter_non_tissue(
    as.data.frame(ccle_crc_5$counts),
    ID = ccle_crc_5$samples$cancer,
    type = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_true(all(non_tissue_genes_1 %in% NK_markers$HGNC_Symbol))

  ## test matrix object
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5$counts,
    ID = ccle_crc_5$samples$cancer,
    type = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test Matrix object
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5$counts |> Matrix::Matrix(),
    ID = ccle_crc_5$samples$cancer,
    type = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test DGEList object
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5,
    ID = "cancer",
    type = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(
    counts = ccle_crc_5$counts,
    meta.data = ccle_crc_5$samples,
    markers = NK_markers$HGNC_Symbol
  )
  non_tissue_genes_2 <- filter_non_tissue(
    data_seurat,
    ID = "cancer",
    type = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  non_tissue_genes_2 <- filter_non_tissue(
    data_sce,
    ID = "cancer",
    type = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test without providing markers
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5,
    ID = "cancer",
    type = "CRC"
  )
  expect_true(all(non_tissue_genes_1 %in% union(non_tissue_genes_2,
                                                NK_markers$HGNC_Symbol)))
})
