test_that("filter_non_tissue works", {
  data("NK_markers", "ccle_crc_5")

  ## test default 'CCLE'
  non_tissue_genes <- filter_non_tissue(
    # data = "CCLE",
    group_col = "primary_disease",
    target_group = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_true(all(non_tissue_genes %in% NK_markers$HGNC_Symbol))

  ## test data.frame object
  non_tissue_genes_1 <- filter_non_tissue(
    as.data.frame(ccle_crc_5$counts),
    group_col = ccle_crc_5$samples$cancer,
    target_group = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_true(all(non_tissue_genes_1 %in% NK_markers$HGNC_Symbol))

  ## test matrix object
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5$counts,
    group_col = ccle_crc_5$samples$cancer,
    target_group = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test Matrix object
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5$counts |> Matrix::Matrix(),
    group_col = ccle_crc_5$samples$cancer,
    target_group = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test DGEList object
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5,
    group_col = "cancer",
    target_group = "CRC",
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
    group_col = "cancer",
    target_group = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  non_tissue_genes_2 <- filter_non_tissue(
    data_sce,
    group_col = "cancer",
    target_group = "CRC",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes_1, non_tissue_genes_2)

  ## test without providing markers
  non_tissue_genes_2 <- filter_non_tissue(
    ccle_crc_5,
    group_col = "cancer",
    target_group = "CRC"
  )
  expect_true(all(non_tissue_genes_1 %in% union(non_tissue_genes_2,
                                                NK_markers$HGNC_Symbol)))
})
