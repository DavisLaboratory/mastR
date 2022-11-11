test_that("filter_non_tissue works", {
  ## load CCLE
  CCLE <- depmap::depmap_TPM()
  CCLE <- tidyr::pivot_wider(
    CCLE[, c("gene_name", "cell_line", "rna_expression")],
    names_from = "cell_line",
    values_from = "rna_expression"
  ) |> data.frame()
  rownames(CCLE) <- CCLE$gene_name
  CCLE$gene_name <- NULL
  CCLE_meta <- depmap::depmap_metadata()
  cell_lines <- CCLE_meta$primary_disease[match(
    colnames(CCLE),
    CCLE_meta$cell_line
  )]
  rm(CCLE_meta)

  ## test CCLE
  non_tissue_genes <- filter_non_tissue(
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_true(all(non_tissue_genes %in% NK_markers$HGNC_Symbol))

  ## test data.frame object
  non_tissue_genes_2 <- filter_non_tissue(CCLE,
    ID = cell_lines,
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes, non_tissue_genes_2)

  ## test matrix object
  non_tissue_genes_2 <- filter_non_tissue(CCLE |> as.matrix(),
    ID = cell_lines,
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes, non_tissue_genes_2)

  ## test Matrix object
  non_tissue_genes_2 <- filter_non_tissue(CCLE |> as.matrix() |> Matrix::Matrix(),
    ID = cell_lines,
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes, non_tissue_genes_2)

  ## test DGEList object
  DGE <- edgeR::DGEList(counts = CCLE, group = cell_lines)
  non_tissue_genes_2 <- filter_non_tissue(DGE,
    ID = "group",
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes, non_tissue_genes_2)

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(
    counts = DGE$counts,
    meta.data = DGE$samples,
    markers = NK_markers$HGNC_Symbol
  )
  non_tissue_genes_2 <- filter_non_tissue(data_seurat,
    ID = "group",
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes, non_tissue_genes_2)

  ## test sce object
  data_sce <- Seurat::as.SingleCellExperiment(data_seurat)
  DEGs <- filter_non_tissue(data_sce,
    ID = "group",
    type = "colorectal",
    markers = NK_markers$HGNC_Symbol
  )
  expect_setequal(non_tissue_genes, non_tissue_genes_2)

  ## test without providing markers
  non_tissue_genes_2 <- filter_non_tissue(CCLE,
    ID = cell_lines,
    type = "colorectal"
  )
  expect_true(all(non_tissue_genes %in% union(non_tissue_genes_2,
                                              NK_markers$HGNC_Symbol)))
})
