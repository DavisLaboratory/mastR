test_that("sig_boxplot works", {
  data("im_data_6", "NK_markers")
  ## test eSet object
  ### test expression boxplot
  p <- sig_boxplot(
    im_data_6, sigs = NK_markers$HGNC_Symbol[10:30],
    group_col = "celltype:ch1", target_group = "NK",
    type = "expression",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test score boxplot
  p <- sig_boxplot(
    im_data_6, sigs = NK_markers$HGNC_Symbol[10:30],
    group_col = "celltype:ch1", target_group = "NK",
    type = "score",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  p <- sig_boxplot(
    dge, sigs = NK_markers$HGNC_Symbol[10:30],
    group_col = "group", target_group = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  p <- sig_boxplot(
    data_seurat, sigs = NK_markers$HGNC_Symbol[10:30],
    group_col = "group", target_group = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list objects
  p <- sig_boxplot(
    list(A = im_data_6, B = im_data_6),
    sigs = NK_markers$HGNC_Symbol[10:30],
    group_col = "celltype:ch1", target_group = "NK",
    type = c("score", "expression"),
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
