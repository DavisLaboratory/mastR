test_that("sig_heatmap works", {
  data("im_data_6", "NK_markers")
  ## test matrix object
  p <- sig_heatmap(
    im_data_6@assayData$exprs,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = im_data_6$`celltype:ch1`,
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test eSet object
  ### test expression boxplot
  p <- sig_heatmap(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test score boxplot
  p <- sig_heatmap(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  p <- sig_heatmap(
    dge, sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "group",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  p <- sig_heatmap(
    data_seurat,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "group",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list objects
  p <- sig_heatmap(
    list(A = im_data_6, B = im_data_6),
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
