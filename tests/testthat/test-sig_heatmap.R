test_that("sig_heatmap works", {
  data("im_data_6", "nk_markers")
  ## test matrix object
  p <- sig_heatmap(
    Biobase::exprs(im_data_6),
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = im_data_6$`celltype:ch1`,
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test eSet object
  ### test expression boxplot
  p <- sig_heatmap(
    im_data_6,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is(p, "Heatmap"))

  ## test DGEList object
  dge <- edgeR::DGEList(
    counts = Biobase::exprs(im_data_6),
    group = im_data_6$`celltype:ch1`
  )
  p <- sig_heatmap(
    dge,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "group",
    gene_id = "ENSEMBL"
  )
  expect_true(is(p, "Heatmap"))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(
    counts = Biobase::exprs(im_data_6),
    meta.data = dge$samples
  )
  p <- sig_heatmap(
    data_seurat,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "group",
    gene_id = "ENSEMBL"
  )
  expect_true(is(p, "Heatmap"))

  ## test list objects
  p <- sig_heatmap(
    list(A = im_data_6, B = im_data_6),
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list sigs
  p <- sig_heatmap(
    im_data_6,
    sigs = list(
      A = nk_markers$HGNC_Symbol[1:10],
      B = nk_markers$HGNC_Symbol[21:40]
    ),
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL",
    scale = "row"
  )
  expect_true(is(p, "Heatmap"))
})
