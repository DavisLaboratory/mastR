test_that("pca_matrix_plot works", {
  data("im_data_6", "nk_markers")
  ## test without loadings
  p <- pca_matrix_plot(
    im_data_6,
    features = nk_markers$HGNC_Symbol,
    group_by = "celltype:ch1", loading = FALSE,
    gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))

  ## test with loadings
  p <- pca_matrix_plot(
    im_data_6,
    features = nk_markers$HGNC_Symbol,
    group_by = "celltype:ch1", loading = TRUE,
    n_loadings = 3, gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))

  ## test matrix
  p <- pca_matrix_plot(
    Biobase::exprs(im_data_6),
    features = nk_markers$HGNC_Symbol,
    group_by = im_data_6$`celltype:ch1`,
    gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(
    counts = Biobase::exprs(im_data_6),
    samples = Biobase::pData(im_data_6)
  )
  p <- pca_matrix_plot(
    dge,
    features = nk_markers$HGNC_Symbol,
    group_by = "celltype.ch1",
    gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))

  ## test seurat object
  data_seurat <- SeuratObject::CreateSeuratObject(
    counts = Biobase::exprs(im_data_6),
    meta.data = dge$samples
  )
  p <- pca_matrix_plot(
    data_seurat,
    features = nk_markers$HGNC_Symbol,
    group_by = "celltype.ch1",
    gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))
})
