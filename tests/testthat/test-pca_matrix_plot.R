test_that("pca_matrix_plot works", {
  data("im_data_6", "NK_markers")
  ## test without loadings
  p <- pca_matrix_plot(im_data_6,
    features = NK_markers$HGNC_Symbol,
    group_by = "celltype:ch1", loading = FALSE,
    gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))

  ## test with loadings
  p <- pca_matrix_plot(im_data_6,
    features = NK_markers$HGNC_Symbol,
    group_by = "celltype:ch1", loading = TRUE,
    n_loadings = 3, gene_id = "ENSEMBL", scale = FALSE
  )
  expect_true(ggplot2::is.ggplot(p))
})
