test_that("sig_rankdensity_plot works", {
  data("im_data_6", "NK_markers")
  ## test matrix object
  ### test aggregate = FALSE
  p <- sig_rankdensity_plot(
    im_data_6@assayData$exprs,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = im_data_6$`celltype:ch1`,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test aggregate = TRUE
  p <- sig_rankdensity_plot(
    im_data_6@assayData$exprs,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = im_data_6$`celltype:ch1`,
    gene_id = "ENSEMBL",
    aggregate = TRUE
  )
  expect_true(is.ggplot(p))

  ## test eSet object
  ### test expression boxplot
  p <- sig_rankdensity_plot(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test score boxplot
  p <- sig_rankdensity_plot(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  p <- sig_rankdensity_plot(
    dge, sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "group",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  p <- sig_rankdensity_plot(
    data_seurat,
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "group",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list objects
  p <- sig_rankdensity_plot(
    list(A = im_data_6, B = im_data_6),
    sigs = NK_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
