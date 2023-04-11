test_that("sig_rankdensity_plot works", {
  data("im_data_6", "nk_markers")
  ## test matrix object
  ### test aggregate = FALSE
  p <- sig_rankdensity_plot(
    Biobase::exprs(im_data_6),
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = im_data_6$`celltype:ch1`,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test aggregate = TRUE
  p <- sig_rankdensity_plot(
    Biobase::exprs(im_data_6),
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = im_data_6$`celltype:ch1`,
    gene_id = "ENSEMBL",
    aggregate = TRUE
  )
  expect_true(is.ggplot(p))

  ## test eSet object
  ### test expression boxplot
  p <- sig_rankdensity_plot(
    im_data_6,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test score boxplot
  p <- sig_rankdensity_plot(
    im_data_6,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(
    counts = Biobase::exprs(im_data_6),
    group = im_data_6$`celltype:ch1`
  )
  p <- sig_rankdensity_plot(
    dge,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "group",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- SeuratObject::CreateSeuratObject(
    counts = Biobase::exprs(im_data_6),
    meta.data = dge$samples
  )
  p <- sig_rankdensity_plot(
    data_seurat,
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "group",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list objects
  p <- sig_rankdensity_plot(
    list(A = im_data_6, B = im_data_6),
    sigs = nk_markers$HGNC_Symbol[1:20],
    group_col = "celltype:ch1",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
