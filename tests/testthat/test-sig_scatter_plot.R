test_that("multiplication works", {
  data("NK_markers")
  data("im_data_6", "NK_markers")
  ## test eSet object
  p <- sig_scatter_plot(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol,
    group_col = "celltype:ch1", target_group = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(
    counts = Biobase::exprs(im_data_6),
    group = im_data_6$`celltype:ch1`
  )
  p <- sig_scatter_plot(
    dge,
    sigs = NK_markers$HGNC_Symbol,
    group_col = "group", target_group = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(
    counts = Biobase::exprs(im_data_6),
    meta.data = dge$samples
  )
  p <- sig_scatter_plot(
    data_seurat,
    sigs = NK_markers$HGNC_Symbol,
    group_col = "group", target_group = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
