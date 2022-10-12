test_that("multiplication works", {
  ## test eSet object
  p <- sig_scatter_plot(
    im_data_6, sigs = NK_markers$HGNC_Symbol,
    ID = "celltype:ch1", type = "NK",
    gene_id = "ENSEMBL"
    )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  p <- sig_scatter_plot(
    dge, sigs = NK_markers$HGNC_Symbol,
    ID = "group", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  p <- sig_scatter_plot(
    data_seurat, sigs = NK_markers$HGNC_Symbol,
    ID = "group", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
