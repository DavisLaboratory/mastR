test_that("sig_boxplot works", {
  ## test eSet object
  ### test expression boxplot
  p <- sig_boxplot(
    im_data_6, sigs = NK_markers$HGNC_Symbol[10:30],
    ID = "celltype:ch1", type = "NK",
    plot.score = FALSE,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test score boxplot
  p <- sig_boxplot(
    im_data_6, sigs = NK_markers$HGNC_Symbol[10:30],
    ID = "celltype:ch1", type = "NK",
    plot.score = TRUE,
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  p <- sig_boxplot(
    dge, sigs = NK_markers$HGNC_Symbol[10:30],
    ID = "group", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  p <- sig_boxplot(
    data_seurat, sigs = NK_markers$HGNC_Symbol[10:30],
    ID = "group", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list objects
  p <- sig_boxplot(
    list(A = im_data_6, B = im_data_6),
    sigs = NK_markers$HGNC_Symbol[10:30],
    ID = "celltype:ch1", type = "NK",
    plot.score = c(TRUE, FALSE),
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
