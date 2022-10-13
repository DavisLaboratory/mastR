test_that("filter_subset_sig works", {
  ## test matrix object
  sig <- filter_subset_sig(
    im_data_6@assayData$exprs,
    ID = im_data_6$`celltype:ch1`,
    type = "NK",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% NK_markers$HGNC_Symbol))

  ## test eSet object
  sig <- filter_subset_sig(
    im_data_6,
    ID = "celltype:ch1", type = "NK",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% NK_markers$HGNC_Symbol))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  sig <- filter_subset_sig(
    dge,
    ID = "group", type = "NK",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% NK_markers$HGNC_Symbol))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  sig <- filter_subset_sig(
    data_seurat,
    ID = "group", type = "NK",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% NK_markers$HGNC_Symbol))

  ## test list objects
  p <- filter_subset_sig(
    list(A = im_data_6, B = im_data_6),
    ID = "celltype:ch1", type = "NK",
    markers = NK_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% NK_markers$HGNC_Symbol))
})
