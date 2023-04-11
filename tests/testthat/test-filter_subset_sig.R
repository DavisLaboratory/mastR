test_that("filter_subset_sig works", {
  data("im_data_6", "nk_markers")
  ## test matrix object
  sig <- filter_subset_sig(
    Biobase::exprs(im_data_6),
    group_col = im_data_6$`celltype:ch1`,
    target_group = "NK",
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% nk_markers$HGNC_Symbol))

  ## test eSet object
  sig <- filter_subset_sig(
    im_data_6,
    group_col = "celltype:ch1", target_group = "NK",
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% nk_markers$HGNC_Symbol))

  ## test DGEList object
  dge <- edgeR::DGEList(
    counts = Biobase::exprs(im_data_6),
    group = im_data_6$`celltype:ch1`
  )
  sig <- filter_subset_sig(
    dge,
    group_col = "group", target_group = "NK",
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% nk_markers$HGNC_Symbol))

  ## test seurat object
  data_seurat <- SeuratObject::CreateSeuratObject(
    counts = Biobase::exprs(im_data_6),
    meta.data = dge$samples
  )
  sig <- filter_subset_sig(
    data_seurat,
    group_col = "group", target_group = "NK",
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% nk_markers$HGNC_Symbol))

  ## test list objects
  p <- filter_subset_sig(
    list(A = im_data_6, B = im_data_6),
    group_col = "celltype:ch1", target_group = "NK",
    markers = nk_markers$HGNC_Symbol,
    gene_id = "ENSEMBL"
  )
  expect_true(all(sig %in% nk_markers$HGNC_Symbol))
})
