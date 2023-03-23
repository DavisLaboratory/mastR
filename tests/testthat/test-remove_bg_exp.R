test_that("remove_bg_exp works", {
  data("nk_markers", "ccle_crc_5", "im_data_6")
  ## test matrix
  m1 <- remove_bg_exp(
    bg_data = ccle_crc_5$counts,
    sig_data = Biobase::exprs(im_data_6),
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = c("ENSEMBL", "SYMBOL")
  )
  expect_true(is.vector(m1))

  m1 <- remove_bg_exp(
    bg_data = ccle_crc_5$counts,
    sig_data = Biobase::exprs(im_data_6),
    b_group_col = ccle_crc_5$samples$cancer,
    b_target_group = "CRC",
    s_group_col = im_data_6$`celltype:ch1`,
    s_target_group = "NK",
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = c("ENSEMBL", "SYMBOL")
  )
  expect_true(is.vector(m1))

  ## test DGEList
  m2 <- remove_bg_exp(
    bg_data = ccle_crc_5,
    sig_data = Biobase::exprs(im_data_6),
    b_group_col = "cancer",
    b_target_group = "CRC",
    s_group_col = im_data_6$`celltype:ch1`,
    s_target_group = "NK",
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = c("ENSEMBL", "SYMBOL")
  )
  expect_true(is.vector(m2))

  expect_setequal(m1, m2)

  m <- remove_bg_exp(
    bg_data = Biobase::exprs(im_data_6),
    sig_data= ccle_crc_5,
    b_group_col = im_data_6$`celltype:ch1`,
    b_target_group = "NK",
    s_group_col = "cancer",
    s_target_group = "CRC",
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = c("SYMBOL", "ENSEMBL")
  )
  expect_true(is.vector(m))

  ## test eSet
  m <- remove_bg_exp(
    bg_data = im_data_6,
    sig_data = Biobase::exprs(im_data_6),
    b_group_col = "celltype:ch1",
    b_target_group = "CD4",
    s_group_col = im_data_6$`celltype:ch1`,
    s_target_group = "NK",
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = "ENSEMBL"
  )
  expect_true(is.vector(m))

  ## test seurat
  data_seurat <- Seurat::CreateSeuratObject(
    counts = ccle_crc_5$counts,
    meta.data = ccle_crc_5$samples
  )
  m <- remove_bg_exp(
    bg_data = data_seurat,
    sig_data = Biobase::exprs(im_data_6),
    b_group_col = "cancer",
    b_target_group = "CRC",
    s_group_col = im_data_6$`celltype:ch1`,
    s_target_group = "NK",
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = c("ENSEMBL", "SYMBOL")
  )
  expect_true(is.vector(m))

  ## test 'CCLE'
  m <- remove_bg_exp(
    # "CCLE",
    sig_data = Biobase::exprs(im_data_6),
    b_group_col = "primary_disease",
    b_target_group = "colorectal",
    s_group_col = im_data_6$`celltype:ch1`,
    s_target_group = "NK",
    markers = nk_markers$HGNC_Symbol[40:50],
    gene_id = c("ENSEMBL", "SYMBOL"),
    filter = c(1, 10),
    ignore.case = TRUE
  )
  expect_true(is.vector(m))
})
