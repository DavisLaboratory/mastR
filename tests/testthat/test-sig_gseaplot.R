test_that("sig_gseaplot works", {
  ## test matrix object
  p <- sig_gseaplot(
    im_data_6@assayData$exprs,
    sigs = NK_markers$HGNC_Symbol[1:20],
    ID = im_data_6$`celltype:ch1`,
    type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test eSet object
  ### test expression boxplot
  p <- sig_gseaplot(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol[1:20],
    ID = "celltype:ch1", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ### test score boxplot
  p <- sig_gseaplot(
    im_data_6,
    sigs = NK_markers$HGNC_Symbol[1:20],
    ID = "celltype:ch1", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(counts = im_data_6@assayData$exprs,
                        group = im_data_6$`celltype:ch1`)
  p <- sig_gseaplot(
    dge, sigs = NK_markers$HGNC_Symbol[1:20],
    ID = "group", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(counts = im_data_6@assayData$exprs,
                                            meta.data = dge$samples)
  p <- sig_gseaplot(
    data_seurat,
    sigs = NK_markers$HGNC_Symbol[1:20],
    ID = "group", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))

  ## test list objects with list of signatures
  p <- sig_gseaplot(
    list(A = im_data_6, B = im_data_6),
    sigs = list(a = NK_markers$HGNC_Symbol[1:20],
                b = NK_markers$HGNC_Symbol[21:35]),
    ID = "celltype:ch1", type = "NK",
    gene_id = "ENSEMBL"
  )
  expect_true(is.ggplot(p))
})
