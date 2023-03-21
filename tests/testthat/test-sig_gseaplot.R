test_that("sig_gseaplot works", {
  data("im_data_6", "NK_markers")
  ## test matrix object
  suppressWarnings(
    p <- sig_gseaplot(
      im_data_6@assayData$exprs,
      sigs = NK_markers$HGNC_Symbol[1:20],
      group_col = im_data_6$`celltype:ch1`,
      target_group = "NK",
      gene_id = "ENSEMBL"
    )
  )
  expect_true(is.ggplot(p))

  ## test eSet object
  ### test expression boxplot
  suppressWarnings(
    p <- sig_gseaplot(
      im_data_6,
      sigs = NK_markers$HGNC_Symbol[1:20],
      group_col = "celltype:ch1", target_group = "NK",
      gene_id = "ENSEMBL"
    )
  )
  expect_true(is.ggplot(p))

  ## test DGEList object
  dge <- edgeR::DGEList(
    counts = im_data_6@assayData$exprs,
    group = im_data_6$`celltype:ch1`
  )
  suppressWarnings(
    p <- sig_gseaplot(
      dge,
      sigs = NK_markers$HGNC_Symbol[1:20],
      group_col = "group", target_group = "NK",
      gene_id = "ENSEMBL"
    )
  )
  expect_true(is.ggplot(p))

  ## test seurat object
  data_seurat <- Seurat::CreateSeuratObject(
    counts = im_data_6@assayData$exprs,
    meta.data = dge$samples
  )
  suppressWarnings(
    p <- sig_gseaplot(
      data_seurat,
      sigs = NK_markers$HGNC_Symbol[1:20],
      group_col = "group", target_group = "NK",
      gene_id = "ENSEMBL"
    )
  )
  expect_true(is.ggplot(p))

  ## test list objects with list of signatures
  suppressWarnings(
    p <- sig_gseaplot(
      list(A = im_data_6, B = im_data_6),
      sigs = list(
        a = NK_markers$HGNC_Symbol[1:20],
        b = NK_markers$HGNC_Symbol[21:35]
      ),
      group_col = "celltype:ch1", target_group = "NK",
      gene_id = "ENSEMBL"
    )
  )
  expect_true(is.ggplot(p))

  ## test gseaplot2
  suppressWarnings(
    p <- sig_gseaplot(
      im_data_6,
      sigs = NK_markers$HGNC_Symbol[1:20],
      group_col = "celltype:ch1", target_group = "NK",
      gene_id = "ENSEMBL", method = "gseaplot"
    )
  )
  expect_true(is.ggplot(p))
})
