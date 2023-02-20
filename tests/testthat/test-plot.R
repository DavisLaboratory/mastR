# library(vdiffr)

test_that("plot_density works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  p <- function() plot_density(dge = dge, ID = "group",
                               keep = 1:1e4, counts = TRUE)

  # expect_doppelganger("basic density plot", p)
  expect_silent(p())
})

test_that("plot_rle works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  p <- function() plot_rle(dge = dge, ID = "group", keep = 1:1e4, counts = TRUE)

  expect_true(is(rle(dge$counts), "matrix"))
  # expect_doppelganger("basic RLE boxplot", p)
  expect_silent(p())
})

test_that("plot_MDS works", {
  data("im_data_6")
  dge <- edgeR::DGEList(counts = Biobase::exprs(im_data_6),
                        group = im_data_6$`celltype:ch1`)

  p <- function() plot_MDS(dge = dge, ID = "group", keep = 1:1e4, counts = TRUE)

  # expect_doppelganger("basic MDS plot", p)
  expect_silent(p())
})

test_that("scatter_plot_init works", {

  data("im_data_6", "NK_markers")
  p <- scatter_plot_init(expr = Biobase::exprs(im_data_6),
                         sigs = NK_markers$HGNC_Symbol,
                         type = "NK",
                         by = im_data_6$`celltype:ch1`,
                         gene_id = "ENSEMBL")

  expect_true(is.ggplot(p))
})
