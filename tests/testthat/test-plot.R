# library(vdiffr)

test_that("plot_density_init works", {
  data("im_data_6")
  data <- data.frame(
    logcounts = rnorm(100),
    Group = rep(c("a", "b"), 50),
    Sample = rep(1:4, each = 25)
  )

  p <- plot_density_init(data1 = data, data2 = data[1:80, ])

  # expect_doppelganger("basic density plot", p)
  expect_true(is.ggplot(p))
})

test_that("plot_rle_init works", {
  data("im_data_6")

  p <- function() {
    plot_rle_init(
      expr = Biobase::exprs(im_data_6)[1:1000, ],
      group_col = im_data_6$`celltype:ch1`
    )
  }
  expect_silent(p())
})

test_that("plot_MDS_init works", {
  data("im_data_6")

  p <- function() {
    plot_MDS_init(
      expr1 = Biobase::exprs(im_data_6)[1:1000, ],
      expr2 = Biobase::exprs(im_data_6)[1001:2000, ],
      group_col = im_data_6$`celltype:ch1`
    )
  }

  # expect_doppelganger("basic MDS plot", p)
  expect_silent(p())
})

test_that("scatter_plot_init works", {
  data("im_data_6", "NK_markers")
  p <- scatter_plot_init(
    expr = Biobase::exprs(im_data_6),
    sigs = NK_markers$HGNC_Symbol,
    target_group = "NK",
    by = im_data_6$`celltype:ch1`,
    gene_id = "ENSEMBL"
  )

  expect_true(is.ggplot(p))
})
