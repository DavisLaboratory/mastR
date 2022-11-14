test_that("pseudo_samples works", {
  # test matrix
  counts <- matrix(abs(rnorm(10000, 10, 10)), 100)
  rownames(counts) <- 1:100
  colnames(counts) <- 1:100
  meta <- data.frame(subset = rep(c("A", "B"), 50),
                     level = rep(1:4, each = 25))
  rownames(meta) <- 1:100
  pb <- pseudo_samples(counts, by = meta,
                       min.cells = 10, max.cells = 20)

  expect_true(ncol(pb) == 8)

  # test seurat
  ## generate scRNA data
  scRNA <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  pb <- pseudo_samples(scRNA, by = c("subset", "level"),
                       min.cells = 10, max.cells = 20)

  expect_true(ncol(pb) == 8)

  # test sce
  pb <- pseudo_samples(scRNA |> Seurat::as.SingleCellExperiment(),
                       by = c("subset", "level"),
                       min.cells = 10, max.cells = 20)

  expect_true(ncol(pb) == 8)
})
