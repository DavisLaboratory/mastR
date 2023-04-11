#' @title Aggregate single cells to pseudo-samples according to specific factors
#'
#' @description Gather cells for each group according to specified factors,
#'   then randomly assign and aggregate cells to each pseudo-samples with randomized
#'   cell size. (min.cells <= size <= max.cells)
#'
#' @param data a matrix or Seurat/SCE object containing expression and metadata
#' @param by a vector of group names or dataframe for aggregation
#' @param fun chr, methods used to aggregate cells, could be 'sum' or 'mean',
#'            default 'sum'
#' @param scale a num or NULL, if to multiply a scale to the average expression
#' @param min.cells num, default 300, the minimum size of cells aggregating
#'                  to each pseudo-sample
#' @param max.cells num, default 600, the maximum size of cells aggregating
#'                  to each pseudo-sample
#' @param slot chr, specify which slot of seurat object to aggregate, can be
#'             'counts', 'data', 'scale.data'..., default is 'counts'
#'
#' @return An expression matrix after aggregating cells on specified factors
#'
#' @examples
#' counts <- matrix(abs(rnorm(10000, 10, 10)), 100)
#' rownames(counts) <- 1:100
#' colnames(counts) <- 1:100
#' meta <- data.frame(
#'   subset = rep(c("A", "B"), 50),
#'   level = rep(1:4, each = 25)
#' )
#' rownames(meta) <- 1:100
#' scRNA <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
#' pseudo_samples(scRNA,
#'   by = c("subset", "level"),
#'   min.cells = 10, max.cells = 20
#' )
#'
#' @export
setGeneric(
  "pseudo_samples",
  function(data,
           by,
           fun = c("sum", "mean"),
           scale = NULL,
           min.cells = 0,
           max.cells = Inf,
           slot = "counts") {
    standardGeneric("pseudo_samples")
  }
)

#' @rdname pseudo_samples
setMethod(
  "pseudo_samples", signature(
    data = "matrix",
    by = "data.frame"
  ),
  function(data,
           by,
           fun = c("sum", "mean"),
           scale = NULL,
           min.cells = 0,
           max.cells = Inf,
           slot = "counts") {
    fun <- match.arg(fun)
    stopifnot(
      is.character(fun), is.numeric(min.cells),
      is.numeric(max.cells), is.character(slot)
    )

    l <- pseudo_sample_list(
      data = data, by = by,
      min.cells = min.cells,
      max.cells = max.cells
    )
    p_samples <- unlist(lapply(names(l), function(x) {
      paste(x, names(l[[x]]), sep = "_")
    }))
    p_samples <- p_samples[!grepl("_$", p_samples)]
    pb <- lapply(l, function(i) {
      vapply(i, function(j) {
        if (fun == "mean") {
          Matrix::rowMeans(data[, j], na.rm = TRUE)
        } else if (fun == "sum") Matrix::rowSums(data[, j], na.rm = TRUE)
      }, FUN.VALUE = rep(1, nrow(data)))
    })
    pb <- do.call(cbind, pb)
    pb <- pb[, which(!is.na(colnames(pb)))]
    colnames(pb) <- p_samples
    pb <- apply(pb, c(1, 2), as.numeric)
    if (fun == "mean" && !is.null(scale)) {
      pb <- pb * scale
    }
    return(pb)
  }
)

#' @rdname pseudo_samples
setMethod(
  "pseudo_samples", signature(
    data = "matrix",
    by = "vector"
  ),
  function(data,
           by,
           fun = c("sum", "mean"),
           scale = NULL,
           min.cells = 0,
           max.cells = Inf,
           slot = "counts") {
    stopifnot(
      is.character(fun), is.numeric(min.cells),
      is.numeric(max.cells), is.character(slot)
    )

    pb <- pseudo_samples(
      data = data, by = as.data.frame(by), fun = fun,
      scale = scale, min.cells = min.cells,
      max.cells = max.cells, slot = slot
    )
    return(pb)
  }
)

#' @rdname pseudo_samples
setMethod(
  "pseudo_samples", signature(
    data = "Seurat",
    by = "character"
  ),
  function(data,
           by,
           fun = c("sum", "mean"),
           scale = NULL,
           min.cells = 0,
           max.cells = Inf,
           slot = "counts") {
    stopifnot(
      is.character(fun), is.numeric(min.cells),
      is.numeric(max.cells), is.character(slot)
    )

    expr <- as.matrix(SeuratObject::GetAssayData(data, slot = slot))
    coldata <- slot(data, "meta.data")[, by]
    rm(data)

    pb <- pseudo_samples(
      data = expr, by = coldata, fun = fun,
      scale = scale, min.cells = min.cells,
      max.cells = max.cells, slot = slot
    )
    return(pb)
  }
)

#' @rdname pseudo_samples
setMethod(
  "pseudo_samples", signature(
    data = "SummarizedExperiment",
    by = "character"
  ),
  function(data,
           by,
           fun = c("sum", "mean"),
           scale = NULL,
           min.cells = 0,
           max.cells = Inf,
           slot = "counts") {
    stopifnot(
      is.character(fun), is.numeric(min.cells),
      is.numeric(max.cells), is.character(slot)
    )

    expr <- as.matrix(SummarizedExperiment::assay(data, slot))
    coldata <- as.data.frame(SummarizedExperiment::colData(data)[, by])
    rm(data)

    pb <- pseudo_samples(
      data = expr, by = coldata, fun = fun,
      scale = scale, min.cells = min.cells,
      max.cells = max.cells, slot = slot
    )
    return(pb)
  }
)

#' Split cells according to specific factors
#'
#' Gathering cells to make the pool according to specific factors, and randomly
#' assign the cells from the pool to pseudo-sample with the randomized cell
#' size. (min.cells <= size <= max.cells)
#'
#' @param data matrix or data.frame or other single cell expression object
#' @param by a vector or data.frame contains factor(s) for aggregation
#' @param min.cells num, default 0, the minimum size of cells aggregating
#'                  to each pseudo-sample
#' @param max.cells num, default Inf, the maximum size of cells aggregating
#'                  to each pseudo-sample
#'
#' @return A list of cell names for each pseudo-sample
#' @export
#'
#' @examples
#' counts <- matrix(abs(rnorm(10000, 10, 10)), 100)
#' rownames(counts) <- 1:100
#' colnames(counts) <- 1:100
#' meta <- data.frame(
#'   subset = rep(c("A", "B"), 50),
#'   level = rep(1:4, each = 25)
#' )
#' rownames(meta) <- 1:100
#' scRNA <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
#' pseudo_sample_list(scRNA,
#'   by = c("subset", "level"),
#'   min.cells = 10, max.cells = 20
#' )
pseudo_sample_list <- function(data,
                               by,
                               min.cells = 0,
                               max.cells = Inf) {
  ## split cell groups based on given factor(s)
  l <- split(colnames(data), by)

  psamples <- lapply(l, function(cells) {
    i <- 1
    psample <- list()
    ## spread cells into multiple pseudo samples with cell counts from min.cells
    ## to max.cells when max.cells != Inf
    while (length(cells) > max.cells) {
      psample[[i]] <- sample(cells, runif(1, min.cells, max.cells), replace = FALSE)
      cells <- setdiff(cells, psample[[i]])
      i <- i + 1
    }
    ## only keep pseudo samples with cell size > min.cells
    if (length(cells) >= min.cells) psample[[i]] <- cells
    names(psample) <- seq_along(psample)
    return(psample)
  })
  return(psamples)
}
