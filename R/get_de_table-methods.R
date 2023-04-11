#' @include DE_functions.R plot.R
NULL

#' @rdname get_de_table
setMethod(
  "get_de_table", signature(
    data = "DGEList",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           slot = "counts",
           ...) {
    stopifnot("slot must be character!" = is.character(slot))

    DGE <- data
    DGE$samples$group <- DGE$samples[[group_col]]
    DGE$original <- DGE$counts
    DGE$counts <- DGE[[slot]]
    rm(data)

    # ## standard DE analysis with edgeR and limma::voom pipeline
    # voom_res <- de_analysis(
    #   dge = DGE,
    #   group_col = group_col,
    #   target_group = target_group,
    #   ...
    # )
    proc_data <- process_data(
      data = DGE,
      group_col = group_col,
      target_group = target_group,
      ...
    )

    ## save DE result tables into list
    DE_table <- list()
    for (i in seq_len(ncol(proc_data$tfit))) {
      ## use limma::topTreat() to get statistics of DEA
      DE_table[[i]] <- na.omit(limma::topTreat(proc_data$tfit,
        coef = i,
        number = Inf
      ))
    }
    names(DE_table) <- colnames(proc_data$tfit)

    return(DE_table)
  }
)

#' @rdname get_de_table
setMethod(
  "get_de_table", signature(
    data = "matrix",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           ...) {
    DGE <- edgeR::DGEList(counts = data, group = group_col)
    group_col <- "group"
    rm(data)

    DE_table <- get_de_table(
      data = DGE, group_col = group_col,
      target_group = target_group,
      slot = "counts", ...
    )

    return(DE_table)
  }
)

#' @rdname get_de_table
setMethod(
  "get_de_table", signature(
    data = "Matrix",
    group_col = "vector",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           ...) {
    DGE <- edgeR::DGEList(counts = data, group = group_col)
    group_col <- "group"
    rm(data)

    DE_table <- get_de_table(
      data = DGE, group_col = group_col,
      target_group = target_group,
      slot = "counts", ...
    )

    return(DE_table)
  }
)

#' @rdname get_de_table
setMethod(
  "get_de_table", signature(
    data = "ExpressionSet",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           ...) {
    expr <- Biobase::exprs(data)
    coldata <- Biobase::pData(data)

    DGE <- edgeR::DGEList(
      counts = expr,
      samples = coldata,
      group = coldata[[group_col]]
    )
    group_col <- make.names(group_col)
    rm(data, expr, coldata)

    DE_table <- get_de_table(
      data = DGE, group_col = group_col,
      target_group = target_group,
      slot = "counts", ...
    )

    return(DE_table)
  }
)

#' @rdname get_de_table
setMethod(
  "get_de_table", signature(
    data = "SummarizedExperiment",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           slot = "counts",
           ...) {
    expr <- SummarizedExperiment::assay(data, slot)
    coldata <- SummarizedExperiment::colData(data)

    DGE <- edgeR::DGEList(
      counts = expr,
      samples = coldata,
      group = coldata[[group_col]]
    )
    group_col <- make.names(group_col)
    rm(data, expr, coldata)

    DE_table <- get_de_table(
      data = DGE, group_col = group_col,
      target_group = target_group,
      slot = "counts", ...
    )

    return(DE_table)
  }
)

#' @rdname get_de_table
setMethod(
  "get_de_table", signature(
    data = "Seurat",
    group_col = "character",
    target_group = "character"
  ),
  function(data,
           group_col,
           target_group,
           slot = "counts",
           ...) {
    expr <- SeuratObject::GetAssayData(data, slot = slot)
    coldata <- data@meta.data

    DGE <- edgeR::DGEList(
      counts = expr,
      samples = coldata,
      group = coldata[[group_col]]
    )
    group_col <- make.names(group_col)
    rm(data, expr, coldata)

    DE_table <- get_de_table(
      data = DGE, group_col = group_col,
      target_group = target_group,
      slot = "counts", ...
    )

    return(DE_table)
  }
)
