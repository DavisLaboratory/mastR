#' @import ggplot2 ggpubr
NULL

tableau_20 <- c(
  "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
  "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
  "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
  "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
)

## plot functions for comparing unprocessed and processed data
## ---------------------------------------------------------------
# helper: calculate RLE
rle <- function(expr_mat, log = FALSE) {
  ## compute RLE
  if (log) {
    geomeans <- Matrix::rowMeans(log1p(expr_mat), na.rm = TRUE)
  } else {
    geomeans <- Matrix::rowMeans(expr_mat, na.rm = TRUE)
  }
  if (log) {
    dev <- log1p(expr_mat) - geomeans
  } else {
    dev <- expr_mat - geomeans
  }
  ## compute boxplot
  rledf <- t(apply(dev, 2, function(x) grDevices::boxplot.stats(x)$stats))
  rledf <- as.data.frame(rledf)
  colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")

  return(rledf)
}

# helper: Mean-Variance plot after voom
plot_voom <- function(vfit, span = 0.5) {
  stopifnot("vfit must be EList class!" = is(vfit, "EList"))

  ## Fit linear model to data
  fit <- limma::lmFit(vfit$E, design = vfit$design)
  ## calculate plot
  sx <- fit$Amean + mean(log2(vfit$targets$lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(vfit$E) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)

  ggplot(data.frame(x = sx, y = sy)) +
    geom_point(aes(x = x, y = y), size = 0.01) +
    geom_line(
      data = as.data.frame(l),
      aes(x = x, y = y),
      col = "red", linewidth = 0.5
    ) +
    theme_classic() +
    labs(
      x = "Average log-expression",
      y = "sqrt(sigma)",
      title = "voom:Mean-variance trend"
    )
}

plot_sa <- function(fit,
                    xlab = "Average log-expression",
                    ylab = "sqrt(sigma)",
                    title = "Final model: Mean-variance trend",
                    zero.weights = FALSE,
                    col = c("black", "red")) {
  if (!is(fit, "MArrayLM")) {
    stop("fit must be an MArrayLM object")
  }
  x <- fit$Amean
  y <- sqrt(fit$sigma)
  if (!(is.null(fit$weights) || zero.weights)) {
    allzero <- rowSums(fit$weights > 0, na.rm = TRUE) == 0
    y[allzero] <- NA
  }
  colv <- rep_len("Normal", nrow(fit))
  if (length(fit$df.prior) > 1L) {
    df2 <- max(fit$df.prior)
    s2 <- fit$sigma^2 / fit$s2.prior
    pdn <- pf(s2, df1 = fit$df.residual, df2 = df2)
    pup <- pf(s2, df1 = fit$df.residual, df2 = df2, lower.tail = FALSE)
    FDR <- p.adjust(2 * pmin(pdn, pup), method = "BH")
    colv[FDR <= 0.5] <- "Outlier"
  }

  p <- ggplot(data.frame(x = x, y = y, col = colv)) +
    geom_point(aes(x = x, y = y, col = col), size = 0.01) +
    scale_color_manual(values = c("Normal" = col[1], "Outlier" = col[2])) +
    theme_classic() +
    labs(
      x = xlab,
      y = ylab,
      title = title
    )

  if (!is.null(fit$s2.prior)) {
    if (length(fit$s2.prior) == 1L) {
      p <- p + geom_hline(yintercept = sqrt(sqrt(fit$s2.prior)), col = "blue")
    } else {
      o <- order(x)
      p <- p + geom_line(
        data = data.frame(
          x = x[o],
          y = sqrt(sqrt(fit$s2.prior[o]))
        ),
        aes(x = x, y = y),
        col = "blue", linewidth = 0.5
      )
    }
  }
  if (length(fit$df.prior) <= 1L) {
    p <- p + guides(col = "none")
  }

  return(p)
}


# helper: logCPM density plot between the 2 datasets
plot_density_init <- function(data1, data2, abl = 2,
                              expr = "logcounts", col = "Group",
                              sample = "Sample") {
  (ggplot(data1, aes(x = !!sym(expr), col = !!sym(col), group = !!sym(sample))) +
    geom_density() +
    geom_vline(xintercept = abl, lty = 2) +
    ggtitle("Original data") + ## original data1
    ggplot(data2, aes(x = logcounts, col = Group, group = Sample)) +
    geom_density() +
    geom_vline(xintercept = abl, lty = 2) +
    ggtitle("Processed data")) * ## processed data2
    theme_classic() *
    scale_color_manual(values = tableau_20) +
    patchwork::plot_layout(guides = "collect")
}

# helper: RLE plot between 2 datasets
plot_rle_init <- function(expr, group_col) {
  # constant - sample size at which standard plot becomes dense
  dense_thresh <- 50

  rledf <- rle(expr, log = FALSE)
  ## annotate samples
  rledf$MtchID <- rownames(rledf)
  rledf$Group <- group_col
  ## reorder samples
  rledf <- rledf[order(group_col), ]
  rledf$x <- seq_len(nrow(rledf))

  # build plot
  p1 <- ggplot2::ggplot(rledf, aes(
    x = x, y = middle, group = x,
    col = Group
  )) +
    ggplot2::geom_boxplot(
      aes(
        ymin = ymin, ymax = ymax, upper = upper,
        middle = middle, lower = lower
      ),
      stat = "identity"
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
    ggplot2::labs(x = "Samples", y = "Relative log expression") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = element_blank())

  # update plot if too many samples are plot
  if (nrow(rledf) > dense_thresh) {
    ## geom_point will inherit relevant aesthetics from top `aes`
    ## include y=middle
    p1 <- p1 + ggplot2::geom_point()
  }

  return(p1)
}

# helper: MDS plot between 2 datasets
plot_MDS_init <- function(expr1, expr2, group_col) {
  ## get MDS data
  mds1 <- limma::plotMDS(expr1, plot = FALSE)
  mds2 <- limma::plotMDS(expr2, plot = FALSE)

  mds1 <- data.frame(
    Sample = colnames(expr1),
    Group = group_col,
    "logFC_dim_1" = mds1$x,
    "logFC_dim_2" = mds1$y
  )
  mds2 <- data.frame(
    Sample = colnames(expr2),
    Group = group_col,
    "logFC_dim_1" = mds2$x,
    "logFC_dim_2" = mds2$y
  )
  p1 <- ggplot(mds1, aes(x = logFC_dim_1, y = logFC_dim_2, col = Group)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = tableau_20) +
    ggtitle("Original data") +
    theme_classic()
  p2 <- ggplot(mds2, aes(x = logFC_dim_1, y = logFC_dim_2, col = Group)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = tableau_20) +
    ggtitle("Processed data") +
    theme_classic()
  p1 + p2 + patchwork::plot_layout(guides = "collect")
}

## matrix plot function for PCA
## ---------------------------------------------------------------

## helper: Scale loadings coordinates to PCs
calculate_coor_end <- function(x, r, f = 1.5) {
  xend <- x * r * f
  return(xend)
}

## helper
#' @title Single PCA plot function
#'
#' @param prcomp prcomp object generated by [stats::prcomp()]
#' @param loading logical, if to plot and label loadings of PCA, default 'FALSE'
#' @param n_loadings num, top n loadings to plot; or a vector of gene IDs;
#'                   only work when `loading = TRUE`
#' @param dims a vector of 2 elements, specifying PCs to plot
#' @param group_by character, specify the column to be grouped and colored,
#'                 default NULL
#'
#' @return ggplot of PCA
plotPCAbiplot <- function(prcomp,
                          loading = FALSE,
                          n_loadings = 10,
                          dims = c(1, 2),
                          group_by = NULL) {
  stopifnot(
    "dims can only be numeric vector of 2 elements!" =
      is.numeric(dims) && length(dims) == 2
  )
  ## basic PCA plot
  p <- ggplot() +
    geom_point(aes(
      x = prcomp$x[, dims[1]],
      y = prcomp$x[, dims[2]],
      col = group_by
    )) +
    ggplot2::scale_color_manual(values = tableau_20) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = paste0(
        "PC", dims[1], "(",
        round(summary(prcomp)$importance[2, dims[1]] * 100, 2),
        "%)"
      ),
      y = paste0(
        "PC", dims[2], "(",
        round(summary(prcomp)$importance[2, dims[2]] * 100, 2),
        "%)"
      )
    )

  ## add loadings
  if (loading == TRUE) {
    ## extract top loadings
    if (is.numeric(n_loadings)) {
      pc1_genes <- -abs(prcomp$rotation[, dims[1]]) |>
        sort() |>
        head(n_loadings) |>
        names()
      pc2_genes <- -abs(prcomp$rotation[, dims[2]]) |>
        sort() |>
        head(n_loadings) |>
        names()
    } else if (is.character(n_loadings)) {
      stopifnot(
        "Not all n_loadings are in input prcomp data" =
          n_loadings %in% rownames(prcomp$rotation)
      )
      pc1_genes <- pc2_genes <- rownames(prcomp$rotation[n_loadings, ])
    }

    genes2plot <- unique(pc1_genes, pc2_genes)

    ## convert coordinates for loadings
    r <- min(
      (max(prcomp$x[, dims[1]]) - min(prcomp$x[, dims[1]]) /
        (max(prcomp$rotation[, dims[1]]) - min(prcomp$rotation[, dims[1]]))),
      (max(prcomp$x[, dims[2]]) - min(prcomp$x[, dims[2]]) /
        (max(prcomp$rotation[, dims[2]]) - min(prcomp$rotation[, dims[2]])))
    )

    loadings2plot <- as.data.frame(prcomp$rotation[genes2plot, dims])
    loadings2plot <- dplyr::mutate(
      loadings2plot,
      Gene = genes2plot,
      x_end = calculate_coor_end((!!sym(paste0("PC", dims[1]))), r),
      y_end = calculate_coor_end((!!sym(paste0("PC", dims[2]))), r)
    )

    ## add lodings
    p <- p +
      geom_segment(
        data = loadings2plot,
        aes(x = 0, y = 0, xend = x_end, yend = y_end),
        arrow = arrow(length = unit(1 / 2, "picas"), ends = "last"),
        color = "grey",
        linewidth = .5,
        alpha = 1,
        show.legend = NA
      ) +
      ggrepel::geom_text_repel(
        data = loadings2plot,
        aes(label = Gene, x = x_end, y = y_end, hjust = 0),
        color = "grey",
        size = 3
      )
  }
  return(p)
}

#' @title Make a matrix plot of PCA with top PCs
#'
#' @inheritParams plotPCAbiplot
#' @param data expression matrix
#' @param features vector of gene symbols or 'all', specify the genes used for
#'                 PCA, default 'all'
#' @param scale logical, if to scale data for PCA, default TRUE
#' @param n num, specify top n PCs to plot
#' @param gene_id character, specify which column of IDs used to calculate TPM,
#'                also indicate the ID type of expression data's rowname,
#'                could be one of 'ENSEMBL', 'SYMBOL', 'ENTREZ'...,
#'                default 'SYMBOL'
#'
#' @return matrix plot of PCA
pca_matrix_plot_init <- function(data,
                                 features = "all",
                                 group_by = NULL,
                                 scale = TRUE,
                                 n = 4,
                                 loading = FALSE,
                                 n_loadings = 10,
                                 gene_id = "SYMBOL") {
  stopifnot(
    is.vector(features),
    is.null(group_by) | is.vector(group_by),
    is.logical(scale), is.numeric(n),
    is.logical(loading), is.character(gene_id)
  )

  ## scale and do PCA on selected features
  if (length(features) == 1 && features[1] == "all") {
    tmp <- stats::prcomp(Matrix::t(data), scale = scale)
  } else {
    id <- na.omit(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      features, gene_id, "SYMBOL"
    ))
    id <- id[id %in% rownames(data)]
    tmp <- stats::prcomp(Matrix::t(data[id, ]), scale = scale)
    rownames(tmp$rotation) <- names(id)
  }

  ## screeplot of importance of PCs
  imp <- as.data.frame(t(summary(tmp)$importance))
  p0 <- ggplot2::ggplot(
    data = imp,
    aes(
      x = reorder(rownames(imp), seq_len(nrow(imp))),
      y = `Proportion of Variance`,
      group = 1
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "PCs") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  ## top n PCs matrix plot
  p1 <- lapply(seq_len(n), function(x) {
    ggplot2::ggplot() +
      ggplot2::geom_text(
        ggplot2::aes(
          x = 1, y = 1,
          label = paste0(
            "PC", x, "\n",
            round(summary(tmp)$importance[2, x] * 100, 2), "%"
          )
        ),
        size = 5
      ) +
      ggplot2::theme_void()
  })

  dims <- lapply(seq_len(n - 1), seq_len)
  names(dims) <- seq_len(n)[-1]
  dims <- apply(stack(dims), c(1, 2), as.numeric)
  dims <- dims[, 2:1]
  ## plot multi-PCs PCA
  p2 <- lapply(seq_len(nrow(dims)), \(i) plotPCAbiplot(
    prcomp = tmp,
    loading = loading,
    n_loadings = n_loadings,
    dims = dims[i, ],
    group_by = group_by
  ) + theme(legend.position = "none"))

  p3 <- ggpubr::as_ggplot(ggpubr::get_legend(
    ggplot() +
      geom_point(aes(x = 1, y = 1, col = group_by)) +
      scale_color_manual(values = tableau_20) +
      theme_classic()
  )) ## get legend

  layout_mat <- matrix(NA, n, n)
  layout_mat[upper.tri(layout_mat)] <- seq_along(p2)
  layout_mat[(ceiling(n / 2) + 1):n, seq_len(floor(n / 2))] <- length(p2) + 1
  diag(layout_mat) <- (length(p2) + 2):(length(p2) + 1 + n)

  p <- gridExtra::arrangeGrob(
    grobs = c(p2, list(p0), p1),
    layout_matrix = layout_mat
  )
  p <- ggpubr::as_ggplot(p)

  if (!is.null(group_by)) {
    p <- p + p3 + patchwork::plot_layout(widths = c(5, 1))
  }

  return(p)
}


## scatter_plot function
## ---------------------------------------------------------------
scatter_plot_init <- function(expr, sigs, target_group, by,
                              xint = 1, yint = 1,
                              gene_id = "SYMBOL") {
  stopifnot(is.numeric(xint), is.numeric(yint), is.character(gene_id))

  ## convert result gene symbols into gene IDs mapped to data rownames
  sigs <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    sigs,
    columns = gene_id,
    keytype = "SYMBOL"
  )
  sigs <- sigs[!duplicated(sigs[[gene_id]]), gene_id]
  sigs <- intersect(rownames(expr), sigs)

  ## calculate median expression of sigs by aggregating cells based on by group
  group <- ifelse(grepl(target_group, by), target_group, by)
  expr <- dplyr::group_by(
    as.data.frame(Matrix::t(expr[sigs, ])),
    factor(group)
  ) |>
    dplyr::summarise_all(.funs = mean, na.rm = TRUE) |>
    Matrix::t()
  colnames(expr) <- expr[1, ]
  expr <- as.data.frame(apply(expr[-1, ], MARGIN = c(1, 2), FUN = as.numeric))

  ## biplot of specified target_group vs all other types
  p <- tidyr::pivot_longer(expr, -target_group,
    names_to = "celltype",
    values_to = "Median_Expression"
  ) |>
    dplyr::mutate(
      flag = ifelse(!!sym(target_group) > yint & Median_Expression < xint,
        "High", "Low"
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x = Median_Expression,
      y = !!sym(target_group),
      col = flag
    )) +
    ggplot2::geom_point(alpha = 0.4, shape = 1) +
    ggplot2::geom_vline(xintercept = xint, lty = 2) +
    ggplot2::geom_hline(yintercept = yint, lty = 2) +
    ggplot2::facet_wrap(~celltype) +
    ggplot2::labs(x = "Median of z-scored Expression", col = "Specificity") +
    ggplot2::theme_classic()
  return(p)
}


## boxplot of expression for signatures
## ---------------------------------------------------------------
exp_boxplot_init <- function(expr, sigs, target_group, by,
                             method = "t.test", gene_id = "SYMBOL") {
  stopifnot(is.character(gene_id))

  ## convert result gene symbols into gene IDs mapped to data rownames
  sigs <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    sigs,
    columns = gene_id,
    keytype = "SYMBOL"
  )
  sigs <- sigs[!duplicated(sigs[[gene_id]]), gene_id]
  sigs <- intersect(rownames(expr), sigs)

  ## calculate median expression of sigs by aggregating samples on group
  group <- ifelse(grepl(target_group, by), target_group, as.character(by))
  expr <- dplyr::group_by(
    as.data.frame(Matrix::t(expr[sigs, ])),
    factor(group)
  ) |>
    dplyr::summarise_all(.funs = median, na.rm = TRUE) |>
    Matrix::t()
  colnames(expr) <- expr[1, ]
  expr <- as.data.frame(apply(expr[-1, ], MARGIN = c(1, 2), FUN = as.numeric))
  expr$Gene <- rownames(expr)

  ## boxplot of signature in specified target_group vs all other types
  p <- tidyr::pivot_longer(expr, -Gene,
    names_to = "Group",
    values_to = "Median Expression"
  ) |>
    ggplot2::ggplot(ggplot2::aes(
      x = !!sym("Group"),
      y = !!sym("Median Expression")
    )) +
    ggplot2::geom_boxplot(ggplot2::aes(col = Group)) +
    ggplot2::geom_point(ggplot2::aes(col = Group)) +
    ggplot2::geom_line(ggplot2::aes(group = Gene), col = "grey", alpha = .5) +
    ggpubr::stat_compare_means(
      label = "p.signif",
      method = method,
      ref.group = target_group
    ) +
    ggplot2::labs(
      title = "Signatures Median Expression across Groups",
      x = "Groups", y = "Expression"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  return(p)
}


## score boexplot
## ---------------------------------------------------------------
# helper: calculate singscores
singscore_init <- function(expr, sigs, by,
                           gene_id = "SYMBOL") {
  ## singscore samples using given signature genes
  rank_data <- singscore::rankGenes(expr)
  UP <- na.omit(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    sigs,
    column = gene_id,
    keytype = "SYMBOL"
  ))

  scores <- singscore::simpleScore(rank_data, upSet = UP)
  return(scores)
}

# helper: boxplot of comparing scores
score_boxplot_init <- function(scores, by, target_group, method) {
  ## set Group for samples by given factor
  scores$Group <- ifelse(grepl(target_group, by),
    target_group,
    as.character(by)
  )

  ## plot
  p <- ggplot2::ggplot(
    data = scores,
    ggplot2::aes(x = Group, y = TotalScore)
  ) +
    ggplot2::geom_boxplot(ggplot2::aes(col = Group)) +
    ggplot2::geom_jitter(ggplot2::aes(col = Group), alpha = 0.2) +
    ggpubr::stat_compare_means(
      label = "p.signif", method = method,
      ref.group = target_group, label.y.npc = 1
    ) +
    ggplot2::labs(title = "Scores across Groups") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  return(p)
}


## GSEA plot for signatures
## ---------------------------------------------------------------
## gse analysis
gsea_analysis <- function(tDEG, gsets, gene_id = "SYMBOL", digits = 2) {
  gse <- lapply(names(tDEG), function(n) {
    glist <- tDEG[[n]][, "logFC"]
    names(glist) <- rownames(tDEG[[n]])
    glist <- sort(glist, decreasing = TRUE)

    ## pick common genes in glist and gsets
    gsets <- gsets[gsets[, gene_id] %in% names(glist), ]

    gse <- clusterProfiler::GSEA(glist, TERM2GENE = gsets[, c("set", gene_id)])
    if (nrow(gse) == 0) {
      ms <- paste("No term was found enriched in", n, ".")
      message(ms)
      invisible()
    } else {
      ## only keep 2 digits for pvalue table
      slot(gse, "result")$p.adjust <- signif(slot(gse, "result")$p.adjust, digits = digits)
      slot(gse, "result")$pvalue <- signif(slot(gse, "result")$pvalue, digits = digits)
      ## save comparison info
      slot(gse, "result")$group <- n
      return(gse)
    }
  })
  names(gse) <- names(tDEG)
  return(gse)
}
## gseaplot
gsea_plot_init <- function(gse, pvalue_table = FALSE) {
  p <- lapply(names(gse), function(n) {
    if (is.null(gse[[n]])) {
      ms <- paste("No term was found enriched in", n, ".")
      message(ms)
      p <- ggpubr::as_ggplot(grid::textGrob(ms))
    } else {
      p <- enrichplot::gseaplot2(gse[[n]],
        geneSetID = seq_len(nrow(gse[[n]])),
        pvalue_table = pvalue_table,
        title = n
      )
      p <- patchwork::wrap_elements(patchwork::wrap_plots(p, ncol = 1))
      # ## add pvalue_table
      # p <- p - patchwork::inset_element(gridExtra::tableGrob(
      #   gse[[n]][,c("pvalue", "p.adjust")] |> signif(digits),
      #   theme = gridExtra::ttheme_default(base_size = 5)
      # ),
      # left = 0.6, bottom = 0.8, right = 1, top = 1)
    }
  })
  names(p) <- names(gse)
  return(p)
}
## dotplot
gsea_dotplot_init <- function(gse,
                              col = "-log10(p.adjust)",
                              size = "enrichmentScore") {
  res <- do.call(rbind, lapply(gse, \(x) {
    if (is.null(x)) {
      invisible()
    } else {
      return(slot(x, "result"))
    }
  }))

  if (is.null(res)) {
    ms <- "No term was found enriched in any comparison!"
    message(ms)
    p <- ggpubr::as_ggplot(grid::textGrob(ms))
  } else {
    res[["-log10(p.adjust)"]] <- -log10(res$p.adjust)
    p <- ggplot(res) +
      geom_point(aes(
        x = ID, y = group,
        col = !!sym(col),
        size = !!sym(size)
      )) +
      labs(x = "Comparison", y = "Signature") +
      theme_bw() +
      scale_color_viridis_c() +
      theme(axis.text.x = element_text(angle = 90))
  }

  return(p)
}


## heatmap for comparing signature with original markers pool
## ---------------------------------------------------------------
heatmap_init <- function(expr, sigs, by, markers,
                         scale = c("none", "row", "column"),
                         gene_id = "SYMBOL", ranks_plot = FALSE,
                         ...) {
  scale <- match.arg(scale)
  ## name sigs list
  if (is.null(names(sigs))) {
    names(sigs) <- seq_along(sigs)
  }
  idx <- which(names(sigs) == "")
  names(sigs)[idx] <- idx

  ## plot rank instead of expression
  if (ranks_plot == TRUE) {
    expr <- log2(apply(expr, 2, rank))
  }
  ## scale
  if (scale == "row") {
    expr <- t(scale(t(expr)))
  } else if (scale == "column") {
    expr <- scale(expr)
  }
  ## convert sigs
  sigs <- utils::stack(sigs)
  sigs[[gene_id]] <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    sigs$values, gene_id, "SYMBOL"
  )
  ## kepp common genes
  idx <- sigs[[gene_id]] %in% rownames(expr)
  if (any(idx == FALSE)) {
    message(paste("Gene", sigs$values[!idx], "is not in data\n"))
  }
  sigs <- sigs[idx, ]

  ## subset expr
  expr <- expr[sigs[[gene_id]], ]
  rownames(expr) <- sigs$values

  ## annotation
  top_ann <- ComplexHeatmap::columnAnnotation(Group = by)
  right_ann <- ComplexHeatmap::rowAnnotation(Sigs = sigs$ind)

  ## plot
  p <- ComplexHeatmap::Heatmap(
    expr,
    column_split = by,
    row_split = sigs$ind,
    top_annotation = top_ann,
    right_annotation = right_ann,
    ...
  )

  return(p)
}


## rankdensity plot of signature
## ---------------------------------------------------------------
# helper: aggregate samples of the same group to one sample with mean expression
agg_mean <- function(expr, by) {
  expr <- vapply(split(colnames(expr), by),
    FUN = function(x) Matrix::rowMeans(expr[, x], na.rm = TRUE),
    FUN.VALUE = rep(1, nrow(expr))
  )
  return(expr)
}

# helper: arrange plot matrix
# add plot_spacer() at the end for samples less than max
add_spacer <- function(x, ns) {
  if (length(x) < ns) {
    spacer <- paste(
      c(
        "list(",
        paste(
          rep(
            "patchwork::plot_spacer()",
            ns - length(x)
          ),
          collapse = ","
        ),
        ")"
      ),
      collapse = ""
    )
    x <- append(x, eval(parse(text = spacer)), length(x))
  }
  x
}

# helper: make matrix plot of rankdensity
rankdensity_init <- function(expr, sigs, by,
                             aggregate = FALSE,
                             gene_id = "SYMBOL") {
  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }

  ## calculate ranks of genes using singscore
  rank_data <- singscore::rankGenes(expr)
  UP <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    sigs,
    columns = gene_id,
    keytype = "SYMBOL"
  )
  UP <- intersect(UP[!duplicated(UP[[gene_id]]), gene_id], rownames(expr))

  data <- data.frame(t(rank_data[UP, ] / nrow(rank_data)), Group = by)
  data$Sample <- rownames(data)

  if (aggregate == TRUE) {
    # update y-positions for barcode of up-set
    data <- data |>
      tidyr::pivot_longer(-c("Sample", "Group"),
        names_to = "Gene",
        values_to = "Rank"
      ) |>
      dplyr::group_by(Group) |>
      dplyr::mutate(max_den = max(density(Rank)$y))

    p <- ggplot(data, aes(x = Rank, col = Group, group = Group)) +
      geom_density() +
      geom_segment(aes(y = max_den, xend = Rank, yend = max_den * 1.1)) +
      facet_wrap(~Group, scales = "free") +
      theme_bw()
  } else {
    # update y-positions for barcode of up-set
    data <- data |>
      tidyr::pivot_longer(-c("Sample", "Group"),
        names_to = "Gene",
        values_to = "Rank"
      ) |>
      dplyr::group_by(Group) |>
      dplyr::mutate(y_s = as.numeric(factor(Sample))) |>
      dplyr::mutate(y_s = -y_s / max(y_s), y_e = y_s - 0.1)

    p <- ggplot(data, aes(x = Rank, col = Group, group = Sample)) +
      geom_density() +
      geom_segment(aes(y = y_s, xend = Rank, yend = y_e)) +
      facet_wrap(~Group, scales = "free") +
      theme_bw()
  }

  return(p)
}



## scatter plot of HDBSCAN cluster result
## ---------------------------------------------------------------
scatter_hdb_cl <- function(sig_matrix, markers_list) {
  stopifnot(
    "sig_matrix must be a matrix or data.frame!" = is.matrix(sig_matrix) | is.data.frame(sig_matrix),
    "markers_list is not a list!" = is.list(markers_list)
  )

  ml <- utils::stack(markers_list)
  ml <- split(ml$ind, ml$values)

  p <- lapply(rownames(sig_matrix), function(m) {
    p <- ggplot2::ggplot(
      data.frame(
        Gene = m,
        Expression = as.numeric(sig_matrix[m, ]),
        Type = ifelse(colnames(sig_matrix) %in% ml[[m]],
          colnames(sig_matrix),
          ""
        )
      ),
      ggplot2::aes(
        x = Gene, y = Expression,
        col = (Type != "")
      )
    ) +
      ggrepel::geom_label_repel(ggplot2::aes(label = Type)) +
      ggplot2::geom_point() +
      ggplot2::labs(col = "is.marker") +
      ggplot2::facet_wrap(~Gene) +
      ggplot2::scale_color_manual(values = c(
        "TRUE" = "red",
        "FALSE" = "black"
      )) +
      ggplot2::theme_classic()
  })

  names(p) <- rownames(sig_matrix)

  return(p)
}

utils::globalVariables(c(
  "median", "Gene", "Expression", "Group", "Type",
  "target_group", "group_col", "TotalScore",
  "Proportion of Variance", "Rank", "Sample",
  "max_den", "y_s", "y_e", "x_end", "y_end",
  "group", "ID", "NES", "logFC_dim_1",
  "logFC_dim_2", "logcounts", "x", "y", "middle",
  "ymin", "ymax", "upper", "lower", "flag", "Median_Expression"
))
