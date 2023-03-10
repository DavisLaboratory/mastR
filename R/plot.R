#' @import ggplot2 ggpubr
NULL

tableau_20 <-  c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
                 "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
                 "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
                 "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")

## plot functions for comparing unfiltered and filtered data
##---------------------------------------------------------------
#helper: logCPM density plot between filtered and unfiltered
plot_density <- function(dge, group_col, keep = TRUE,
                         normalize = TRUE, filter = 10) {
  ## set colors for each cell type
  col <- dge$samples[[group_col]] |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = tableau_20)(unique(dge$samples[[group_col]]) |> length())

  ## set par for plot to make 2 plots in one page
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))

  ## calculate logCPM
  if(normalize) {
    M <- median(dge$samples$lib.size) * 1e-6
    L <- mean(dge$samples$lib.size) * 1e-6
    lcpm.1 <- edgeR::cpm(
      edgeR::calcNormFactors(dge, method = "TMM"),
      log = TRUE
    )
    lcpm.2 <- edgeR::cpm(
      edgeR::calcNormFactors(dge[keep,, keep.lib.sizes = FALSE]),
      log = TRUE
    )
    abl <- log2(filter/M + 2/L)
  }else {
    # lcpm.1 <- log1p(dge$counts)
    # lcpm.2 <- log1p(dge$counts[keep,])
    # abl <- log1p(filter[1])
    lcpm.1 <- dge$counts
    lcpm.2 <- dge$counts[keep,]
    abl <- filter[1]
  }

  ## plot unfiltered data
  plot(stats::density(lcpm.1[,1]), col = colors[col[1]],
       main = "Unfiltered", xlab = ifelse(normalize, "logCPM", "logcounts"),
       ylim = c(0, max(stats::density(lcpm.1)$y) * 1.2))
  abline(v = abl, lty = 3)
  lapply(2:ncol(dge), function(i) lines(stats::density(lcpm.1[,i]),
                                        col = colors[col[i]]))
  legend("topright", legend = unique(dge$samples[[group_col]]),
         text.col = colors[unique(col)], bty = "n", inset = c(-0.3, 0))

  ## plot filtered data
  plot(stats::density(lcpm.2[,1]), col = colors[col[1]],
       main = "Filtered", xlab = ifelse(normalize, "logCPM", "logcounts"),
       ylim = c(0, max(stats::density(lcpm.2)$y) * 1.2))
  abline(v = abl, lty = 3)
  lapply(2:ncol(dge), function(i) lines(stats::density(lcpm.2[,i]),
                                        col = colors[col[i]]))
  legend("topright", legend = unique(dge$samples[[group_col]]),
         text.col = colors[unique(col)], bty = "n", inset = c(-0.3, 0))

  mtext("Input Dataset", side = 3, line = 0, outer = TRUE)
  ## set par back to original setting
  par(op)
}

#helper: calculate RLE
rle <- function(expr_mat, log = FALSE) {
  if (log) {
    geomeans <- Matrix::rowMeans(log1p(expr_mat), na.rm = TRUE)
  } else {
    geomeans <- Matrix::rowMeans(expr_mat, na.rm = TRUE)
  }
  RLE <- function(cnts) {
    return(cnts - geomeans)
  }
  if (log) {
    dev <- apply(log1p(expr_mat), 2, RLE)
  } else {
    dev <- apply(expr_mat, 2, RLE)
  }
  return(dev)
}

#helper: RLE plot between before and after filtration
plot_rle <- function(dge, group_col, keep = TRUE, normalize = TRUE) {
  ## set colors for each cell type
  col <- dge$samples[[group_col]] |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = tableau_20)(unique(dge$samples[[group_col]]) |>
                                                    length())

  ## calculate logCPM
  if(normalize) {
    lcpm.1 <- edgeR::cpm(
      edgeR::calcNormFactors(dge, method = "TMM"),
      log = TRUE
    )
    lcpm.2 <- edgeR::cpm(
      edgeR::calcNormFactors(dge[keep,, keep.lib.sizes = FALSE]),
      log = TRUE
    )
  }else {
    # lcpm.1 <- log1p(dge$counts)
    # lcpm.2 <- log1p(dge$counts[keep,])
    lcpm.1 <- dge$counts
    lcpm.2 <- dge$counts[keep, ]
  }

  ## set par for boxplot and legend
  op <- par(no.readonly = TRUE)
  par(mfrow = c(3, 1), mar = c(1, 2.5, 2, 1), oma = c(0, 0, 1, 0))

  ## plot unfiltered data
  rle(lcpm.1[, order(col)]) |>
    boxplot(
      outline = FALSE, xaxt = "n", col = sort(col) |> (\(x) colors[x])(),
      main = "RLE plot before filtration"
    )
  abline(h = 0)
  ##plot filtered data
  rle(lcpm.2[, order(col)]) |>
    boxplot(
      outline = FALSE, xaxt = "n", col = sort(col) |> (\(x) colors[x])(),
      main = "RLE plot after filtration"
    )
  abline(h = 0)
  ## add legend
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = unique(dge$samples[[group_col]]) |> sort(),
         fill = colors, xpd = TRUE, ncol = 4
  )
  ## set par back to original
  par(op)
}

#helper: MDS plot before and after filtration
plot_MDS <- function(dge, group_col, keep = TRUE, normalize = TRUE) {
  ## set colors for each cell type
  col <- dge$samples[[group_col]] |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = tableau_20)(unique(dge$samples[[group_col]]) |>
                                                    length())

  ## calculate logCPM
  if(normalize) {
    lcpm.1 <- edgeR::cpm(
      edgeR::calcNormFactors(dge, method = "TMM"),
      log = TRUE
    )
    lcpm.2 <- edgeR::cpm(
      edgeR::calcNormFactors(dge[keep,, keep.lib.sizes = FALSE]),
      log = TRUE
    )
  }else {
    # lcpm.1 <- log1p(dge$counts)
    # lcpm.2 <- log1p(dge$counts[keep,])
    lcpm.1 <- dge$counts
    lcpm.2 <- dge$counts[keep, ]
  }

  ## set par to make 2 plots in one page
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 3), oma = c(0, 0, 3, 0))
  limma::plotMDS(lcpm.1,
                 pch = 1,
                 col = colors[col],
                 main = "MDS plot before filtration"
  )
  limma::plotMDS(lcpm.2,
                 pch = 1,
                 col = colors[col],
                 main = "MDS plot after filtrattion"
  )
  ## add legend
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = unique(dge$samples[[group_col]]) |> sort(),
         fill = colors, xpd = TRUE
  )
  mtext("Input Data", side = 3, line = 0, outer = TRUE)
  ## set par back to original
  par(op)
  return(NULL)
}

#helper: Mean-Variance plot after voom
plot_voom <- function(vfit, span = 0.5) {

  stopifnot("vfit must be EList class!" = is(vfit, "EList"))

  ## Fit linear model to data
  fit <- limma::lmFit(vfit$E, design = vfit$design)
  ## Fit lowess trend to sqrt-sigma by log-expression
  l <- lowess(fit$Amean, sqrt(fit$sigma), f = span)

  limma::plotSA(fit)
  title("voom: Mean-variance trend")
  lines(l, col="red")
}


#helper: density plot init
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

#helper: rle plot init
plot_rle_init <- function(expr1, expr2, group_col) {

  ## set colors for each cell type
  col <- group_col |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = tableau_20)(length(unique(group_col)))
  ## set par for boxplot and legend
  op <- par(no.readonly = TRUE)
  par(mfrow = c(3, 1), mar = c(1, 2.5, 2, 1), oma = c(0, 0, 1, 0))

  ## plot unfiltered data
  rle(expr1[, order(col)]) |>
    boxplot(
      outline = FALSE, xaxt = "n", col = colors[sort(col)],
      main = "RLE plot before filtration"
    )
  abline(h = 0)
  ##plot filtered data
  rle(expr2[, order(col)]) |>
    boxplot(
      outline = FALSE, xaxt = "n", col = colors[sort(col)],
      main = "RLE plot after filtration"
    )
  abline(h = 0)
  ## add legend
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center",
         legend = unique(group_col) |> sort(),
         fill = colors, xpd = TRUE, ncol = 4
  )
  ## set par back to original
  par(op)
  return(NULL)
}

#helper: MDS plot init
plot_MDS_init <- function(expr1, expr2, group_col) {
  ## get MDS data
  mds1 <- limma::plotMDS(expr1, plot = FALSE)
  mds2 <- limma::plotMDS(expr2, plot = FALSE)

  mds1 <- data.frame(Sample = colnames(expr1),
                     Group = group_col,
                     "logFC_dim_1" = mds1$x,
                     "logFC_dim_2" = mds1$y)
  mds2 <- data.frame(Sample = colnames(expr2),
                     Group = group_col,
                     "logFC_dim_1" = mds2$x,
                     "logFC_dim_2" = mds2$y)
  p1 <- ggplot(mds1, aes(x = logFC_dim_1, y = logFC_dim_2, col = Group)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = tableau_20) +
    theme_classic()
  p2 <- ggplot(mds2, aes(x = logFC_dim_1, y = logFC_dim_2, col = Group)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = tableau_20) +
    theme_classic()
  p1 + p2 + patchwork::plot_layout(guides = "collect")
}

## matrix plot function for PCA
##---------------------------------------------------------------

##helper: Scale loadings coordinates to PCs
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
  stopifnot("dims can only be numeric vector of 2 elements!" =
              is.numeric(dims) && length(dims) == 2)
  ## basic PCA plot
  p <- ggplot() +
    geom_point(aes(x = prcomp$x[,dims[1]],
                   y = prcomp$x[,dims[2]],
                   col = group_by)) +
    ggplot2::scale_color_manual(values = tableau_20) +
    ggplot2::theme_classic()

  ## add loadings
  if (loading == TRUE) {
    ## extract top loadings
    if (is.numeric(n_loadings)) {
      pc1_genes <- -abs(prcomp$rotation[,dims[1]]) |>
        sort() |>
        head(n_loadings) |>
        names()
      pc2_genes <- -abs(prcomp$rotation[,dims[2]]) |>
        sort() |>
        head(n_loadings) |>
        names()
    } else if (is.character(n_loadings)) {
      stopifnot("Not all n_loadings are in input prcomp data" =
                  n_loadings %in% rownames(prcomp$rotation))
      pc1_genes <- pc2_genes <- prcomp$rotation[n_loadings, ] |>
        rownames()
    }

    genes2plot <- unique(pc1_genes, pc2_genes)

    ## convert coordinates for loadings
    r <- min(
      (max(prcomp$x[, dims[1]]) - min(prcomp$x[, dims[1]]) /
         (max(prcomp$rotation[, dims[1]]) - min(prcomp$rotation[, dims[1]]))),
      (max(prcomp$x[, dims[2]]) - min(prcomp$x[, dims[2]]) /
         (max(prcomp$rotation[, dims[2]]) - min(prcomp$rotation[, dims[2]])))
    )

    loadings2plot <- prcomp$rotation[genes2plot, dims] |>
      as.data.frame() |>
      dplyr::mutate(
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
      ) +
      ggplot2::labs(
        x = paste0("PC", dims[1], "(",
                   round(summary(prcomp)$importance[2, dims[1]] * 100, 2),
                   "%)"),
        y = paste0("PC", dims[2], "(",
                   round(summary(prcomp)$importance[2, dims[2]] * 100, 2),
                   "%)")
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
#' @param normalize logical, TRUE indicates raw counts data to normalize and
#'                  calculate logCPM before PCA
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
                                 normalize = FALSE,
                                 group_by = NULL,
                                 scale = TRUE,
                                 n = 4,
                                 loading = FALSE,
                                 n_loadings = 10,
                                 gene_id = "SYMBOL") {

  stopifnot(is.logical(normalize), is.vector(features),
            is.null(group_by) | is.vector(group_by),
            is.logical(scale), is.numeric(n),
            is.logical(loading), is.character(gene_id))

  ## calculate logCPM for raw counts data
  if(normalize) {
    data <- edgeR::DGEList(counts = data)
    data <- edgeR::calcNormFactors(data)
    data <- edgeR::cpm(data, log = TRUE)
  }

  ## scale and do PCA on selected features
  if (length(features) == 1 && features[1] == "all") {
    tmp <- stats::prcomp(Matrix::t(data), scale = scale)
  }else {
    id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                features, gene_id, "SYMBOL") |>
      na.omit()
    id <- id[id %in% rownames(data)]
    tmp <- stats::prcomp(Matrix::t(data[id,]), scale = scale)
    rownames(tmp$rotation) <- names(id)
  }

  ## screeplot of importance of PCs
  p0 <- summary(tmp)$importance |>
    Matrix::t() |>
    as.data.frame() |> (\(x) {
      ggplot2::ggplot(data = x,
                      ggplot2::aes(x = reorder(rownames(x), seq_len(nrow(x))),
                                   y = `Proportion of Variance`,
                                   group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "PCs") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    })()

  ## top n PCs matrix plot
  p <- list(NA)
  k <- 1
  p1 <- lapply(seq_len(n), function(x) {
    ggplot2::ggplot() +
      ggplot2::geom_text(ggplot2::aes(
        x = 1, y = 1,
        label = paste0("PC", x, "\n",
                       round(summary(tmp)$importance[2, x] * 100, 2), "%")
      ),
      size = 5) +
      ggplot2::theme_void()
  })

  ## plot multi-PCs PCA
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      p[[k]] <- plotPCAbiplot(
        prcomp = tmp,
        loading = loading,
        n_loadings = n_loadings,
        dims = c(j, i),
        group_by = group_by
      )
      p[[k]] <- p[[k]] + ggplot2::guides(col = "none")
      k <- k + 1
    }
  }

  p3 <- ggpubr::get_legend(
    ggplot() +
      geom_point(aes(x = 1, y = 1, col = group_by)) +
      scale_color_manual(values = tableau_20) +
      theme_classic()
  ) |> ggpubr::as_ggplot()  ## get legend

  ## add color for group_by factor
  if(!is.null(group_by)) {
    p0 <- (p0 + patchwork::inset_element(p3, 0.6, 0.2, 0.8, 0.9)) |>
      patchwork::patchworkGrob()
  }

  layout_mat <- matrix(NA, n, n)
  layout_mat[lower.tri(layout_mat)] <- seq_along(p)
  layout_mat <- Matrix::t(layout_mat)
  layout_mat[(ceiling(n / 2) + 1):n, seq_len(floor(n / 2))] <- length(p) + 1
  diag(layout_mat) <- (length(p) + 2):(length(p) + 1 + n)

  p <- gridExtra::arrangeGrob(grobs = c(p, list(p0), p1),
                              layout_matrix = layout_mat)
  p <- ggpubr::as_ggplot(p)
  return(p)
}


## scatter_plot function
##---------------------------------------------------------------
scatter_plot_init <- function(expr, sigs, target_group, by,
                              normalize = FALSE,
                              xint = 1, yint = 1,
                              gene_id = "SYMBOL") {
  stopifnot(is.logical(normalize), is.numeric(xint),
            is.numeric(yint), is.character(gene_id))

  ## log-normalize data if normalize = T
  if (normalize) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = TRUE)
  } else {
    expr <- expr
  }

  ## convert result gene symbols into gene IDs mapped to data rownames
  sigs <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                sigs,
                                columns = gene_id,
                                keytype = "SYMBOL")
  sigs <- sigs[!duplicated(sigs[[gene_id]]), gene_id]
  sigs <- intersect(rownames(expr), sigs)

  ## calculate median expression of sigs by aggregating cells based on by group
  group <- ifelse(grepl(target_group, by), target_group, by)
  expr <- dplyr::group_by(expr[sigs,] |> Matrix::t() |> as.data.frame(),
                          group |> factor()) |>
    dplyr::summarise_all(.funs = mean, na.rm = TRUE) |> Matrix::t()
  colnames(expr) <- expr[1,]
  expr <- expr[-1, ] |> apply(MARGIN = c(1, 2), FUN = as.numeric) |>
    as.data.frame()

  ## biplot of specified target_group vs all other types
  p <- tidyr::pivot_longer(expr, -target_group,
                           names_to = "celltype",
                           values_to = "Median Expression") |>
    ggplot2::ggplot(ggplot2::aes(x = !!sym("Median Expression"),
                                 y = !!sym(target_group))) +
    ggplot2::geom_point(alpha = 0.4, shape = 1) +
    ggplot2::geom_vline(xintercept = xint, lty = 2) +
    ggplot2::geom_hline(yintercept = yint, lty = 2) +
    ggplot2::facet_wrap(~celltype) +
    ggplot2::theme_classic()
  return(p)
}


## boxplot of expression for signatures
##---------------------------------------------------------------
exp_boxplot_init <- function(expr, sigs, target_group, by, normalize = FALSE,
                             method = "t.test", gene_id = "SYMBOL") {

  stopifnot(is.logical(normalize), is.character(gene_id))

  ## log-normalize data if normalize = T
  if (normalize) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = TRUE)
  } else {
    expr <- expr
  }

  ## convert result gene symbols into gene IDs mapped to data rownames
  sigs <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                sigs,
                                columns = gene_id,
                                keytype = "SYMBOL")
  sigs <- sigs[!duplicated(sigs[[gene_id]]), gene_id]
  sigs <- intersect(rownames(expr), sigs)

  ## calculate median expression of sigs by aggregating samples on group
  group <- ifelse(grepl(target_group, by), target_group, by)
  expr <- dplyr::group_by(expr[sigs, ] |> Matrix::t() |> as.data.frame(),
                          group |> factor()) |>
    dplyr::summarise_all(.funs = median, na.rm = TRUE) |> Matrix::t()
  colnames(expr) <- expr[1, ]
  expr <- expr[-1, ] |> apply(MARGIN = c(1, 2), FUN = as.numeric) |>
    as.data.frame()
  expr$Gene <- rownames(expr)

  ## boxplot of signature in specified target_group vs all other types
  p <- tidyr::pivot_longer(expr, -Gene,
                           names_to = "Group",
                           values_to = "Median Expression") |>
    ggplot2::ggplot(ggplot2::aes(x = !!sym("Group"),
                                 y = !!sym("Median Expression"))) +
    ggplot2::geom_boxplot(ggplot2::aes(col = Group)) +
    ggplot2::geom_point(ggplot2::aes(col = Group)) +
    ggplot2::geom_line(ggplot2::aes(group = Gene), col = "grey", alpha = .5) +
    ggpubr::stat_compare_means(label = "p.signif",
                               method = method,
                               ref.group = target_group) +
    ggplot2::labs(title = "Signatures Median Expression across Groups",
                  x = "Groups", y = "Expression") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  return(p)
}


## score boexplot
##---------------------------------------------------------------
#helper: calculate singscores
singscore_init <- function(expr, sigs, by, normalize = FALSE,
                           gene_id = "SYMBOL") {
  ## calculate log-normalized data if raw counts data
  if (normalize) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = TRUE)
  }

  ## singscore samples using given signature genes
  rank_data <- singscore::rankGenes(expr)
  UP <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                              sigs,
                              column = gene_id,
                              keytype = "SYMBOL") |>
    na.omit()

  scores <- singscore::simpleScore(rank_data, upSet = UP)
  return(scores)
}

#helper: boxplot of comparing scores
score_boxplot_init <- function(scores, by, target_group, method) {
  ## set Group for samples by given factor
  scores$Group <- ifelse(grepl(target_group, by), target_group, by)

  ## plot
  p <- ggplot2::ggplot(data = scores,
                       ggplot2::aes(x = Group, y = TotalScore)) +
    ggplot2::geom_boxplot(ggplot2::aes(col = Group)) +
    ggplot2::geom_jitter(ggplot2::aes(col = Group), alpha = 0.2) +
    ggpubr::stat_compare_means(label = "p.signif", method = method,
                               ref.group = target_group, label.y.npc = 1) +
    ggplot2::labs(title = "Scores across Groups") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}


## GSEA plot for signatures
##---------------------------------------------------------------
## gse analysis
gsea_analysis <- function(tDEG, gsets, gene_id = "SYMBOL", digits = 2) {
  gse <- lapply(names(tDEG), function(n) {
    glist <- tDEG[[n]][,"logFC"]
    names(glist) <- rownames(tDEG[[n]])
    glist <- sort(glist, decreasing = TRUE)

    ## pick common genes in glist and gsets
    gsets <- gsets[gsets[, gene_id] %in% names(glist),]

    gse <- clusterProfiler::GSEA(glist, TERM2GENE = gsets[, c("set", gene_id)])
    if(nrow(gse) == 0) {
      ms <- paste("No term was found enriched in", n, ".")
      message(ms)
      return(NULL)
    } else {
      ## only keep 2 digits for pvalue table
      gse@result$p.adjust <- signif(gse@result$p.adjust, digits = digits)
      gse@result$pvalue <- signif(gse@result$pvalue, digits = digits)
      ## save comparison info
      gse@result$group = n
      return(gse)
    }
  })
  names(gse) <- names(tDEG)
  return(gse)
}
## gseaplot
gsea_plot_init <- function(gse) {
  p <- lapply(names(gse), function(n) {
    if(is.null(gse[[n]])) {
      ms <- paste("No term was found enriched in", n, ".")
      message(ms)
      p <- grid::textGrob(ms) |> ggpubr::as_ggplot()
    } else {
      p <- enrichplot::gseaplot2(gse[[n]], geneSetID = seq_len(nrow(gse[[n]])),
                                 pvalue_table = TRUE,
                                 title = n) |>
        patchwork::wrap_plots(ncol = 1) |>
        patchwork::wrap_elements()
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
gsea_dotplot_init <- function(gse, size = "enrichmentScore") {
  res <- do.call(rbind, lapply(gse, \(x) x@result))

  p <- ggplot(res) +
    geom_point(aes(x = group, y = ID, col = -log10(p.adjust),
                   size = !!sym(size))) +
    labs(x = "Comparison", y = "Signature") +
    theme_bw() +
    scale_color_viridis_c() +
    theme(axis.text.x = element_text(angle = 90))
  return(p)
}


## heatmap for comparing signature with original markers pool
##---------------------------------------------------------------
heatmap_init <- function(expr, sigs, by, markers, normalize = FALSE,
                         scale = c("none", "row", "column"),
                         gene_id = "SYMBOL", ranks_plot = FALSE,
                         ...) {
  scale <- match.arg(scale)
  ## name sigs list
  if(is.null(names(sigs)))
    names(sigs) <- seq_along(sigs)
  idx <- which(names(sigs) == "")
  names(sigs)[idx] <- idx
  ## normalize
  if (normalize == TRUE) {
  expr <- edgeR::DGEList(expr)
  expr <- edgeR::calcNormFactors(expr, method = "TMM")
  expr <- edgeR::cpm(expr, log = TRUE)
  }
  ## plot rank instead of expression
  if (ranks_plot == TRUE)
    expr <- apply(expr, 2, rank) |> log2()
  ## scale
  if(scale == "row") {
    expr <- t(scale(t(expr)))
  } else if(scale == "column") {
    expr <- scale(expr)
  }
  ## convert sigs
  sigs <- utils::stack(sigs)
  sigs[[gene_id]] <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                           sigs$values, gene_id, "SYMBOL")
  ## kepp common genes
  idx <- sigs[[gene_id]] %in% rownames(expr)
  if(any(idx == FALSE))
    message(paste("Gene", sigs$values[!idx], "is not in data\n"))
  sigs <- sigs[idx,]

  ## subset expr
  expr <- expr[sigs[[gene_id]], ]
  rownames(expr) <- sigs$values

  ## annotation
  top_ann <- ComplexHeatmap::columnAnnotation(Group = by)
  left_ann <- ComplexHeatmap::rowAnnotation(Sigs = sigs$ind)

  ## plot
  p <- ComplexHeatmap::Heatmap(
    expr,
    column_split = by,
    row_split = sigs$ind,
    top_annotation = top_ann,
    left_annotation = left_ann,
    ...
  )

  return(p)
}


## rankdensity plot of signature
##---------------------------------------------------------------
#helper: aggregate samples of the same group to one sample with mean expression
agg_mean <- function(expr, by) {
  expr <- split(colnames(expr), by) |>
    sapply(FUN = function(x) Matrix::rowMeans(expr[,x], na.rm = TRUE))
  return(expr)
}

#helper: arrange plot matrix
# add plot_spacer() at the end for samples less than max
add_spacer <- function(x, ns) {
  if(length(x) < ns) {
    x <- append(x, paste(c("list(",
                           paste(rep("patchwork::plot_spacer()",
                                 ns - length(x)),
                                 collapse = ","),
                           ")"),
                         collapse = "") |>
                  (\(x) parse(text = x))() |> eval(),
                length(x))
  }
  x
}

#helper: make matrix plot of rankdensity
rankdensity_init <- function(expr, sigs, by, normalize = FALSE,
                             aggregate = FALSE, gene_id = "SYMBOL") {
  ## calculate log-norm expression
  if (normalize) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = TRUE)
  }else expr <- as.matrix(expr)

  ## calculate ranks of genes using singscore
  rank_data <- singscore::rankGenes(expr)
  UP <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                              sigs,
                              columns = gene_id,
                              keytype = "SYMBOL")
  UP <- UP[!duplicated(UP[[gene_id]]), gene_id] |>
    intersect(rownames(expr))

  data <- t(rank_data[UP,]/nrow(rank_data)) |>
    data.frame(Group = by)
  data$Sample <- rownames(data)

  if(aggregate == TRUE) {
    #update y-positions for barcode of up-set
    data <- data |> tidyr::pivot_longer(-c("Sample", "Group"),
                                        names_to = "Gene",
                                        values_to = "Rank") |>
      dplyr::group_by(Group) |>
      dplyr::mutate(max_den = max(density(Rank)$y))

    p <- ggplot(data, aes(x = Rank, col = Group, group = Group)) +
      geom_density() +
      geom_segment(aes(y = max_den, xend = Rank, yend = max_den*1.1)) +
      facet_wrap(~Group, scales = "free") +
      theme_bw()
  }else {
    #update y-positions for barcode of up-set
    data <- data |> tidyr::pivot_longer(-c("Sample", "Group"),
                                        names_to = "Gene",
                                        values_to = "Rank") |>
      dplyr::group_by(Group) |>
      dplyr::mutate(y_s = as.numeric(factor(Sample))) |>
      dplyr::mutate(y_s = -y_s/max(y_s), y_e = y_s - 0.1)

    p <- ggplot(data, aes(x = Rank, col = Group, group = Sample)) +
      geom_density() +
      geom_segment(aes(y = y_s, xend = Rank, yend = y_e)) +
      facet_wrap(~Group, scales = "free") +
      theme_bw()
  }

  return(p)
}



## scatter plot of HDBSCAN cluster result
##---------------------------------------------------------------
scatter_hdb_cl <- function(sig_matrix, markers_list) {

  stopifnot("sig_matrix must be a matrix or data.frame!" = is.matrix(sig_matrix) | is.data.frame(sig_matrix),
            "markers_list is not a list!" = is.list(markers_list))

  ml <- utils::stack(markers_list) |> (\(x) split(x$ind, x$values))()

  p <- lapply(rownames(sig_matrix), function(m) {
    p <- ggplot2::ggplot(data.frame(
      Gene = m,
      Expression = sig_matrix[m,] |> as.numeric(),
      Type = ifelse(colnames(sig_matrix) %in% ml[[m]],
                    colnames(sig_matrix),
                    "")
      ),
      ggplot2::aes(x = Gene, y = Expression,
                   col = (Type != ""))) +
      ggrepel::geom_label_repel(ggplot2::aes(label = Type)) +
      ggplot2::geom_point() +
      ggplot2::labs(col = "is.marker") +
      ggplot2::facet_wrap(~Gene) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red",
                                             "FALSE" = "black")) +
      ggplot2::theme_classic()
  })

  names(p) <- rownames(sig_matrix)

  return(p)
}

utils::globalVariables(c("median", "Gene", "Expression", "Group", "Type",
                         "target_group", "group_col", "TotalScore",
                         "Proportion of Variance", "Rank", "Sample",
                         "max_den", "y_s", "y_e", "x_end", "y_end",
                         "group", "ID", "NES", "logFC_dim_1",
                         "logFC_dim_2", "logcounts"))
