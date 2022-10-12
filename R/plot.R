#' @import ggplot2 ggpubr ggplotify
NULL

## plot functions for comparing unfiltered and filtered data
##---------------------------------------------------------------
#helper: logCPM density plot between filtered and unfiltered
plot_density <- function(dge, ID, keep = TRUE, counts = TRUE, filter = 10) {
  ## set colors for each cell type
  col <- dge$samples[[ID]] |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = ggthemes::tableau_color_pal("Tableau 20")(20))(unique(dge$samples[[ID]]) |> length())

  ## set par for plot to make 2 plots in one page
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))

  ## calculate logCPM
  if(counts) {
    M <- median(dge$samples$lib.size) * 1e-6
    L <- mean(dge$samples$lib.size) * 1e-6
    lcpm.1 <- edgeR::cpm(dge, log = T)
    lcpm.2 <- edgeR::cpm(dge[keep,, keep.lib.sizes = F], log = T)
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
  plot(density(lcpm.1[,1]), col = colors[col[1]],
       main = "Unfiltered", xlab = ifelse(counts,"logCPM","logcounts"),
       ylim = c(0, max(density(lcpm.1)$y)*1.2))
  abline(v = abl, lty = 3)
  t <- lapply(2:ncol(dge), function(i) lines(density(lcpm.1[,i]),
                                             col = colors[col[i]]))
  legend("topright", legend = unique(dge$samples[[ID]]),
         text.col = colors[unique(col)], bty = "n", inset = c(-0.3, 0))

  ## plot filtered data
  plot(density(lcpm.2[,1]), col = colors[col[1]],
       main = "Filtered", xlab = ifelse(counts,"logCPM","logcounts"),
       ylim = c(0, max(density(lcpm.2)$y)*1.2))
  abline(v = abl, lty = 3)
  t <- lapply(2:ncol(dge), function(i) lines(density(lcpm.2[,i]),
                                             col = colors[col[i]]))
  legend("topright", legend = unique(dge$samples[[ID]]),
         text.col = colors[unique(col)], bty = "n", inset = c(-0.3, 0))

  mtext("Input Dataset", side = 3, line = 0, outer = T)
  ## set par back to original setting
  par(op)
}

#helper: calculate RLE
rle <- function(expr_mat, log = F) {
  if (log) {
    geomeans <- Matrix::rowMeans(log1p(expr_mat), na.rm = T)
  } else {
    geomeans <- Matrix::rowMeans(expr_mat, na.rm = T)
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
plot_rle <- function(dge, ID, keep = TRUE, counts = TRUE) {
  ## set colors for each cell type
  col <- dge$samples[[ID]] |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = ggthemes::tableau_color_pal("Tableau 20")(20))(unique(dge$samples[[ID]]) |> length())

  ## calculate logCPM
  if(counts) {
    lcpm.1 <- edgeR::cpm(dge, log = T)
    lcpm.2 <- edgeR::cpm(dge[keep,, keep.lib.sizes = F], log = T)
  }else {
    # lcpm.1 <- log1p(dge$counts)
    # lcpm.2 <- log1p(dge$counts[keep,])
    lcpm.1 <- dge$counts
    lcpm.2 <- dge$counts[keep,]
  }

  ## set par for boxplot and legend
  op <- par(no.readonly = T)
  par(mfrow = c(3, 1), mar = c(1, 2.5, 2, 1), oma = c(0, 0, 1, 0))

  ## plot unfiltered data
  rle(lcpm.1[, order(col)]) |>
    boxplot(
      outline = F, xaxt = "n", col = sort(col) |> (\(x)colors[x])(),
      main = "RLE plot before filtration"
    )
  abline(h = 0)
  ##plot filtered data
  rle(lcpm.2[, order(col)]) |>
    boxplot(
      outline = F, xaxt = "n", col = sort(col) |> (\(x)colors[x])(),
      main = "RLE plot after filtration"
    )
  abline(h = 0)
  ## add legend
  plot(0, type = "n", axes = F, xlab = "", ylab = "")
  legend("center",
         legend = unique(dge$samples[[ID]]) |> sort(),
         fill = colors, xpd = T, ncol = 4
  )
  ## set par back to original
  par(op)
}

#helper: MDS plot before and after filtration
plot_MDS <- function(dge, ID, keep = TRUE, counts = TRUE) {
  ## set colors for each cell type
  col <- dge$samples[[ID]] |> as.factor() |> as.numeric()
  colors <- colorRampPalette(colors = ggthemes::tableau_color_pal("Tableau 20")(20))(unique(dge$samples[[ID]]) |> length())

  ## calculate logCPM
  if(counts) {
    lcpm.1 <- edgeR::cpm(dge, log = T)
    lcpm.2 <- edgeR::cpm(dge[keep,, keep.lib.sizes = F], log = T)
  }else {
    # lcpm.1 <- log1p(dge$counts)
    # lcpm.2 <- log1p(dge$counts[keep,])
    lcpm.1 <- dge$counts
    lcpm.2 <- dge$counts[keep,]
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
  plot(0, type = "n", axes = F, xlab = "", ylab = "")
  legend("center",
         legend = unique(dge$samples[[ID]]) |> sort(),
         fill = colors, xpd = T
  )
  mtext("Input Data", side = 3, line = 0, outer = T)
  ## set par back to original
  par(op)
}


## matrix plot function for PCA
##---------------------------------------------------------------
pca_matrix_plot_init <- function(data,
                                 features = "all",
                                 is.counts = TRUE,
                                 group_by = NULL,
                                 scale = TRUE,
                                 n = 4,
                                 loading = FALSE,
                                 gene_id = "SYMBOL") {

  stopifnot(is.logical(is.counts), is.vector(features),
            is.null(group_by) | is.vector(group_by),
            is.logical(scale), is.numeric(n),
            is.logical(loading), is.character(gene_id))

  ## calculate logCPM for raw counts data
  if(is.counts){
    data <- edgeR::DGEList(counts = data)
    data <- edgeR::calcNormFactors(data)
    data <- edgeR::cpm(data, log = T)
  }

  ## scale and do PCA on selected features
  if (length(features) == 1 & features[1] == "all") {
    tmp <- prcomp(Matrix::t(data), scale = scale)
  }else{
    features <- AnnotationDbi::select(org.Hs.eg.db, features, gene_id, "SYMBOL")
    id <- match(features[[gene_id]], rownames(data)) |> na.omit()
    tmp <- prcomp(Matrix::t(data[id,]), scale = scale)
  }

  ## screeplot of importance of PCs
  p0 <- summary(tmp)$importance |>
    Matrix::t() |>
    as.data.frame() |> (\(x) {
      ggplot2::ggplot(data = x,
                      ggplot2::aes(x = reorder(rownames(x), 1:nrow(x)),
                                   y = `Proportion of Variance`, group = 1)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "PCs") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    })()

  ## top n PCs matrix plot
  p <- list(NA)
  k <- 1
  p1 <- lapply(1:n, function(x) {
    ggplot2::ggplot() +
      ggplot2::geom_text(ggplot2::aes(
        x = 1, y = 1,
        label = paste("PC", x, round(summary(tmp)$importance[2, x] * 100, 2), "%")
      ), size = 5) +
      ggplot2::theme_void()
  })
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (is.null(group_by)) {
        Group <- data.frame(group = rep("1", ncol(data)))
      } else {
        Group <- data.frame(group = group_by)
      }
      p[[k]] <- ggplot2::autoplot(tmp, loadings = loading,
                                  loadings.label = loading,
                                  loadings.col = "grey",
                                  loadings.label.col = "black",
                                  loadings.label.size  = 2,
                                  data = Group, colour = 'group',
                                  x = j, y = i) +
        ggplot2::labs(
          x = paste0("PC", j, "(", round(tmp$importance[2, j] * 100, 2), "%)"),
          y = paste0("PC", i, "(", round(tmp$importance[2, i] * 100, 2), "%)")
        ) +
        ggthemes::scale_color_tableau("Tableau 20") +
        ggplot2::theme_classic()
      p3 <- ggpubr::get_legend(p[[k]]) |> ggpubr::as_ggplot()
      p[[k]] <- p[[k]] + ggplot2::guides(col = "none")
      k <- k + 1
    }
  }

  ## add color for group_by factor
  if(!is.null(group_by))
    p0 <- (p0 + patchwork::inset_element(p3, 0.6, 0.2, 0.8, 0.9)) |>
    patchwork::patchworkGrob()

  layout_mat <- matrix(NA, n, n)
  layout_mat[lower.tri(layout_mat)] <- 1:length(p)
  layout_mat <- Matrix::t(layout_mat)
  layout_mat[(ceiling(n / 2) + 1):n, 1:floor(n / 2)] <- length(p) + 1
  diag(layout_mat) <- (length(p) + 2):(length(p) + 1 + n)

  p <- gridExtra::arrangeGrob(grobs = c(p, list(p0), p1),
                              layout_matrix = layout_mat)
  p <- ggplotify::as.ggplot(p)
  return(p)
}


## scatter_plot function
##---------------------------------------------------------------
scatter_plot_init <- function(expr, sigs, type, by,
                              counts = TRUE,
                              xint = 1, yint = 1,
                              gene_id = "SYMBOL") {
  stopifnot(is.logical(counts), is.numeric(xint),
            is.numeric(yint), is.character(gene_id))

  ## log-normalize data if counts = T
  if (counts) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = T)
  } else {
    expr <- expr
  }

  ## convert result gene symbols into gene IDs mapped to data rownames
  sigs <- AnnotationDbi::select(org.Hs.eg.db, sigs,
                                columns = gene_id,
                                keytype = "SYMBOL")
  sigs <- sigs[!duplicated(sigs[[gene_id]]), gene_id]
  sigs <- intersect(rownames(expr), sigs)

  ## calculate median expression of sigs by aggregating cells based on by group
  group <- ifelse(grepl(type, by), type, by)
  expr <- dplyr::group_by(expr[sigs,] |> Matrix::t() |> as.data.frame(),
                          group |> factor()) |>
    dplyr::summarise_all(.funs = mean, na.rm = T) |> Matrix::t()
  colnames(expr) <- expr[1,]
  expr <- expr[-1,] |> apply(MARGIN = c(1,2), FUN = as.numeric) |>
    as.data.frame()

  ## biplot of specified type vs all other types
  p <- tidyr::gather(expr, "celltype", "Median Expression", -type) |>
    ggplot2::ggplot(ggplot2::aes_q(x = as.name("Median Expression"),
                                   y = as.name(type))) +
    ggplot2::geom_point(alpha = 0.4, shape = 1) +
    ggplot2::geom_vline(xintercept = xint, lty = 2) +
    ggplot2::geom_hline(yintercept = yint, lty = 2) +
    ggplot2::facet_wrap(~celltype) +
    ggplot2::theme_classic()
  return(p)
}


## boxplot of expression for signatures
##---------------------------------------------------------------
exp_boxplot_init <- function(expr, sigs, type, by, counts = TRUE,
                             method = 't.test', gene_id = "SYMBOL") {

  stopifnot(is.logical(counts), is.character(gene_id))

  ## log-normalize data if counts = T
  if (counts) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = T)
  } else {
    expr <- expr
  }

  ## convert result gene symbols into gene IDs mapped to data rownames
  sigs <- AnnotationDbi::select(org.Hs.eg.db, sigs,
                                columns = gene_id,
                                keytype = "SYMBOL")
  sigs <- sigs[!duplicated(sigs[[gene_id]]), gene_id]
  sigs <- intersect(rownames(expr), sigs)

  ## calculate median expression of sigs by aggregating samples based on by group
  group <- ifelse(grepl(type, by), type, by)
  expr <- dplyr::group_by(expr[sigs,] |> Matrix::t() |> as.data.frame(),
                          group |> factor()) |>
    dplyr::summarise_all(.funs = median, na.rm = T) |> Matrix::t()
  colnames(expr) <- expr[1,]
  expr <- expr[-1,] |> apply(MARGIN = c(1,2), FUN = as.numeric) |>
    as.data.frame() |> tibble::rownames_to_column(var = "Gene")

  ## boxplot of signature in specified type vs all other types
  p <- tidyr::gather(expr, "Group", "Median Expression", -Gene) |>
    ggplot2::ggplot(ggplot2::aes_q(x = as.name("Group"),
                                   y = as.name("Median Expression"))) +
    ggplot2::geom_boxplot(ggplot2::aes(col = Group)) +
    ggplot2::geom_point(ggplot2::aes(col = Group)) +
    ggplot2::geom_line(ggplot2::aes(group = Gene), col = "grey", alpha = .5) +
    ggpubr::stat_compare_means(label = "p.signif",
                               method = method,
                               ref.group = type) +
    ggplot2::labs(title = "Signatures Median Expression across Groups",
                  x = "Types", y = "Expression") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  return(p)
}


## score boexplot
##---------------------------------------------------------------
#helper: calculate singscores
singscore_init <- function(expr, sigs, by, counts = TRUE,
                           gene_id = "SYMBOL") {
  ## calculate log-normalized data if raw counts data
  if (counts) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = T)
  }

  ## singscore samples using given signature genes
  rank_data <- singscore::rankGenes(expr)
  UP <- AnnotationDbi::select(org.Hs.eg.db, sigs,
                              columns = gene_id,
                              keytype = "SYMBOL")

  UP <- UP[!duplicated(UP[[gene_id]]), gene_id]
  scores <- singscore::simpleScore(rank_data, upSet = UP)
  return(scores)
}

#helper: boxplot of comparing scores
score_boxplot_init <- function(scores, by, type, method) {
  ## set Group for samples by given factor
  scores$Group <- ifelse(grepl(type, by), type, by)

  ## plot
  p <- ggplot2::ggplot(data = scores,
                       ggplot2::aes(x = Group, y = TotalScore)) +
    ggplot2::geom_boxplot(ggplot2::aes(col = Group)) +
    ggplot2::geom_jitter(ggplot2::aes(col = Group), alpha = 0.2) +
    ggpubr::stat_compare_means(label = "p.signif", method = method,
                               ref.group = type, label.y.npc = 1) +
    ggplot2::labs(title = "Scores across Groups") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}


## GSEA plot for signatures
##---------------------------------------------------------------
gsea_plot_init <- function(tDEG, gsets, gene_id = "SYMBOL", digits = 2) {
  p <- lapply(names(tDEG), function(n){
    glist <- tDEG[[n]][,"logFC"]
    names(glist) <- rownames(tDEG[[n]])
    glist <- sort(glist, decreasing = T)

    gse <- clusterProfiler::GSEA(glist, TERM2GENE = gsets[,c("set", gene_id)])
    if(nrow(gse) == 0) {
      cat(paste("No term was found enriched in", n, "."))
      p <- grid::textGrob(paste("No term was found enriched in", n, "."))
    } else {
      p <- enrichplot::gseaplot2(gse, geneSetID = 1:nrow(gse),
                                 pvalue_table = F,
                                 title = paste("GSEA with signatures in", n))
      p <- p - patchwork::inset_element(gridExtra::tableGrob(
        gse[,c("pvalue", "p.adjust")] |> signif(digits),
        theme = gridExtra::ttheme_default(base_size = 5)
      ),
      left = 0.6, bottom = 0.8, right = 1, top = 1)
    }
  })
  names(p) <- names(tDEG)
  return(p)
}



## heatmap for comparing signature with original markers pool
##---------------------------------------------------------------
heatmap_init <- function(expr, sigs, by, markers, counts = TRUE,
                         scale = "none", min_max = FALSE,
                         gene_id = "SYMBOL", ranks_plot = FALSE,
                         col = grDevices::colorRampPalette(colors = c(rgb(1, 1, 0),
                                                                      "red"))(256)) {
  if (counts) {
    expr <- edgeR::DGEList(expr)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = T)
    expr <- expr[,order(by)] ## order samples based on by factor
  } else {
    expr <- expr[, order(by)]
  }

  ## convert markers and sigs into gene IDs matched to expr rownames
  tmp <- AnnotationDbi::select(org.Hs.eg.db, rownames(expr),
                               columns = "SYMBOL", keytype = gene_id)

  ## get index of matched genes between markers/sigs and expr data
  index.1 <- tmp[[gene_id]][tmp$SYMBOL %in% markers] |>
    unique() |> sort() |> match(rownames(expr)) |> na.omit()
  index.2 <- tmp[[gene_id]][tmp$SYMBOL %in% sigs] |>
    unique() |> sort() |> match(rownames(expr)) |> na.omit()

  ## set color for annotation bar
  anno_col <- data.frame(cell_type = sort(by))
  rownames(anno_col) <- colnames(expr)
  ann_colors <- list(cell_type = colorRampPalette(
    ggthemes::tableau_color_pal("Tableau 20")(min(20,
                                                  length(unique(by)))))(length(unique(by)))
  )
  names(ann_colors$cell_type) <- sort(unique(by))

  ## if to convert expression matrix to rank matrix
  if (ranks_plot)
    expr <- apply(expr, 2, rank) |> log2()
  ## scale genes expression across samples using min_max(z-score) way
  if (min_max)
    expr <- apply(expr, 1, function(x) (x - min(x)) / (max(x) - min(x))) |> Matrix::t()
  ## heatmap for markers pool genes
  p1 <- pheatmap::pheatmap(expr[c(setdiff(index.1, index.2), index.2),],
                           silent = T, scale = scale,
                           border_color = NA, annotation_col = anno_col,
                           show_rownames = F,
                           gaps_col = table(by) |> cumsum(),
                           show_colnames = F, main = "Original Markers Pool",
                           gaps_row = length(index.1) - length(index.2),
                           cluster_cols = F, cluster_rows = F, color = col,
                           annotation_colors = ann_colors
  ) |> ggplotify::as.ggplot()
  ## heatmap for final signature genes
  p2 <- pheatmap::pheatmap(expr[index.2,],
                           silent = T, scale = scale,
                           border_color = NA, annotation_col = anno_col,
                           show_rownames = F,
                           gaps_col = table(by) |> cumsum(),
                           show_colnames = F, main = "Screened Signature",
                           cluster_cols = F, cluster_rows = F, color = col,
                           annotation_colors = ann_colors
  ) |> ggplotify::as.ggplot()

  return(list(p1, p2))
}


## rankdensity plot of signature
##---------------------------------------------------------------
#helper: aggregate samples of the same group to one sample with mean expression
agg_mean <- function(expr, by) {
  expr <- split(colnames(expr), by) |>
    sapply(FUN = function(x) Matrix::rowMeans(expr[,x], na.rm = T))
  return(expr)
}

#helper: arrange plot matrix, add plot_spacer() at the end for samples less than max
add_spacer <- function(x, ns){
  if(length(x) < ns) {
    x <- append(x, paste(c("list(",
                           paste(rep("patchwork::plot_spacer()",ns-length(x)),
                                 collapse = ","),")"),collapse = "") |>
                  (\(x) parse(text = x))() |> eval(),
                length(x))
  }
  x
}

#helper: make matrix plot of rankdensity
rankdensity_init <- function(expr, sigs, by, counts = TRUE,
                             aggregate = FALSE, gene_id = "SYMBOL") {
  ## calculate log-norm expression
  if (counts) {
    expr <- edgeR::DGEList(expr,
                           group = by)
    expr <- edgeR::calcNormFactors(expr, method = "TMM")
    expr <- edgeR::cpm(expr, log = T)
  }

  ## aggregate samples based on by group
  if(aggregate) {
    expr <- agg_mean(expr = expr, by = by)
    by <- colnames(expr)
    # expr <- expr[, order(by)]
  }else expr <- expr[, order(by)]

  ## calculate ranks of genes using singscore
  rank_data <- singscore::rankGenes(expr)
  UP <- AnnotationDbi::select(org.Hs.eg.db, sigs,
                              columns = gene_id,
                              keytype = "SYMBOL")
  UP <- UP[!duplicated(UP[[gene_id]]), gene_id]
  NK_scores <- singscore::simpleScore(rank_data, upSet = UP)

  ## plot rankdensity
  p <- list()
  for (j in 1:ncol(rank_data)) {
    p[[j]] <- singscore::plotRankDensity(rank_data[, j, drop = F],
                                         upSet = UP,
                                         textSize = 1
    )
  }
  tags <- lapply(sort(by) |> unique(),
                 function(x) {
                   ggplot2::ggplot() +
                     ggplot2::geom_text(ggplot2::aes(x = 1, y = 1, label = x),
                                        size = 5) +
                     ggplot2::theme_void()
                 })
  # layout_mat <- split(1:ncol(expr), data[[i]]@phenoData@data[[ID[i]]] %>% sort()) %>% do.call(rbind, args = .)
  # layout_mat[t(apply(layout_mat, 1, duplicated))] <- "#"
  # layout_mat <- cbind(1:length(unique(data[[i]]@phenoData@data[[ID[i]]])) + ncol(expr), layout_mat)
  # layout_mat <- apply(layout_mat, 1, paste, collapse = "") %>%
  #   paste(collapse = "\n")  ## convert matrix to design for patchwork

  ## set num of column for matrix plot
  ns <- table(by) |> max()

  ## plot
  p <- split(p, sort(by)) |>
    lapply(add_spacer, ns = ns) |>
    do.call(what = c)
  p <- patchwork::wrap_plots(p, ncol = ns)
  tags <- patchwork::wrap_plots(tags, ncol = 1)
  p <- (tags | p) + patchwork::plot_layout(widths = c(1,ns))

  return(p)
}



## scatter plot of HDBSCAN cluster result
##---------------------------------------------------------------
scatter_hdb_cl <- function(sig_matrix, minPts = 2, ...) {

  stopifnot(is.matrix(sig_matrix) | is.data.frame(sig_matrix))

  p <- BiocParallel::bplapply(rownames(sig_matrix), function(m) {
    ## cluster cell types based on each gene profile
    cl <- dbscan::hdbscan(cbind(sig_matrix[m,], sig_matrix[m,]),
                          minPts = minPts, ...)

    flag <- cl$cluster[which.max(sig_matrix[m,])]  ## which cluster is the max expression
    if(flag == "0"){
      ## if max is noise point, assign this gene to the unique cell type
      type <- colnames(sig_matrix)[which.max(sig_matrix[m,])]
    }else{
      ## if max is a cluster, assign this gene to all cell types in the cluster
      type <- colnames(sig_matrix)[cl$cluster == flag]
    }

    p <- ggplot2::ggplot(data.frame(Gene = m,
                                    Expression = sig_matrix[m,],
                                    Cluster = cl$cluster,
                                    Type = ifelse(colnames(sig_matrix) %in% type,
                                                  colnames(sig_matrix),
                                                  "")),
                         ggplot2::aes(x = Gene , y = Expression,
                                      col = (Type != ""))) +
      ggrepel::geom_label_repel(aes(label = Type)) +
      ggplot2::geom_point() +
      ggplot2::labs(col = "is.marker") +
      ggplot2::facet_wrap(~Gene) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      ggplot2::theme_classic()
  })

  names(p) <- rownames(sig_matrix)

  return(p)
}

