## =========================================================================##
## =========================================================================##
##                       Generic methods definitions                        ##
## =========================================================================##
## =========================================================================##


## ===========================================================================
## Generic for data processing
## ---------------------------------------------------------------------------
#' @title process data
#'
#' @description filter low expression genes, normalize data by 'TMM' and apply
#'     [limma::voom()], [limma::lmFit()] and [limma::treat()] on normalized data
#'
#' @inheritParams de_analysis
#' @param data expression object
#' @param slot character, specify which slot to use only for DGEList, sce or
#'             seurat object, optional, default 'counts'
#' @param ... params for [voom_fit_treat()]
#'
#' @return A DGEList containing vfit by [limma::voom()] (if normalize = TRUE)
#'         and tfit by [limma::treat()]
#'
#' @include DE_functions.R plot.R
#' @export
#'
#' @examples
#' data("im_data_6")
#' proc_data <- process_data(
#'   im_data_6,
#'   group_col = "celltype:ch1",
#'   target_group = "NK"
#' )
setGeneric(
  "process_data",
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           filter = c(10, 10),
           lfc = 0,
           p = 0.05,
           markers = NULL,
           gene_id = "SYMBOL",
           slot = "counts",
           ...) {
    standardGeneric("process_data")
  }
)


## ===========================================================================
## Generic for getting DEGs
## ---------------------------------------------------------------------------
#' @title Get differentially expressed genes by comparing specified groups
#'
#' @description This function uses edgeR and limma to get 'UP' and 'DOWN' DEG
#'   lists, for multiple comparisons, DEGs can be obtained from intersection of
#'   all DEGs or by using product of p value ranks for multiple
#'   comparisons. Filter out low expressed genes and extract DE genes by using
#'   limma::voom and limma::treat, and also create an object `proc_data` to
#'   store processed data.
#'
#' @inheritParams process_data
#' @inheritParams select_sig
#' @param data expression object
#' @param group_col vector or character, specify the group factor or column name of
#'           coldata for DE comparisons
#' @param slot character, specify which slot to use only for DGEList, sce or
#'             seurat object, optional, default 'counts'
#' @param batch vector of column name(s) or dataframe, specify the batch effect
#'              factor(s), default NULL
#' @param ... params for [process_data()] and [select_sig()]
#'
#' @return A list of 'UP', 'DOWN' gene set of all differentially expressed
#'         genes, and a DGEList 'proc_data' containing data after process
#'         (filtration, normalization and voom fit). Both 'UP' and 'DOWN' are
#'         ordered by rank product or 'Rank' variable if keep.top is NULL
#'
#' @examples
#' data("im_data_6")
#' DEGs <- get_degs(im_data_6,
#'   group_col = "celltype:ch1",
#'   target_group = "NK", gene_id = "ENSEMBL"
#' )
#'
#' @include DE_functions.R plot.R
#' @export
setGeneric(
  "get_degs",
  function(data,
           group_col,
           target_group,
           normalize = TRUE,
           feature_selection = c("auto", "rankproduct", "none"),
           slot = "counts",
           batch = NULL,
           ...) {
    standardGeneric("get_degs")
  }
)


## ===========================================================================
## Generic for getting DEG tables
## ---------------------------------------------------------------------------
#' @title Get DE analysis result table(s) with statistics
#'
#' @description This function uses edgeR and limma to get DE analysis results
#'   lists for multiple comparisons. Filter out low expressed genes and obtain
#'   DE statistics by using limma::voom and limma::treat, and also create an
#'   object `proc_data` to store processed data.
#'
#' @inheritParams get_degs
#' @param ... params for function [de_analysis()]
#'
#' @return A list of DE result table of all comparisons.
#'
#' @include DE_functions.R plot.R
#' @examples
#' data("im_data_6")
#' DE_tables <- get_de_table(im_data_6, group_col = "celltype:ch1", target_group = "NK")
#'
#' @export
setGeneric(
  "get_de_table",
  function(data,
           group_col,
           target_group,
           slot = "counts",
           ...) {
    standardGeneric("get_de_table")
  }
)


## ===========================================================================
## Generic for wrapper functions
## ---------------------------------------------------------------------------

#' Filter specific cell type signature genes against other subsets.
#'
#' Specify the signature of the subset matched 'target_group' against other subsets,
#' either "union", "intersect" or "RRA" can be specified when input is a list
#' of datasets to integrate the signatures into one.
#'
#' @inheritParams process_data
#' @inheritParams select_sig
#' @param data An expression data or a list of expression data objects
#' @param group_col vector or character, specify the group factor or column name of
#'           coldata for DE comparisons
#' @param dir character, could be 'UP' or 'DOWN' to use only up- or
#'            down-expressed genes
#' @param comb 'RRA' or Fun for combining sigs from multiple datasets, keep all
#'             passing genes or only intersected genes, could be `union` or
#'             `intersect` or `setdiff` or customized Fun, or could be 'RRA' to
#'             use Robust Rank Aggregation method for integrating multi-lists
#'             of sigs, default 'union'
#' @param filter (list of) vector of 2 numbers, filter condition to remove low
#'               expression genes, the 1st for min.counts (if normalize = TRUE)
#'               or CPM/TPM (if normalize = FALSE), the 2nd for samples size 'large.n'
#' @param s_thres num, threshold of score if comb = 'RRA'
#' @param ... other params for [get_degs()]
#'
#' @return a vector of gene symbols
#'
#' @examples
#' data("im_data_6", "nk_markers")
#' sigs <- filter_subset_sig(im_data_6, "celltype:ch1", "NK",
#'   markers = nk_markers$HGNC_Symbol,
#'   gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric(
  "filter_subset_sig",
  function(data,
           group_col,
           target_group,
           markers = NULL,
           normalize = TRUE,
           dir = "UP",
           gene_id = "SYMBOL",
           feature_selection = c("auto", "rankproduct", "none"),
           comb = union,
           filter = c(10, 10),
           s_thres = 0.05,
           ...) {
    standardGeneric("filter_subset_sig")
  }
)

#' @title Remove markers with high signal in background data.
#'
#' @description Specify signatures against specific tissues or cell lines by
#'   removing genes with high expression in the background.
#'
#' @param sig_data log-transformed expression object, can be matrix or DGEList,
#'                 as signal data
#' @param bg_data 'CCLE' or log-transformed expression object as background data
#' @param markers vector, a vector of gene names, listed the gene symbols to be
#'                filtered. Must be gene SYMBOLs
#' @param s_group_col vector or character, to specify the group of signal
#'                    target_groups, or column name of group, default NULL
#' @param s_target_group pattern, specify the target group of interest in
#'                       sig_data, default NULL
#' @param b_group_col vector or character, to specify the group of background
#'                     target_groups, or column name of [depmap::depmap_metadata()],
#'                     e.g. 'primary_disease', default NULL
#' @param b_target_group pattern, specify the target_group of interest in bg_data,
#'                       e.g. 'colorectal', default NULL
#' @param snr num, the cutoff of SNR to screen markers which are not or lowly
#'            expressed in bg_data
#' @param ... params for [grep()] to find matched cell lines in bg_data
#' @param filter NULL or a vector of 2 num, filter condition to remove low
#'               expression genes in bg_data, the 1st for logcounts,
#'               the 2nd for samples size
#' @param gene_id character, specify the gene ID type of rownames of expression
#'                data, could be one of 'ENSEMBL', 'SYMBOL', 'ENTREZ'...,
#'                default 'SYMBOL'
#' @param s_slot character, specify which slot to use of DGEList, sce or seurat
#'               object for sig_data, optional, default 'counts'
#' @param b_slot character, specify which slot to use of DGEList, sce or seurat
#'               object for bg_data, optional, default 'counts'
#' @param ccle_tpm ccle_tpm data from [depmap::depmap_TPM()], only used when
#'                 data = 'CCLE', default NULL
#' @param ccle_meta ccle_meta data from [depmap::depmap_metadata()], only used
#'                  when data = 'CCLE', default NULL
#'
#' @return a vector of genes after filtration
#'
#' @examples
#' data("im_data_6", "nk_markers", "ccle_crc_5")
#' remove_bg_exp(
#'   sig_data = Biobase::exprs(im_data_6),
#'   bg_data = ccle_crc_5,
#'   im_data_6$`celltype:ch1`, "NK", ## for sig_data
#'   "cancer", "CRC", ## for bg_data
#'   markers = nk_markers$HGNC_Symbol[40:50],
#'   filter = c(1, 2),
#'   gene_id = c("ENSEMBL", "SYMBOL")
#' )
#'
#' @export
setGeneric(
  "remove_bg_exp",
  function(sig_data,
           bg_data = "CCLE",
           markers,
           s_group_col = NULL,
           s_target_group = NULL,
           b_group_col = NULL,
           b_target_group = NULL,
           snr = 1,
           ...,
           filter = NULL,
           gene_id = "SYMBOL",
           s_slot = "counts",
           b_slot = "counts",
           ccle_tpm = NULL,
           ccle_meta = NULL) {
    standardGeneric("remove_bg_exp")
  }
)


## ===========================================================================
## Generic for subsetting (MSigDB) geneset collection
## ---------------------------------------------------------------------------
#' Collect genes from MSigDB or provided GeneSetCollection.
#'
#' Collect gene sets from MSigDB or given GeneSetCollection, of which the gene-set
#' names are matched to the given regex pattern by using [grep()] function.
#' By setting cat and subcat, matching can be constrained in the union of given
#' categories and subcategories if gsc = 'msigdb'.
#'
#' @param gsc 'msigdb' or GeneSetCollection to be searched
#' @param pattern pattern pass to [grep()], to match the MsigDB gene-set name of
#'                interest, e.g. 'NATURAL_KILLER_CELL_MEDIATED'
#' @param cat character, stating the category(s) to be retrieved.
#'            The category(s) must be one from [msigdb::listCollections()],
#'            see details in [msigdb::subsetCollection()]
#' @param subcat character, stating the sub-category(s) to be retrieved.
#'               The sub-category(s) must be one from
#'               [msigdb::listSubCollections()], see details in
#'               [msigdb::subsetCollection()]
#' @param species character, species of interest, can be 'hs' or 'mm'
#' @param id a character, representing the ID type to use ("SYM" for gene
#'           SYMBOLs and "EZID" for ENTREZ IDs)
#' @param version a character, stating the version of MSigDB to be retrieved
#'                (should be >= 7.2). See [msigdb::getMsigdbVersions()].
#' @param ... params for [grep()], used to match pattern to gene-set names
#'
#' @return A GeneSet object containing all matched gene-sets in MSigDB
#' @export
#'
#' @examples
#' data("msigdb_gobp_nk")
#' get_gsc_sig(
#'   gsc = msigdb_gobp_nk,
#'   pattern = "natural_killer_cell_mediated",
#'   subcat = "GO:BP",
#'   ignore.case = TRUE
#' )
setGeneric(
  "get_gsc_sig",
  function(gsc = "msigdb",
           pattern,
           cat = NULL,
           subcat = NULL,
           species = c("hs", "mm"),
           id = c("SYM", "EZID"),
           version = msigdb::getMsigdbVersions(),
           ...) {
    standardGeneric("get_gsc_sig")
  }
)

## ===========================================================================
## Generic for gene-set conversion
## ---------------------------------------------------------------------------
#' Convert gene-set list into GeneSetCollection
#'
#' @param ... vector of genes or list of genes
#'
#' @return GeneSetCollection
#' @export
#'
#' @examples
#' data("msigdb_gobp_nk")
#' gls2gsc(GSEABase::geneIds(msigdb_gobp_nk[1:3]))
setGeneric(
  "gls2gsc",
  function(...) {
    standardGeneric("gls2gsc")
  }
)


## ===========================================================================
## Generic for pseudo bulking
## ---------------------------------------------------------------------------
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


## ===========================================================================
## Generic for visualization
## ---------------------------------------------------------------------------

#' @title Make a matrix plot of PCA with top PCs
#'
#' @inheritParams pca_matrix_plot_init
#' @param data expression data, can be matrix, eSet, seurat...
#' @param slot character, specify the slot name of expression to be used,
#'             optional
#'
#' @return matrix plot of PCA
#'
#' @examples
#' data("im_data_6")
#' pca_matrix_plot(data = im_data_6, scale = FALSE)
#'
#' @include plot.R
#' @export
setGeneric(
  "pca_matrix_plot",
  function(data,
           features = "all",
           slot = "counts",
           group_by = NULL,
           scale = TRUE,
           n = 4,
           loading = FALSE,
           n_loadings = 10,
           gene_id = "SYMBOL") {
    standardGeneric("pca_matrix_plot")
  }
)

#' Boxplot of median expression or scores of signature
#'
#' Make boxplot and show expression or score level of signature across subsets.
#'
#' @param data expression data, can be matrix, DGEList, eSet, seurat, sce...
#' @param sigs a vector of signature (Symbols)
#' @param group_col character or vector, specify the column name to compare in coldata
#' @param target_group pattern, specify the group of interest as reference
#' @param type one of "score" and "expression", to plot score or expression of
#'             the signature
#' @param method a character string indicating which method to be used for
#'               `stat_compare_means()` to compare the means across groups,
#'               could be "t.test", 'wilcox.test', 'anova'..., default "t.test"
#' @param slot character, indicate which slot used as expression, optional
#' @param gene_id character, indicate the ID type of rowname of expression data's ,
#'                could be one of 'ENSEMBL', 'SYMBOL', ... default 'SYMBOL'
#'
#' @return patchwork or ggplot of boxplot
#'
#' @examples
#' data("im_data_6", "nk_markers")
#' p <- sig_boxplot(
#'   im_data_6,
#'   sigs = nk_markers$HGNC_Symbol[1:30],
#'   group_col = "celltype:ch1", target_group = "NK",
#'   gene_id = "ENSEMBL"
#' )
#'
#' @include plot.R
#' @import ggplot2
#' @export
setGeneric(
  "sig_boxplot",
  function(data,
           sigs,
           group_col,
           target_group,
           type = c("score", "expression"),
           method = "t.test",
           slot = "counts",
           gene_id = "SYMBOL") {
    standardGeneric("sig_boxplot")
  }
)

#' Visualize GSEA result with input list of gene symbols.
#'
#' Visualize GSEA result with multiple lists of genes by using `clusterProfiler`.
#'
#' @inheritParams sig_boxplot
#' @param sigs a vector of signature (Symbols) or a list of signatures
#' @param method one of "gseaplot" and "dotplot", how to plot GSEA result
#' @param col column name of [clusterProfiler::GSEA()] result, used for dot
#'            col when method = "dotplot"
#' @param size column name of [clusterProfiler::GSEA()] result, used for dot
#'             size when method = "dotplot"
#' @param pvalue_table logical, if to add p value table if method = "gseaplot"
#' @param digits num, specify the number of significant digits of pvalue table
#' @param rank_stat character, specify which metric used to rank for GSEA,
#'                  default "logFC"
#' @param ... params for function [get_de_table()] and function [enrichplot::gseaplot2()]
#'
#' @return patchwork object for all comparisons
#'
#' @examples
#' data("im_data_6", "nk_markers")
#' sig_gseaplot(
#'   sigs = list(
#'     A = nk_markers$HGNC_Symbol[1:15],
#'     B = nk_markers$HGNC_Symbol[20:40],
#'     C = nk_markers$HGNC_Symbol[60:75]
#'   ),
#'   data = im_data_6, group_col = "celltype:ch1",
#'   target_group = "NK", gene_id = "ENSEMBL"
#' )
#'
#' @include plot.R
#' @export
setGeneric(
  "sig_gseaplot",
  function(data,
           sigs,
           group_col,
           target_group,
           gene_id = "SYMBOL",
           slot = "counts",
           method = c("dotplot", "gseaplot"),
           col = "-log10(p.adjust)",
           size = "enrichmentScore",
           pvalue_table = FALSE,
           digits = 2,
           rank_stat = "logFC",
           ...) {
    standardGeneric("sig_gseaplot")
  }
)

#' Heatmap original markers and screened signature
#'
#' Compare the heatmap before and after screening.
#'
#' @inheritParams sig_gseaplot
#' @inheritParams sig_boxplot
#' @param markers a vector of gene names, listed the gene symbols of original
#'                markers pool
#' @param scale could be one of 'none' (default), 'row' or 'column'
#' @param ranks_plot logical, if to use ranks instead of expression of genes to
#'                   draw heatmap
#' @param ... params for [ComplexHeatmap::Heatmap()]
#'
#' @return patchwork object of heatmap
#'
#' @examples
#' data("im_data_6", "nk_markers")
#' sig_heatmap(
#'   data = im_data_6, sigs = nk_markers$HGNC_Symbol[1:10],
#'   group_col = "celltype:ch1",
#'   gene_id = "ENSEMBL"
#' )
#'
#' @include plot.R
#' @import ggplot2 patchwork
#' @export
setGeneric(
  "sig_heatmap",
  function(data,
           sigs,
           group_col,
           markers,
           scale = c("none", "row", "column"),
           gene_id = "SYMBOL",
           ranks_plot = FALSE,
           slot = "counts",
           ...) {
    standardGeneric("sig_heatmap")
  }
)

#' Plot rank density
#'
#' Show the rank density of given signature in the given comparison.
#'
#' @inheritParams sig_boxplot
#' @param aggregate logical, if to aggregate expression according to `group_col`,
#'                  default FALSE
#'
#' @return ggplot or patchwork
#'
#' @examples
#' data("im_data_6", "nk_markers")
#' sig_rankdensity_plot(
#'   data = im_data_6, sigs = nk_markers$HGNC_Symbol[1:10],
#'   group_col = "celltype:ch1", gene_id = "ENSEMBL"
#' )
#'
#' @export
setGeneric(
  "sig_rankdensity_plot",
  function(data,
           sigs,
           group_col,
           aggregate = FALSE,
           slot = "counts",
           gene_id = "SYMBOL") {
    standardGeneric("sig_rankdensity_plot")
  }
)

#' @title Scatter plot of signature for specific subset vs others
#'
#' @description  Scatter plot depicts mean expression for each signature gene in
#'   the specific subset against other cell types.
#'
#' @inheritParams sig_boxplot
#' @param xint intercept of vertical dashed line, default 1
#' @param yint intercept of horizontal dashed line, default 1
#'
#' @return patchwork or ggplot of scatter plot of median expression
#'
#' @examples
#' data("im_data_6", "nk_markers")
#' sig_scatter_plot(
#'   sigs = nk_markers$HGNC_Symbol, data = im_data_6,
#'   group_col = "celltype:ch1", target_group = "NK",
#'   gene_id = "ENSEMBL"
#' )
#'
#' @include plot.R
#' @import ggplot2
#' @export
setGeneric(
  "sig_scatter_plot",
  function(data,
           sigs,
           group_col,
           target_group,
           slot = "counts",
           xint = 1,
           yint = 1,
           gene_id = "SYMBOL") {
    standardGeneric("sig_scatter_plot")
  }
)
