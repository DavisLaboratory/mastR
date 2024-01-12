#' @import ggplot2
#' @import methods
#' @import stats
#' @import graphics
#' @import patchwork
#' @import limma
#' @import edgeR
#' @import ggpubr
#' @import org.Hs.eg.db
#' @import msigdb
#' @import GSEABase
#' @importFrom utils head stack globalVariables modifyList
#' @importFrom SeuratObject GetAssayData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay colData
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr group_by summarise_all summarise arrange inner_join
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom Biobase exprs pData rowMedians
NULL

#' Screen Immune Cells Signature for Specific Cancer or Tissue Type
#'
#' `mastR` This package enables automated screening of group specific signature
#'   for specific tissues. The package is developed for generating refined lists
#'   of signature genes from multiple group comparisons based on the results
#'   from edgeR and limma differential expression (DE) analysis workflow. It
#'   also takes into account the background expression of tissue-specificity,
#'   which is often ignored by other markers generation tools. This package also
#'   provides pseudo bulking function to deal with scRNA-seq data. Multiple
#'   visualization functions are implemented in this package.
#'
#' @author Jinjin Chen \email{chen.j@@wehi.edu.au}
#' @name mastR_Package
#' @docType package
#' @aliases mastR mastR_Package
#' @keywords internal
#' @return Automated screened signature
#'
NULL

#' LM7 matrix for CIBERSORT.
#'
#' A dataset containing 375 marker genes expression of 7 immune subsets which
#' is generated for CIBERSORT.
#'
#' @format A data frame with 375 rows 9 variables:
#' \describe{
#'   \item{Gene}{gene symbols}
#'   \item{Subset}{immune subset of the marker gene}
#'   \item{B cells}{gene median expression in B cells}
#'   \item{T CD4}{gene median expression in T CD4 cells}
#'   \item{T CD8}{gene median expression in T CD8 cells}
#'   \item{T gamma delta}{gene median expression in T gamma delta cells}
#'   \item{NK}{gene median expression in NK cells}
#'   \item{MoMaDC}{gene median expression in MoMaDC cells}
#'   \item{granulocytes}{gene median expression in granulocytes}
#' }
#' @usage data(lm7)
#' @return data frame
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5384348/}
"lm7"

#' LM22 matrix for CIBERSORT.
#'
#' A dataset containing 547 marker genes expression of 22 immune subsets which
#' is generated for CIBERSORT.
#'
#' @format A data frame with 547 rows 23 variables:
#' \describe{
#'   \item{Gene}{gene symbols}
#'   \item{B cells naive}{0 or 1, represents if the gene is significantly
#'         up-regulated in the subset}
#'   \item{B cells memory}{0 or 1}
#'   \item{Plasma cells}{0 or 1}
#'   \item{T cells CD8}{0 or 1}
#'   \item{T cells CD4 naive}{0 or 1}
#'   \item{T cells CD4 memory resting}{0 or 1}
#'   \item{T cells CD4 memory activated}{0 or 1}
#'   \item{T cells follicular helper}{0 or 1}
#'   \item{T cells regulatory (Tregs)}{0 or 1}
#'   \item{T cells gamma delta}{0 or 1}
#'   \item{NK cells resting}{0 or 1}
#'   \item{NK cells activated}{0 or 1}
#'   \item{Monocytes}{0 or 1}
#'   \item{Macrophages M0}{0 or 1}
#'   \item{Macrophages M1}{0 or 1}
#'   \item{Macrophages M2}{0 or 1}
#'   \item{Dendritic cells resting}{0 or 1}
#'   \item{Dendritic cells activated}{0 or 1}
#'   \item{Mast cells resting}{0 or 1}
#'   \item{Mast cells activated}{0 or 1}
#'   \item{Eosinophils}{0 or 1}
#'   \item{Neutrophils}{0 or 1}
#' }
#' @usage data(lm22)
#' @return data frame
#' @source \url{https://cibersort.stanford.edu/}
"lm22"

#' NK cell markers combination.
#'
#' A dataset containing 114 NK cell markers from LM22,
#' LM7 and human orthologs in mice.
#'
#' @format A data frame with 114 rows and  at least 4 variables:
#' \describe{
#'   \item{HGNC_Symbol}{gene symbols}
#'   \item{LM22}{if included in LM22}
#'   \item{LM7}{if included in LM7}
#'   \item{Huntington}{if included in orthologs}
#' }
#' @usage data(nk_markers)
#' @return data frame
#' @source \url{https://cancerimmunolres.aacrjournals.org/content/7/7/1162.long}
"nk_markers"

#' RNA-seq TMM normalized counts data of 6 sorted immune subsets.
#'
#' An ExpressionSet objects containing 6 immune subsets (B-cells, CD4,
#' CD8, Monocytes, Neutrophils, NK) from healthy individuals.
#'
#' @format An ExpressionSet objects of 6*4 samples.
#' @usage data(im_data_6)
#' @return ExpressionSet
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424}
"im_data_6"

#' RNA-seq TPM data of 5 CRC cell line samples from CCLE.
#'
#' A test DGEList object with RNA-seq RSEM quantified TPM data of 5 CRC cell
#' line samples from CCLE [depmap::depmap_TPM()].
#'
#' @format A DGEList of `r nrow(mastR::ccle_crc_5)` genes * 5 samples.
#' @usage data(ccle_crc_5)
#' @return DGEList
#' @source [depmap::depmap_TPM()]
"ccle_crc_5"

#' Sub-collection of MSigDB gene sets.
#'
#' A small GeneSetCollection object, contains gene sets with gene set name matched
#' to 'NATURAL_KILLER' from GO:BP MSigDB v7.4 database.
#'
#' @format A GeneSetCollection of `r length(mastR::msigdb_gobp_nk)` gene sets.
#' @usage data(msigdb_gobp_nk)
#' @return GeneSetCollection
#' @source [msigdb::getMsigdb()]
"msigdb_gobp_nk"
