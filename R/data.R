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
#'   \item{T γδ}{gene median expression in T γδ cells}
#'   \item{NK}{gene median expression in NK cells}
#'   \item{MoMaDC}{gene median expression in MoMaDC cells}
#'   \item{granulocytes}{gene median expression in granulocytes}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5384348/}
"LM7"

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
#' @source \url{https://cibersort.stanford.edu/}
"LM22"

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
#' @source \url{https://cancerimmunolres.aacrjournals.org/content/7/7/1162.long}
"NK_markers"

#' RNA-seq TMM normalized counts data of 6 sorted immune subsets.
#'
#' A list of ExpressionSet objects containing 6 immune subsets (B-cells, CD4,
#' CD8, Monocytes, Neutrophils, NK) from healthy individuals.
#'
#' @format A list of ExpressionSet objects of 6*4 samples:
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424}
"im_data_6"

#' RNA-seq TPM data of 5 CRC cell line samples from CCLE.
#'
#' A test DGEList object with RNA-seq RSEM quantified TPM data of 5 CRC cell
#' line samples from CCLE [depmap::depmap_TPM()].
#'
#' @format A DGEList of `r nrow(ccle_crc_5)` genes * 5 samples:
#' @source [depmap::depmap_TPM()]
"ccle_crc_5"

