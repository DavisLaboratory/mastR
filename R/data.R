#' LM7 matrix for CIBERSORT.
#'
#' A dataset containing 375 marker genes expression of 7 immune subsets which
#' is generated for CIBERSORT.
#'
#' @format A data frame with 375 rows 9 variables:
#' \describe{
#'   \item{Gene}{gene symbols}
#'   \item{B cells}{gene median expression in B cells}
#'   ...
#'   \item{Subset}{immune subset of the marker gene}
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
#'   ...
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
#'   ...
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

