#' @import ggplot2
#' @import methods
#' @import stats
#' @import utils
#' @import graphics
#' @import patchwork
#' @import ggfortify
#' @import ggpubr
#' @import org.Hs.eg.db
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr group_by summarise_all arrange_
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom Biobase exprs pData
#' @importFrom depmap depmap_TPM depmap_metadata
NULL

#' Screen Immune Cells Signature for Specific Cancer or Tissue Type
#'
#' `TILsignaturesScreen` This package enables screening immune cells signatures
#'   for specific cancer or tissue type. Imported NK cell markers are based on
#'   LM22 and LM7 and human orthologs in mice. MSigDB genesets of relevant terms
#'   and Panglao DB markers can also be retrieved by this package. Users can
#'   choose to use sorted immune cell bulk RNA-seq data and CCLE adherent cell
#'   line RNA-seq data (CCLE) to screen signatures or use customized data.
#'   This package also provides pseudo bulking function to deal with scRNA-seq data. Multiple
#'   visualization functions are implemented in this package.
#'
#' @author Jinjin Chen \email{chen.j@@wehi.edu.au}
#' @name TILsignaturesScreen_Package
#' @docType package
#' @aliases TILsignaturesScreen TILsignaturesScreen_Package
#' @keywords internal
#'
NULL
