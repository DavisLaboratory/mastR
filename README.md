# mastR: Markers Automated Screening Tool in R

[![R-CMD-check](https://github.com/DavisLaboratory/mastR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/DavisLaboratory/mastR/actions)
[![BioC status](https://bioconductor.org/shields/years-in-bioc/mastR.svg)](https://bioconductor.org/packages/mastR/)

An R package for automatically screening group specific signature for specific tissues.

mastR is an R package designed for automated screening of signatures of interest for specific research questions. The package is developed for generating refined lists of signature genes from multiple group comparisons based on the results from edgeR and limma differential expression (DE) analysis workflow.

It also takes into account the background expression of tissue-specificity, which is often ignored by other markers generation tools. The package allows users to input their expression data containing the groups of interest, along with group labels, to obtain a list of marker genes (signature) for the target group, and remove the genes with background expression from the signature.

This package is particularly useful for the identification of group markers in various biological and medical applications, including cancer research and developmental biology.

[**Check out the standard demonstration.**](https://davislaboratory.github.io/mastR/articles/mastR_Demo.html)

## Installation

mastR can be installed from Bioconductor directly as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mastR")
```
