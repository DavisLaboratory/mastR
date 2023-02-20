## code to prepare `ccle_crc_5` dataset goes here

### download CCLE TPM data from depmap
CCLE <- depmap::depmap_TPM()
### convert long data into wide data
CCLE <- CCLE[,c("gene_name", "cell_line", "rna_expression")] |>
  tidyr::pivot_wider(names_from = "cell_line",
                     values_from = "rna_expression") |>
  data.frame()
rownames(CCLE) <- CCLE$gene_name
CCLE$gene_name <- NULL

### create DGEList with 5 CRC cell line samples
ccle_crc_5 <- edgeR::DGEList(
  counts = CCLE[, grep("LARGE_INTESTINE", colnames(CCLE))[1:5]]
)
ccle_crc_5$samples$cancer <- "CRC"

### save data
usethis::use_data(ccle_crc_5, overwrite = TRUE)
