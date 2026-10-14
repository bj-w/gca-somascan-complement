library(tidyverse)
library(scuttle)
library(edgeR)
library(sparseMatrixStats)
library(DelayedArray)
library(arrow)
library(anndataR)

# function to format scRNAseq data and produce nTPM matrix ----------------
scrna_nTPM <- function(sce) {
  # use deontXcounts - ambient RNA corrected counts
  # setting this as the count assay which is used in downstream SCE functions
  names(assays(sce))[which(names(assays(sce)) == "decontXcounts")] <- "counts"
  # keep protein coding genes
  sce <- sce[rowData(sce)$feature_type == "protein_coding", ]
  message(paste0(nrow(colData(sce)), " cells"))
  ### pTPM
  message("Calculating total counts per gene in each cell type")
  gene_counts_per_cell_type <- assay(summarizeAssayByGroup(
    sce,
    ids = colData(sce)$cell_type,
    statistics = "sum"
  ))
  message("Calculating total counts per cell type")
  total_counts_per_cell_type <- colSums(gene_counts_per_cell_type)
  gene_counts_per_cell_type <- gene_counts_per_cell_type[, match(
    names(total_counts_per_cell_type),
    colnames(gene_counts_per_cell_type)
  )]
  stopifnot(identical(
    colnames(gene_counts_per_cell_type),
    names(total_counts_per_cell_type)
  ))
  # pTPM matrix
  message("Generating pTPM matrix")
  pTPM <- sweep(
    gene_counts_per_cell_type,
    MARGIN = 2,
    STATS = total_counts_per_cell_type,
    FUN = "/"
  )
  pTPM <- pTPM * 1e+6
  ### TMM
  # use in-built function from edgeR to extract per-cell normalisation factors
  tmm <- DGEList(counts = pTPM)
  message("Calculating TMM normalisation factors")
  tmm <- calcNormFactors(tmm, method = "TMM")
  # divide values in each column by the normalisation factor
  message("Generating nTPM matrix")
  tmm <- sweep(tmm$counts, 2, tmm$samples$norm.factors, "/")
  return(tmm)
}


# vasculature -------------------------------------------------------------
# https://datasets.cellxgene.cziscience.com/c542113b-cdb7-4a12-9c91-a7363256b872.h5ad

vasc <- read_h5ad(
  "raw_data/tissue_enrichment/TSv2_vasculature_cellxgene_20250826.h5ad",
  as = "SingleCellExperiment"
)

vasc <- scrna_nTPM(sce = vasc)
# 42650 cells
vasc %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  write_parquet("01_differential_abundance/data/vasculature_nTPM.parquet")
gc()

# blood -------------------------------------------------------------------
# https://datasets.cellxgene.cziscience.com/988defcd-7e39-4d07-91b9-a9853af1e769.h5ad
blood <- read_h5ad(
  "raw_data/tissue_enrichment/TSv2_blood_cellxgene_20250904.h5ad",
  as = "SingleCellExperiment"
)
blood <- scrna_nTPM(sce = blood)
# 85233 cells
blood %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  write_parquet("01_differential_abundance/data/blood_nTPM.parquet")
gc()
