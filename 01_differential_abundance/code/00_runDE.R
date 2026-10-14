library(arrow)
library(limma)
library(edgeR)
library(broom)
library(tidyverse)
source("functions.R")

# DATA IMPORT -------------------------------------------------------------

# meta
meta <- readRDS("00_setup/data/metadata.samples.postQC.rds")
analytes <- readRDS("00_setup/data/metadata.analytes.rds")
analytes <- analytes[, c(
  "apt_id",
  "TargetFullName",
  "Target",
  "UniProt",
  "ensembl_id",
  "EntrezGeneID",
  "EntrezGeneSymbol",
  "Dilution",
  "SeqId",
  "apt_label"
)]

# counts
tidy_counts <- read_parquet("00_setup/data/count_tidy.postQC.parquet")
counts <- tidy_counts %>%
  pivot_wider(
    id_cols = apt_id,
    names_from = study_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()
stopifnot(identical(colnames(counts), meta$study_id))

# RUN LIMMA ---------------------------------------------------------------
# reference diagnosis
design <- model.matrix(
  ~ 0 +
    ref_dx +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = meta
)
colnames(design) <- make.names(colnames(design))
contrast_matrix <- makeContrasts(
  cGCA = ref_dxcGCA - ref_dxnot_cGCA,
  levels = colnames(design)
)
ref_dx <- runLimma(
  metadata = meta,
  counts = counts,
  design = design,
  contrastMatrix = contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id",
  outputCsvPath = "01_differential_abundance/data/limma.ref_dx.csv"
)
ref_dx$results %>% filter(adj.P.Val < 0.05)


# CRP
crp_cor <- tidy_counts |>
  left_join(select(meta, study_id, baseline_crp), by = "study_id") |>
  group_by(apt_id) |>
  do(tidy(cor.test(
    ~ log2_RFU + baseline_crp,
    data = .,
    method = "spearman"
  ))) |>
  ungroup(apt_id)
saveRDS(crp_cor, "01_differential_abundance/data/soma_baseline_crp_spearman.rds")

# ESR
esr_cor <- tidy_counts |>
  left_join(select(meta, study_id, baseline_esr), by = "study_id") |>
  group_by(apt_id) |>
  do(tidy(cor.test(
    ~ log2_RFU + baseline_esr,
    data = .,
    method = "spearman"
  ))) |>
  ungroup(apt_id)
saveRDS(esr_cor, "01_differential_abundance/data/soma_baseline_esr_spearman.rds")
# platelets
platelet_cor <- tidy_counts |>
  left_join(select(meta, study_id, platelets), by = "study_id") |>
  group_by(apt_id) |>
  do(tidy(cor.test(~ log2_RFU + platelets, data = ., method = "spearman"))) |>
  ungroup(apt_id)
saveRDS(platelet_cor, "01_differential_abundance/data/soma_platelets_spearman.rds")
