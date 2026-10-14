library(here)
library(SomaDataIO)
library(variancePartition)
library(arrow)
library(tidyverse)
library(ggtext)
require(ggrepel)
library(broom)
library(patchwork)


source("functions.R")
theme_set(figTheme())
dir.create("00_setup/data/qc/", showWarnings = FALSE)

# import ------------------------------------------------------------------
meta <- as.data.frame(readRDS(here("00_setup/data/metadata.samples.preQC.rds")))
row.names(meta) <- meta$study_id # stops variance partition warning
analytes <- readRDS(here("00_setup/data/metadata.analytes.rds"))
counts <- read_parquet(here("00_setup/data/count_tidy.preQC.parquet"))
adat <- read_adat(
  "raw_data/ICL-SS-2343676_2025-02-21/SS-2343676_v4.1_CitratePlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
)


# exclude EP-01-104 due to leak
meta <- meta %>% filter(!study_id == "EP-01-104")
counts <- counts |> filter(!study_id == "EP-01-104")
adat <- adat %>% filter(!SampleId == "EP-01-104")

# outlier detection -------------------------------------------------------
# first method uses PCA
# calculating PCs
prcomp_obj_pre_qc <- counts %>%
  pivot_wider(
    id_cols = study_id,
    names_from = apt_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("study_id") %>%
  prcomp(center = TRUE, scale. = TRUE)
# export these PCs in case I need to plot them later
saveRDS(prcomp_obj_pre_qc, "00_setup/data/qc/prcompObject_preQC.rds")

# run PCA outlier detection. Will flag samples +/- 5 SDs from the mean PC1 and
# PC2 values.
outlier_pca <- outlierPCA(prcompObject = prcomp_obj_pre_qc, sdLimit = 5)
outlier_pca$plot
saveRDS(outlier_pca, "00_setup/data/qc/outlier.pca.rds")


# second method uses abundance median and IQR. Will flag samples +/- 5 SDs from
# the abundance median and IQR.
outlier_median_iqr <- outlierMedianIQR(
  counts_df = counts,
  fillGroup = outlier,
  sdLimit = 5,
  metadata = meta
)
outlier_median_iqr$plot
saveRDS(outlier_median_iqr, "00_setup/data/qc/outlier.median_iqr.rds")


# defining outlying samples as those flagged by the two methods above.
outliers_to_remove <- union(
  outlier_pca$outliers,
  outlier_median_iqr$outliers
)

# remove the outliers from counts and metadata, then export.
counts <- counts %>% filter(!study_id %in% outliers_to_remove)
write_parquet(counts, "00_setup/data/count_tidy.postQC.parquet")

meta <- meta %>% filter(!study_id %in% outliers_to_remove)
saveRDS(meta, "00_setup/data/metadata.samples.postQC.rds")

# also export as count matrix
mat <- counts %>%
  pivot_wider(id_cols = apt_id, names_from = study_id, values_from = log2_RFU)

write_parquet(mat, "00_setup/data/count_matrix.log2.postQC.parquet")

# export residualised count matrix
# lod ---------------------------------------------------------------------
# using buffer samples to estimate the limit of detection.
elod <- adat %>%
  filter(SampleType == "Buffer") %>%
  calc_eLOD()
# adding on SOMAmer info
elod <- counts %>%
  left_join(analytes[, c("apt_id", "AptName")], by = "apt_id") %>%
  left_join(elod, by = c("AptName" = "SeqId")) %>%
  select(apt_id, study_id, RFU, eLOD)
# assessing the proportion of SOMAmers above the limit of detection in each sample.
elod <- elod %>%
  mutate(eLOD_threshold = as.integer(RFU > eLOD)) %>%
  reframe(
    n_total = n_distinct(apt_id),
    n_above = sum(eLOD_threshold),
    prop = n_above / n_total,
    .by = study_id
  )
summary(elod$prop)
write_parquet(elod, "00_setup/data/qc/elod.parqet")


# PCA ---------------------------------------------------------------------
# recalculating PCs after outlier removal.
prcomp_obj <- counts %>%
  pivot_wider(
    id_cols = study_id,
    names_from = apt_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("study_id") %>%
  prcomp(center = TRUE, scale. = TRUE)
saveRDS(prcomp_obj, "00_setup/data/qc/prcompObject_postQC.rds")

# exemplar PCA biplots.
plotPCABiplot(
  prcompObject = prcomp_obj,
  metadata = meta,
  color_var = "row_check",
  x_pc = "PC1",
  y_pc = "PC2",
  sample_id_col = "study_id",
  marginalDensity = TRUE
)
plotPCABiplot(
  prcompObject = prcomp_obj,
  metadata = meta,
  color_var = "ref_dx",
  x_pc = "PC1",
  y_pc = "PC2",
  sample_id_col = "study_id",
  marginalDensity = TRUE
)

# PCA association ---------------------------------------------------------
# measures the association between potential covariates and abundance PCs.

# create dataframe containing the potential sources of variation
variation_metrics <- meta %>%
  select(
    study_id,
    Age = age_at_dx,
    Sex = sex,
    Diagnosis = ref_dx,
    `SomaScan plate` = plate,
    `Glucocorticoid treatment duration` = gc_days,
    `Sample storage time` = freezer_months,
    `Median RFU` = median_RFU,
    `Median log2 RFU` = median_log2_RFU,
    `Total RFU` = total_RFU,
    `Total log2 RFU` = total_log2_RFU,
    `IQR RFU` = iqr_RFU,
    `IQR log2 RFU` = iqr_log2_RFU
  )

# function to perform the association
pca_assoc <- lmPCA(
  metricDF = variation_metrics,
  prcompObject = prcomp_obj,
  uniqueIdColumn = "study_id",
  PC = 10,
  title = NULL
)
saveRDS(pca_assoc, "00_setup/data/qc/pca_covariate_associations.rds")

# variance partition ------------------------------------------------------
# specify linear mixed model formula
lmm_formula <- ~ (1 | ref_dx) +
  (1 | sex) +
  age_at_dx +
  iqr_log2_RFU +
  gc_days
# log2 count/rfu matrix as input (rows = probes, columns = samples)
log2_rfu_mat <- counts %>%
  pivot_wider(
    id_cols = apt_id,
    names_from = "study_id",
    values_from = "log2_RFU"
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()
# ensure samples in count matrix and metadata are in the same order.
# will error if not.
stopifnot(identical(colnames(log2_rfu_mat), meta$study_id))
# generate linear mixed models
var_part <- fitExtractVarPartModel(
  exprObj = log2_rfu_mat,
  formula = lmm_formula,
  data = meta
)
# export
saveRDS(var_part, "00_setup/data/qc/variancePartition.rds")
# plot
plotVarPart(sortCols(var_part), label.angle = 60) +
  figTheme() +
  theme(legend.position = "none")
