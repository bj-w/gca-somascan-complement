library(arrow)
library(limma)
library(edgeR)
library(tidyverse)
source("functions.R")
# DATA IMPORT -------------------------------------------------------------
# meta
meta <- readRDS("00_setup/data/metadata.samples.postQC.rds")
# factorise cranial ischaemia complications column
meta$cranial_isch_comp_rm <- factor(meta$cranial_isch_comp_rm)
analytes <- readRDS("00_setup/data/metadata.analytes.rds")
analytes <- analytes[, c(
  "apt_id",
  "TargetFullName",
  "Target",
  "UniProt",
  "EntrezGeneID",
  "EntrezGeneSymbol",
  "Dilution",
  "SeqId"
)]

# counts
counts <- read_parquet("00_setup/data/count_tidy.postQC.parquet") %>%
  pivot_wider(
    id_cols = apt_id,
    names_from = study_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()


# filter for cGCA only
meta <- meta %>% filter(ref_dx == "cGCA")
counts <- counts[, colnames(counts) %in% meta$study_id]


stopifnot(identical(colnames(counts), meta$study_id))

# import histology
histology_meta <- readRDS("00_setup/data/meta.histology.rds") %>%
  as.data.frame()
# restrict to GCA cases
length(unique(histology_meta$study_id))
histology_meta <- histology_meta %>% filter(study_id %in% meta$study_id)
length(unique(histology_meta$study_id))
# exclude samples with no slides
histology_meta <- histology_meta %>% filter(is.na(comment))
length(unique(histology_meta$study_id))
# add on meta info to histolgoy
histology_meta <- unique(histology_meta[, c(
  "study_id",
  "max_barcelona_score_recode",
  "max_media_pattern_recode",
  "max_media_destruction_recode",
  "max_occlusion_grade_recode",
  "max_barcelona_score"
)]) %>%
  inner_join(meta, by = "study_id")
# create separate counts for histology
histology_counts <- counts[, colnames(counts) %in% histology_meta$study_id]
histology_counts <- histology_counts[, match(
  histology_meta$study_id,
  colnames(histology_counts)
)]
stopifnot(identical(colnames(histology_counts), histology_meta$study_id))
row.names(histology_meta) <- histology_meta$study_id

# export sample N
saveRDS(histology_meta, "01_differential_abundance/data/gca_only_meta.rds")

de <- list()


# cranial ischaemia complications -----------------------------------------
meta$cranial_isch_comp_rm %>% table()
de$cranial_isch <- list()
de$cranial_isch$design <- model.matrix(
  ~ 0 +
    cranial_isch_comp_rm +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = meta
)
colnames(de$cranial_isch$design) <- make.names(colnames(de$cranial_isch$design))
de$cranial_isch$contrast_matrix <- makeContrasts(
  isch_comp = cranial_isch_comp_rmYes - cranial_isch_comp_rmNo,
  levels = colnames(de$cranial_isch$design)
)


de$cranial_isch$limma <- runLimma(
  metadata = meta,
  counts = counts,
  design = de$cranial_isch$design,
  contrastMatrix = de$cranial_isch$contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id"
)

# polymyalgia rheumatica symptoms -----------------------------------------
table(meta$pmryesno)
de$pmr <- list()
de$pmr$design <- model.matrix(
  ~ 0 +
    pmryesno +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = meta
)
colnames(de$pmr$design) <- make.names(colnames(de$pmr$design))
de$pmr$contrast_matrix <- makeContrasts(
  pmr = pmryesnoyes - pmryesnono,
  levels = colnames(de$pmr$design)
)
de$pmr$limma <- runLimma(
  meta = meta,
  counts = counts,
  design = de$pmr$design,
  contrastMatrix = de$pmr$contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id",
)


# barcelona score ---------------------------------------------------------

table(histology_meta$max_barcelona_score_recode)
de$barca_score <- list()
de$barca_score$design <- model.matrix(
  ~ 0 +
    max_barcelona_score_recode +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = histology_meta
)
colnames(de$barca_score$design) <- make.names(colnames(de$barca_score$design))
de$barca_score$contrast_matrix <- makeContrasts(
  max_barcelona_score_recodeyes - max_barcelona_score_recodeno,
  levels = colnames(de$barca_score$design)
)

de$barca_score$limma <- runLimma(
  metadata = histology_meta,
  counts = histology_counts,
  design = de$barca_score$design,
  contrastMatrix = de$barca_score$contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id"
)

# media destruction -------------------------------------------------------
de$media_destruction <- list()
de$media_destruction$design <- model.matrix(
  ~ 0 +
    max_media_destruction_recode +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = histology_meta
)
colnames(de$media_destruction$design) <- make.names(colnames(
  de$media_destruction$design
))
de$media_destruction$contrast_matrix <- makeContrasts(
  media_destruction = max_media_destruction_recodeyes -
    max_media_destruction_recodeno,
  levels = colnames(de$media_destruction$design)
)
de$media_destruction$limma <- runLimma(
  metadata = histology_meta,
  counts = histology_counts,
  design = de$media_destruction$design,
  contrastMatrix = de$media_destruction$contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id"
)


# media pattern -----------------------------------------------------------
de$media_pattern <- list()
de$media_pattern$design <- model.matrix(
  ~ 0 +
    max_media_pattern_recode +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = histology_meta
)
colnames(de$media_pattern$design) <- make.names(colnames(
  de$media_pattern$design
))
de$media_pattern$contrast_matrix <- makeContrasts(
  media_pattern = max_media_pattern_recodemultifocal_diffuse -
    max_media_pattern_recodenormal_focal,
  levels = colnames(de$media_pattern$design)
)
de$media_pattern$limma <- runLimma(
  metadata = histology_meta,
  counts = histology_counts,
  design = de$media_pattern$design,
  contrastMatrix = de$media_pattern$contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id"
)
de$media_pattern$limma$results


# occlusion grade ---------------------------------------------------------
de$occlusion_grade <- list()
de$occlusion_grade$design <- model.matrix(
  ~ 0 +
    max_occlusion_grade_recode +
    poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = histology_meta
)
colnames(de$occlusion_grade$design) <- make.names(colnames(
  de$occlusion_grade$design
))
de$occlusion_grade$contrast_matrix <- makeContrasts(
  occlusion_grade = max_occlusion_grade_recodehigh -
    max_occlusion_grade_recodelow,
  levels = colnames(de$occlusion_grade$design)
)
de$occlusion_grade$limma <- runLimma(
  metadata = histology_meta,
  counts = histology_counts,
  design = de$occlusion_grade$design,
  contrastMatrix = de$occlusion_grade$contrast_matrix,
  featureInfo = analytes,
  featureId = "apt_id"
)
de$occlusion_grade$limma$results


# collect all de results objects ------------------------------------------
combined_de <- map(names(de), \(x) pluck(de, x, "limma", "results")) %>%
  set_names(names(de))

saveRDS(combined_de, "01_differential_abundance/data/gca_only_limma_output.rds")
