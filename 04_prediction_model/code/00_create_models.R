library(nestedcv)
library(arrow)
library(future)
library(tidyverse)
source("functions.R")
theme_set(figTheme())
# import ------------------------------------------------------------------
meta <- readRDS("00_setup/data/metadata.samples.postQC.rds")
counts <- read_parquet("00_setup/data/count_tidy.postQC.parquet")

# format counts -----------------------------------------------------------
cm <- counts %>%
  pivot_wider(
    id_cols = "apt_id",
    names_from = study_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()
stopifnot(identical(colnames(cm), meta$study_id))
# transpose
cm <- t(cm)

# model building ----------------------------------------------------------
response <- as.factor(ifelse(meta$ref_dx == "cGCA", 1, 0))

# leave one out cv
set.seed(1997, "L'Ecuyer-CMRG")
res <- nestcv.glmnet(
  x = cm,
  y = response,
  n_inner_folds = 5,
  outer_method = "LOOCV",
  family = "binomial",
  alphaSet = 1, # LASSO
  min_1se = 1, # set lamba to 1se
  finalCV = TRUE
)
summary(res)
# export
saveRDS(res, "04_prediction_model/data/nestedcv.rds")
