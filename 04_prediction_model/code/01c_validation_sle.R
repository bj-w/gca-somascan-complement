library(readxl)
library(tidyverse)
library(data.table)
library(janitor)
library(arrow)
library(nestedcv)


rfu <- read_delim("raw_data/sle/soma-sle-mat.txt", delim = "\t")
analytes <- readRDS("00_setup/data/metadata.analytes.rds")
# metadata
meta <- read_delim("raw_data/sle/soma-sle-pheno.txt")
# just batch 2 + required columns
meta <- meta %>% filter(Batch == 2) %>% select(SampleId, Group)


# RFU
rfu <- rfu %>%
  pivot_longer(cols = !SampleId, names_to = "AptName", values_to = "RFU") %>%
  # just batch 2
  filter(SampleId %in% meta$SampleId) %>%
  left_join(analytes[, c("AptName", "apt_id")], by = "AptName") %>%
  # drop SOMAmers, with no apt_id - will be QC.
  filter(!is.na(apt_id)) %>%
  mutate(log2_RFU = log2(RFU)) %>%
  pivot_wider(
    id_cols = apt_id,
    names_from = SampleId,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()

rfu <- rfu[, match(meta$SampleId, colnames(rfu))]
stopifnot(identical(colnames(rfu), meta$SampleId))

# model
model <- readRDS("04_prediction_model/data/nestedcv.rds")

# transpose rfu matrix
rfu <- t(rfu)
# select model predictors
rfu <- rfu[, colnames(rfu) %in% row.names(model$final_coef)]
# run prediction
prediction <- data.frame(
  prediction = predict(model, newdata = rfu, type = 'response')
) %>%
  rename(prediction = 1) %>%
  rownames_to_column("SampleId")
# add on group info
prediction <- prediction %>% left_join(meta, by = "SampleId")
write_parquet(prediction, "04_prediction_model/data/sle.prediction.parquet")
