library(arrow)
library(nestedcv)
library(tidyverse)

#source("functions.R")
model <- readRDS("04_prediction_model/data/nestedcv.rds")

# GCA (Cunningham) --------------------------------------------------------
cunningham <- list()
cunningham$meta <- readRDS("04_prediction_model/data/cunningham.meta.rds")
cunningham$rfu <- read_parquet(
  "04_prediction_model/data/cunningham.tidy_counts.parquet"
)

# drop those without apt_id (not on cunning somascan panel?)
cunningham$rfu <- cunningham$rfu %>% drop_na(apt_id)

# convert to log2 rfu count matrix
cunningham$cm <- cunningham$rfu %>%
  pivot_wider(
    id_cols = apt_id,
    names_from = unique_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()

stopifnot(identical(colnames(cunningham$cm), cunningham$meta$unique_id))

# transpose for glmnet
cunningham$cm <- t(cunningham$cm)

# filter for final model variables
cunningham$cm <- cunningham$cm[,
  colnames(cunningham$cm) %in% row.names(model$final_coef)
]

cunningham$predictions <- data.frame(
  prediction = predict(model, newdata = cunningham$cm, type = 'response')
) %>%
  rename(prediction = 1) %>%
  rownames_to_column("unique_id") %>%
  left_join(cunningham$meta[, c("unique_id", "dx")], by = "unique_id")
# export predictions
saveRDS(
  cunningham$predictions,
  "04_prediction_model/data/cunningham.predictions.rds"
)

cunningham$predictions %>%
  ggplot(aes(x = dx, y = prediction)) +
  geom_jitter()


# COVID19 (Gisby) ---------------------------------------------------------
gisby <- list()
gisby$meta <- readRDS("04_prediction_model/data/gisby.meta.rds")
gisby$rfu <- read_parquet("04_prediction_model/data/gisby.tidy_counts.parquet")

# drop somamers with no apt_id
gisby$rfu <- gisby$rfu %>% drop_na(apt_id)
# convert tidy rfu df to log2 rfu count matrix
gisby$cm <- gisby$rfu %>%
  pivot_wider(
    id_cols = apt_id,
    names_from = sample_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()
# match count matrix sample order to metadata
gisby$cm <- gisby$cm[, match(gisby$meta$sample_id, colnames(gisby$cm))]
stopifnot(identical(colnames(gisby$cm), gisby$meta$sample_id))

# transpose for glmnet
gisby$cm <- t(gisby$cm)

# filter for final model variables
gisby$cm <- gisby$cm[, colnames(gisby$cm) %in% row.names(model$final_coef)]

gisby$predictions <- data.frame(
  prediction = predict(model, newdata = gisby$cm, type = 'response')
) %>%
  rename(prediction = 1) %>%
  rownames_to_column("sample_id") %>%
  left_join(gisby$meta[, c("sample_id", "WHO_severity")], by = "sample_id")
saveRDS(gisby$predictions, "04_prediction_model/data/gisby.predictions.rds")
gisby$predictions %>%
  ggplot(aes(x = WHO_severity, y = prediction)) +
  geom_jitter()
