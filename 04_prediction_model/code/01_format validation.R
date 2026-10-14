library(readxl)
library(tidyverse)
library(data.table)
library(janitor)
library(arrow)

analytes <- readRDS("00_setup/data/metadata.analytes.rds")
# Cunningham GCA ----------------------------------------------------------
# metadata
cunningham <- list()
cunningham$meta <- read_xlsx(
  "raw_data/gca_cunningham/metadata.xlsx",
  skip = 1
) %>%
  clean_names() %>%
  # exclude parsing error
  filter(!str_detect(subject_id, "aBMI"))

# format metadata
cunningham$meta <- cunningham$meta %>%
  mutate(
    subject_id = str_replace(subject_id, "-", "_"),
    study_group = case_match(
      study_group,
      "Active GCA" ~ "active",
      "Inactive GCA" ~ "inactive",
      .default = study_group
    ),
    gca_classification = case_match(
      gca_classification,
      "Cranial GCA" ~ "Cranial",
      "Large vessel vasculitis" ~ "LV",
      "-" ~ NA,
      .default = gca_classification
    ),
    sex = factor(sex, levels = c("Male", "Female")),
    dx = paste0(study_group, "_", gca_classification)
  ) %>%
  # repeated samples (active/inactive) so create a unique id
  mutate(unique_id = paste0(subject_id, "_", study_group)) %>%
  relocate(unique_id, dx)

# export
saveRDS(cunningham$meta, "04_prediction_model/data/cunningham.meta.rds")


# RFU
cunningham$rfu <- read_xlsx("raw_data/gca_cunningham/rfu.xlsx", skip = 1) %>%
  rename(subject_id = `Subject ID`, study_group = `Study group`) %>%
  # recode
  mutate(
    subject_id = str_replace(subject_id, "-", "_"),
    study_group = case_match(
      study_group,
      "Active GCA" ~ "active",
      "Inactive GCA" ~ "inactive",
      .default = study_group
    )
  ) %>%
  # repeated samples (active/inactive) so create a unique id
  mutate(unique_id = paste0(subject_id, "_", study_group)) %>%
  relocate(unique_id) %>%
  # just use unique ID
  select(-c(subject_id, study_group)) %>%
  pivot_longer(cols = !unique_id, names_to = "seq_id", values_to = "RFU") %>%
  mutate(log2_RFU = log2(RFU))

# match up somamer ids
cunningham$analytes <- analytes %>%
  # remove seq. from AptName to get unique somamer id
  # format: [number.number]
  mutate(AptName = str_remove(AptName, "seq."))
# add somamer id in the format: [number.number]
cunningham$rfu <- cunningham$rfu %>% mutate(AptName = sub(".*_", "", seq_id))
cunningham$rfu <- cunningham$rfu %>%
  left_join(cunningham$analytes[, c("AptName", "apt_id")], by = "AptName")


# abundance metrics - no point adding on as data was quantile normalised prior
# to analysis according to paper
cunningham$rfu %>%
  summarise(
    median_RFU = median(RFU),
    median_log2_RFU = median(log2_RFU),
    total_RFU = sum(RFU),
    total_log2_RFU = sum(log2_RFU),
    iqr_RFU = IQR(RFU),
    iqr_log2_RFU = IQR(log2_RFU),
    .by = unique_id
  )

# export
write_parquet(
  cunningham$rfu,
  "04_prediction_model/data/cunningham.tidy_counts.parquet"
)


# Gisby covid19 -----------------------------------------------------------
gisby <- list()
# metadata is separate from waves 1 and 2. will just bind together
gisby$meta <- map(
  c("1", "2"),
  \(x) {
    read_csv(paste0("raw_data/covid19_gisby/w", x, "_metadata.csv")) %>%
      select(sample_id, sex, calc_age, WHO_severity)
  }
) %>%
  list_rbind() %>%
  # will remove the duplicated controls in both wave metadata files
  distinct() %>%
  # recode WHO severity to be negative if when NA
  mutate(WHO_severity = replace_na(WHO_severity, "negative"))
# format metadata
gisby$meta <- gisby$meta %>%
  mutate(
    sex = case_match(sex, "M" ~ "Male", "F" ~ "Female"),
    sex = factor(sex, levels = c("Male", "Female")),
    WHO_severity = factor(
      WHO_severity,
      levels = c("negative", "mild", "moderate", "severe", "critical")
    )
  )

# RFU
gisby$rfu <- fread("raw_data/covid19_gisby/soma_abundance.csv") %>%
  rename("sample_id" = "V1") %>%
  as_tibble()
# make tidy
gisby$rfu <- gisby$rfu %>%
  pivot_longer(cols = !sample_id, names_to = "AptName", values_to = "RFU")
# log2
gisby$rfu <- gisby$rfu %>% mutate(log2_RFU = log2(RFU))
# SOMAmer is of the format seq.[number].[number]
# equivalent is AptName in analytes df so can easily left join
gisby$rfu <- gisby$rfu %>%
  left_join(analytes[, c("AptName", "apt_id")], by = "AptName")
# export
write_parquet(gisby$rfu, "04_prediction_model/data/gisby.tidy_counts.parquet")


# calculate rfu metrics
gisby$rfu_metrics <- gisby$rfu %>%
  summarise(
    median_RFU = median(RFU),
    median_log2_RFU = median(log2_RFU),
    total_RFU = sum(RFU),
    total_log2_RFU = sum(log2_RFU),
    iqr_RFU = IQR(RFU),
    iqr_log2_RFU = IQR(log2_RFU),
    .by = sample_id
  )
# add rfu metrics to metadata
gisby$meta <- gisby$meta %>% left_join(gisby$rfu_metrics, by = "sample_id")
# just keep samples with proteomic data
gisby$meta <- gisby$meta %>% filter(sample_id %in% gisby$rfu$sample_id)
saveRDS(gisby$meta, "04_prediction_model/data/gisby.meta.rds")
