library(tidyverse)
library(readxl)

# import metadata and restrict to gca cases
# will use raw so Martina can have a look at samples run via ELISA but not
# SomaScan
meta <- readRDS("00_setup/data/metadata.samples.preQC.rds") %>%
  filter(ref_dx == "cGCA")
# import histology data
master <- read_xlsx("raw_data/Olink histology data 06102023.xlsx") %>%
  janitor::clean_names() %>%
  filter(study_id %in% meta$study_id)
# new df for recoding and rename columns for ease
recoded <- master %>%
  rename(
    barcelona_score = barcelona_score_0_normal_1_adventitial_2_adventitial_invasive_3_concentric_bilayer_4_panarteritis,
    media_pattern = media_pattern_0_normal_1_focal_2_multifocal_3_diffuse,
    occlusion_grade = occlusion_grade_0_none_1_50_percent_2_50_75_percent_3_75_100_percent_4_complete
  ) %>%
  mutate(occlusion_grade = as.numeric(occlusion_grade)) # parses as character for some reason?

# fill in GBA present blanks
# if scored as no gca present, subsequent scores are 0
recoded <- recoded %>%
  mutate_at(
    c(
      "barcelona_score",
      "media_pattern",
      "media_destruction",
      "occlusion_grade"
    ),
    ~ replace(., is.na(.), 0)
  )


# take max score from each sample
recoded <- recoded %>%
  mutate(
    max_barcelona_score = max(barcelona_score),
    max_media_pattern = max(media_pattern),
    max_media_destruction = max(media_destruction),
    max_occlusion_grade = max(occlusion_grade),
    .by = study_id
  )


# barcelona score ---------------------------------------------------------
# definition: 0 = normal, 1 = adventitial, 2 = adventitial invase, 3 = concentric bilayer, 4 = panarteritis
# recode: extensive? no <= 2, yes > 2
recodeBarcelonaScore <- function(x) {
  x <- if_else(
    x > 2,
    "yes",
    "no"
  )
  x <- factor(x, levels = c("no", "yes"))
  return(x)
}

recoded <- recoded %>%
  mutate(
    barcelona_score_recode = recodeBarcelonaScore(barcelona_score),
    max_barcelona_score_recode = recodeBarcelonaScore(max_barcelona_score)
  )


# media pattern -----------------------------------------------------------
# dictionary: 0 = normal, 1 = focal, 2 = multifocal, 3 = diffuse
# recoding: 0|1 = normal/focal, 2|3 = multifocal/diffuse
recodeMediaPattern <- function(x) {
  x <- case_match(
    x,
    c(0, 1) ~ "normal_focal",
    c(2, 3) ~ "multifocal_diffuse"
  )
  x <- factor(x, levels = c("normal_focal", "multifocal_diffuse"))
  return(x)
}

recoded <- recoded %>%
  mutate(
    media_pattern_recode = recodeMediaPattern(media_pattern),
    max_media_pattern_recode = recodeMediaPattern(max_media_pattern)
  )


# media destruction -------------------------------------------------------
# recoding: 0 = no, 1 = yes
# 0 = no,  = yes
recodeMediaDestruction <- function(x) {
  x <- case_match(
    x,
    0 ~ "no",
    1 ~ "yes"
  )
  x <- factor(x, levels = c("no", "yes"))
  return(x)
}

recoded <- recoded %>%
  mutate(
    media_destruction_recode = recodeMediaDestruction(media_destruction),
    max_media_destruction_recode = recodeMediaDestruction(max_media_destruction)
  )

# occlusion grade ---------------------------------------------------------
# dictionary: 0 = none, 1 = 50%, 2 = 50-75%, 3 = 75-100%, 4 = complete
# recoding: <50% = low, >50% = high
recodeOcclusionGrade <- function(x) {
  x <- case_match(
    x,
    c(0, 1) ~ "low",
    c(2, 3, 4) ~ "high"
  )
  x <- factor(x, levels = c("low", "high"))
  return(x)
}

recoded <- recoded %>%
  mutate(
    occlusion_grade_recode = recodeOcclusionGrade(occlusion_grade),
    max_occlusion_grade_recode = recodeOcclusionGrade(max_occlusion_grade)
  )


# export ------------------------------------------------------------------
saveRDS(recoded, "00_setup/data/meta.histology.rds")
# tidied up version for Martina
tidy_recoded <- recoded %>%
  select(
    study_id,
    max_barcelona_score_recode,
    max_media_pattern_recode,
    max_media_destruction_recode,
    max_occlusion_grade_recode
  ) %>%
  distinct()
saveRDS(tidy_recoded, "00_setup/data/meta.histology_tidy.rds")
