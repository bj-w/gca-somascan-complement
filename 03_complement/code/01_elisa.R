library(tidyverse)
library(patchwork)
library(broom)
source("functions.R")
theme_set(figTheme())

# elisa concentrations
elisa <- list(
  C3a = read_csv("raw_data/elisa/TABUL_C3a_ELISA.csv") |>
    mutate(concentration = times_th),
  C5a = read_csv("raw_data/elisa/TABUL_C5a_ELISA.csv") |>
    mutate(concentration = times_twenty),
  sC5b9 = read_csv("raw_data/elisa/TABUL_C5b9_ELISA.csv") |>
    mutate(concentration = times_ten),
  CFHR5 = read_csv("raw_data/elisa/TABUL_CFHR5_ELISA_2_2.csv") |>
    mutate(concentration = times_os),
  fB = read_csv("raw_data/elisa/TABUL_fB_ELISA.csv") |>
    mutate(concentration = times_tt),
  fH = read_csv("raw_data/elisa/TABUL_fH_ELISA.csv") |>
    mutate(concentration = times_ot)
)

# metadata
meta <- read_csv("raw_data/TABUL_Soma_Ptnt_Metadata.csv") |>
  janitor::clean_names()
meta$ref_dx
meta <- meta |>
  mutate(
    ref_dx = replace_values(ref_dx, "Not_C-GCA" ~ "non-GCA", "C-GCA" ~ "GCA"),
    ref_dx = factor(ref_dx, levels = c("non-GCA", "GCA"))
  )

# format
elisa <- elisa |>
  list_rbind(names_to = "target") |>
  janitor::clean_names() |>
  rename(elisa_plate = plate) |>
  mutate(study_id = paste0("EP-", study_id)) |>
  mutate(log2_concentration = log2(concentration)) |>
  mutate(
    zscore_log2_concentration = scale(
      log2_concentration,
      center = TRUE,
      scale = TRUE
    ),
    .by = target
  )

# export
saveRDS(elisa, "03_complement/data/elisa.tidy.rds")
saveRDS(meta, "03_complement/data/elisa.metadata.rds")
