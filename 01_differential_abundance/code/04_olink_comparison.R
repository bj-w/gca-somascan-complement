library(arrow)
library(broom)
library(tidyverse)
library(lubridate)
source("functions.R")

# DATA IMPORT -------------------------------------------------------------
# meta
meta <- readRDS("00_setup/data/metadata.samples.postQC.rds")
analytes <- readRDS("00_setup/data/metadata.analytes.rds")


# rfu values
rfu <- read_parquet("00_setup/data/count_tidy.postQC.parquet")


# sample overlap ----------------------------------------------------------
olink_meta <- read_csv("raw_data/lvv_olink/Cohort2_metadata.csv") |>
  janitor::clean_names()
# TAB +
# make sure sample n matches (56 TAB+ C-GCA, 89 not C-GCA)
stopifnot(nrow(filter(olink_meta, subgroup == "C-GCA TAB+")) == 56)
stopifnot(nrow(filter(olink_meta, subgroup == "Not_C-GCA")) == 89)
olink_meta <- olink_meta |> filter(subgroup %in% c("Not_C-GCA", "C-GCA TAB+"))


olink_meta |> glimpse()
soma_meta <- meta |>
  select(
    sample_id = study_id,
    diagnosis = ref_dx,
    age_bin = age_at_dx,
    sex,
    ethnicity = ethnicity_dichotomised,
    freezer_time = freezer_months,
    cranial_isch_comp = cranial_isch_comp_rm,
    pmr = pmryesno,
    tab,
    uss
  ) |>
  mutate(freezer_time = freezer_time / 12) |>
  separate_wider_delim(
    cols = sample_id,
    delim = "-",
    names = c("drop1", "site", "drop2"),
    cols_remove = FALSE
  ) |>
  mutate(site = as.numeric(site)) |>
  select(-c(drop1, drop2)) |>
  mutate(
    diagnosis = replace_values(
      as.character(diagnosis),
      "cGCA" ~ "C-GCA",
      "not_cGCA" ~ "Not_C-GCA"
    ),
    age_bin = cut(age_bin, seq(0, 100, by = 10)),
    age_bin = str_replace(age_bin, "]", ")"),
    age_bin = str_replace(age_bin, ",", ", ")
  ) |>
  mutate_at(
    c("sex", "ethnicity", "cranial_isch_comp", "pmr", "tab", "uss"),
    str_to_title
  )

full_join(
  olink_meta,
  soma_meta,
  by = c(
    "site",
    "diagnosis",
    "sex",
    "age_bin",
    "ethnicity",
    "cranial_isch_comp",
    "pmr",
    "tab",
    "uss"
  )
) |>
  view()

# protein overlap ---------------------------------------------------------
olink_dap <- readxl::read_xlsx(
  "raw_data/lvv_olink/MaughanEtAl2025_supp.xlsx",
  sheet = 14,
  skip = 2
) |>
  janitor::clean_names()

rfu |> distinct(apt_id) |> left_join(analytes, by = "apt_id")

# protein association replication -----------------------------------------
soma_dap <- rfu |>
  filter(
    apt_id %in% filter(analytes, ensembl_id %in% olink_dap$ensembl_id)$apt_id
  ) |>
  left_join(select(meta, study_id, ref_dx), by = "study_id") |>
  group_by(apt_id) |>
  do(tidy(lm(log2_RFU ~ ref_dx, data = .))) |>
  ungroup() |>
  filter(term == "ref_dxcGCA") |>
  select(-term, log2fc = estimate, std_error = std.error, pval = p.value) |>
  mutate(adj_pval = p.adjust(pval, method = "fdr")) |>
  left_join(analytes, by = "apt_id") |>
  arrange(pval)

saveRDS(soma_dap,"01_differential_abundance/data/lvv_olink_soma_replication_lm.rds"
)
