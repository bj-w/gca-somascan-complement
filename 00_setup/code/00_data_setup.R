library(data.table)
library(tidyverse)
library(readxl)
library(janitor)
library(SomaDataIO)
library(arrow)


# data import -------------------------------------------------------------
meta <- read_xlsx("raw_data/TABUL_Soma_Ptnt_Metadata.xlsx") %>%
  clean_names()
adat <- read_adat(
  "raw_data/ICL-SS-2343676_2025-02-21/SS-2343676_v4.1_CitratePlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
)

# metadata ----------------------------------------------------------------

# array sample QC metrics
array_qc <- data.frame(
  study_id = adat$SampleId,
  plate = adat$PlateId,
  plate_position = adat$PlatePosition,
  slide = adat$SlideId,
  row_check = adat$RowCheck
)

# add to meta
meta <- meta %>% left_join(array_qc, by = "study_id")

# extract site from study ID and add as column
meta$site <- str_split_i(meta$study_id, pattern = "-", i = 2)

# estimate freezer time based on sample collection date
meta <- meta %>%
  left_join(
    meta[, c("study_id", "plasma_sample_date")] %>%
      left_join(
        tibble(study_id = adat$SampleId, plate_run_date = adat$PlateRunDate),
        by = "study_id"
      ) %>%
      mutate(plate_run_date = as_datetime(plate_run_date)) %>%
      mutate(
        freezer_months = interval(plasma_sample_date, plate_run_date) %/%
          months(1)
      ) %>%
      select(study_id, freezer_months),
    by = "study_id"
  )


# format
meta <- meta %>%
  mutate(
    sex = factor(sex),
    # recode dx to make syntactically valid downstream
    ref_dx = case_match(ref_dx, "C-GCA" ~ "cGCA", "Not_C-GCA" ~ "not_cGCA"),
    ref_dx = factor(ref_dx, levels = c("not_cGCA", "cGCA")),
    plate = factor(plate),
    row_check = factor(row_check),
    site = factor(site)
  )


# analyte metadata --------------------------------------------------------
analytes <- getAnalyteInfo(adat)
# restrict to human analytes
analytes <- analytes %>% filter(Organism == "Human")
# remove "INTERNAL USE ONLY PROTEINS"
analytes <- analytes %>%
  filter(!str_detect(TargetFullName, "Internal Use Only"))
# remove hybridisation controls
analytes <- analytes %>% filter(!str_detect(SomaId, "HCE"))
# create column combining target with aptname
analytes <- analytes %>% mutate(apt_id = paste0(EntrezGeneSymbol, "_", AptName))
# create column combining gene symbol with numeric aptamer ID for nice labelling on plots
analytes <- analytes %>%
  mutate(apt_label = paste0(EntrezGeneSymbol, " (", SeqId, ")"))

# ensembl id - get from uniprot mapping file
ensembl_ids <- fread(
  "raw_data/20260716_uniprot_HUMAN_9606_idmapping_selected.tab.gz"
)
ensembl_ids <- ensembl_ids %>%
  select(UniProt = V1, ensembl_id = V19) %>%
  # remove ensembl id version
  mutate(
    # keep first recorded ensembl and enztrez gene IDs
    ensembl_id = str_remove(ensembl_id, ";.*"),
    # remove trailing version number from ensembl gene IDs
    ensembl_id = str_remove(ensembl_id, "\\..*")
  )

ensembl_ids <- analytes %>%
  select(AptName, UniProt) %>%
  separate_longer_delim(UniProt, delim = "|") %>%
  left_join(ensembl_ids, by = "UniProt") %>%
  summarise(ensembl_id = paste(ensembl_id, collapse = "|"), .by = AptName)

analytes <- analytes %>%
  left_join(ensembl_ids, by = "AptName") %>%
  relocate(ensembl_id, .after = EntrezGeneSymbol)

# export
saveRDS(analytes, "00_setup/data/metadata.analytes.rds")

# counts ------------------------------------------------------------------
counts <- adat %>%
  as.data.frame() %>%
  filter(SampleType == "Sample") %>%
  select(SampleId, contains("seq.")) %>%
  as_tibble()

counts_tidy <- counts %>%
  pivot_longer(cols = !SampleId, names_to = "AptName", values_to = "RFU") %>%
  inner_join(analytes[, c("AptName", "apt_id")], by = "AptName") # inner join to exclude unecessary targets (e.g. non human proteins)
# count matrix
count_matrix <- counts_tidy %>%
  pivot_wider(id_cols = apt_id, names_from = SampleId, values_from = RFU) %>%
  column_to_rownames("apt_id")
# match sample IDs in count matrix to the order of sample IDs in metadata
count_matrix <- count_matrix[, match(meta$study_id, colnames(count_matrix))]
stopifnot(identical(meta$study_id, names(count_matrix)))

# export count matrix
count_matrix %>%
  rownames_to_column("apt_id") %>%
  write_parquet("00_setup/data/count_matrix.preQC.parquet")
# export log2 transformed count matrix
count_matrix %>%
  log2() %>%
  rownames_to_column("apt_id") %>%
  write_parquet("00_setup/data/count_matrix.log2.preQC.parquet")
# export in tidy format - normalised + log2 transformed
count_matrix %>%
  rownames_to_column("apt_id") %>%
  pivot_longer(cols = !apt_id, names_to = "study_id", values_to = "RFU") %>%
  mutate(log2_RFU = log2(RFU)) %>%
  write_parquet("00_setup/data/count_tidy.preQC.parquet")

# add on abundance metrics to metadata
abundance_metrics <- counts_tidy %>%
  mutate(log2_RFU = log2(RFU)) %>%
  summarise(
    median_RFU = median(RFU),
    median_log2_RFU = median(log2_RFU),
    total_RFU = sum(RFU),
    total_log2_RFU = sum(log2_RFU),
    iqr_RFU = IQR(RFU),
    iqr_log2_RFU = IQR(log2_RFU),
    .by = SampleId
  )
meta <- left_join(meta, abundance_metrics, by = c("study_id" = "SampleId"))
# export sample metadata
saveRDS(meta, "00_setup/data/metadata.samples.preQC.rds")
