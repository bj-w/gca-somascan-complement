library(TissueEnrich)
library(tidyverse)
library(arrow)
library(fgsea)
source("functions.R")
theme_set(figTheme())

output <- list()

# import ------------------------------------------------------------------
res <- read_csv("01_differential_abundance/data/limma.ref_dx.csv") |>
  arrange(P.Value) %>%
  distinct(ensembl_id, .keep_all = TRUE) %>%
  drop_na(ensembl_id)

gs <- formatTissueEnrichGeneSets(limmaResultDf = res)


# consensus (HPA + GTEx) --------------------------------------------------
consensus <- read_tsv(
  "raw_data/tissue_enrichment/rna_tissue_consensus_20251208.tsv"
) %>%
  pivot_wider(id_cols = "Gene", names_from = "Tissue", values_from = "nTPM") %>%
  column_to_rownames("Gene") %>%
  as.matrix()

output$consensus <- runTissueEnrich(
  reference = "custom",
  tpmMatrix = consensus,
  inputGenes = gs$inputGenes,
  backgroundGenes = gs$backgroundGenes,
  categoryIdentifier = "tissue"
)

# secretome ---------------------------------------------------------------
secretome <- read_delim(
  "raw_data/tissue_enrichment/proteinatlas_20250513.tsv"
) %>%
  janitor::clean_names() %>%
  select(ensembl, secretome_location) %>%
  distinct()

output$secretome$ora <- fora(
  pathways = split(secretome$ensembl, secretome$secretome_location),
  genes = filter(res, adj.P.Val < 0.05)$ensembl_id,
  universe = res$ensembl_id
)
output$secretome$categories <- secretome

# export enrichments ------------------------------------------------------
saveRDS(output, "01_differential_abundance/data/tissueEnrichment.rds")
