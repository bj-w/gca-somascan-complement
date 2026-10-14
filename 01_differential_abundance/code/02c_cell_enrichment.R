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

# blood -------------------------------------------------------------------
blood <- read_parquet("01_differential_abundance/data/blood_nTPM.parquet") %>%
  column_to_rownames("gene_id")
output$blood <- runTissueEnrich(
  reference = "custom",
  tpmMatrix = blood,
  inputGenes = gs$inputGenes,
  backgroundGenes = gs$backgroundGenes,
  categoryIdentifier = "cell_type"
)

# vasculature -------------------------------------------------------------
vasculature <- read_parquet(
  "01_differential_abundance/data/vasculature_nTPM.parquet"
) %>%
  column_to_rownames("gene_id")
output$vasculature <- runTissueEnrich(
  reference = "custom",
  tpmMatrix = vasculature,
  inputGenes = gs$inputGenes,
  backgroundGenes = gs$backgroundGenes,
  categoryIdentifier = "cell_type"
)
output$vasculature$source_data |> filter(fdr < 0.05)
# GCA associated proteins overrepresented in vasculature-cell specific genes
output$vasculature$categories |>
  filter(Tissue == "fibroblast") |>
  inner_join(filter(res, adj.P.Val < 0.05), by = c("Gene" = "ensembl_id")) |>
  arrange(P.Value)
# a lot of acute phase + complement appears fibroblast specific - liver derived in plasma. Look for enrichment among
# tissue specific genes
# import tissue specific results
tissue_enrichment <- readRDS(
  "01_differential_abundance/data/tissueEnrichment.rds"
)$consensus$categories |>
  filter(!Tissue == "All")
# ora
output$fibroblast_tissue_ora <- fora(
  pathways = split(tissue_enrichment$Gene, tissue_enrichment$Tissue),
  genes = output$vasculature$categories |>
    filter(
      Tissue == "fibroblast",
      Gene %in% filter(res, adj.P.Val < 0.05)$ensembl_id
    ) |>
    pull(Gene),
  universe = res$ensembl_id
) |>
  as_tibble()

# vasculature fibroblast specific genes overrepresented in liver specific genes
# rerun enrichment removing liver specific genes
liver_specific <- tissue_enrichment |> filter(Tissue == "liver")
gs_excl_liver <- formatTissueEnrichGeneSets(
  limmaResultDf = filter(res, !ensembl_id %in% liver_specific$Gene)
)
output$vasculature_excl_liver <- runTissueEnrich(
  reference = "custom",
  # also exclude from reference matrix
  tpmMatrix = vasculature[!row.names(vasculature) %in% liver_specific$Gene, ],
  inputGenes = gs_excl_liver$inputGenes,
  backgroundGenes = gs_excl_liver$backgroundGenes,
  categoryIdentifier = "cell_type"
)

output$vasculature_excl_liver$categories |>
  filter(Tissue == "fibroblast") |>
  inner_join(filter(res, adj.P.Val < 0.05), by = c("Gene" = "ensembl_id")) |>
  arrange(P.Value)

# export ------------------------------------------------------------------
saveRDS(output, "01_differential_abundance/data/cellEnrichment.rds")
