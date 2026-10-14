library(arrow)
library(msigdbr)
library(fgsea)
library(tidyverse)
set.seed(1997)
conflicted::conflicts_prefer(
  dplyr::filter(),
  dplyr::group_by(),
  dplyr::mutate()
)
source("functions.R")
# import results ----------------------------------------------------------
res <- read_csv("01_differential_abundance/data/limma.ref_dx.csv")

res <- res %>%
  arrange(P.Value) %>%
  distinct(EntrezGeneID, .keep_all = TRUE)


# gene sets ---------------------------------------------------------------
gs <- prepareGeneSets(
  geneIdentifier = "ncbi_gene",
  reactome = TRUE,
  hallmark = TRUE,
  go_bp = TRUE
)


gs <- map(gs, \(x) {
  addComplexEntrezIdsToGeneSet(differentialAbundance = res, geneSet = x)
})

# export a copy of the SOMAmers used
saveRDS(
  res[, "apt_id"],
  "01_differential_abundance/data/somamers_used_gsea.rds"
)

# run GSEA ----------------------------------------------------------------
gsea <- map(names(gs), \(gene_set) {
  runGSEA(results = res, geneSet = gs[[gene_set]])
}) %>%
  set_names(names(gs))
saveRDS(gsea, "01_differential_abundance/data/gsea.rds")
