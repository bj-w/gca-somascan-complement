library(arrow)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(furrr)
set.seed(1997)
# conflicted::conflicts_prefer(dplyr::filter(), dplyr::group_by(), dplyr::mutate())
source("functions.R")
# data import -------------------------------------------------------------
net <- readRDS("02_protein_modules/data/network.rds")
res <- read_csv("01_differential_abundance/data/limma.ref_dx.csv")


# extract module memberships
module_aptamers <- data.frame(module = net$colors) %>%
  rownames_to_column("apt_id") %>%
  mutate(module = paste0("ME", module)) %>%
  # remove module 0
  filter(!module == "ME0") %>%
  # add on
  left_join(res[, c("apt_id", "EntrezGeneID")], by = "apt_id") %>%
  # format so that gene set is 1st column, enztrez is 2nd column
  select(module, EntrezGeneID)
# keep unique SOMAmers
res <- res %>%
  arrange(P.Value) %>%
  distinct(EntrezGeneID, .keep_all = TRUE)


# generic function for module overrepresentation --------------------------
moduleORA <- function(
  deResults,
  module_index,
  geneset_index,
  modules = NULL,
  genesets = NULL,
  nCores = 3
) {
  require(furrr)
  # If not specified, use all modules and genesets
  if (is.null(modules)) {
    modules <- names(module_index)
  }
  if (is.null(genesets)) {
    genesets <- names(geneset_index)
  }

  # Create all combinations
  combinations <- expand.grid(
    module = modules,
    geneset = genesets,
    stringsAsFactors = FALSE
  )

  # Pre-extract EntrezGeneID once for all comparisons
  entrez_ids <- deResults$EntrezGeneID

  # Vectorized batch processing
  plan(multisession, workers = nCores)
  results <- future_pmap_dfr(
    combinations,
    function(module, geneset) {
      module_genes <- module_index[[module]]
      geneset_genes <- geneset_index[[geneset]]

      if (is.null(module_genes) || is.null(geneset_genes)) {
        return(tibble(
          module = module,
          geneset = geneset,
          estimate = NA,
          p.value = NA,
          conf.low = NA,
          conf.high = NA,
          method = NA,
          alternative = NA
        ))
      }

      in_module <- entrez_ids %in% module_genes
      in_geneset <- entrez_ids %in% geneset_genes

      contingency_table <- matrix(
        c(
          sum(in_module & in_geneset),
          sum(in_module & !in_geneset),
          sum(!in_module & in_geneset),
          sum(!in_module & !in_geneset)
        ),
        nrow = 2,
        byrow = TRUE
      )

      fisher_result <- fisher.test(contingency_table)

      tibble(
        module = module,
        geneset = geneset,
        estimate = fisher_result$estimate,
        p.value = fisher_result$p.value,
        conf.low = fisher_result$conf.int[1],
        conf.high = fisher_result$conf.int[2],
        method = fisher_result$method,
        alternative = fisher_result$alternative
      )
    },
    .progress = TRUE
  )
  plan(sequential)
  results <- arrange(results, p.value)
  return(results)
}

# overrepresentation of biological gene sets ------------------------------
# collect gene sets
gs <- prepareGeneSets(
  geneIdentifier = "ncbi_gene",
  hallmark = TRUE,
  reactome = TRUE
)
# add on multiple gene ids mapped to somamer if there's consensus
gs <- map(gs, \(x) {
  addComplexEntrezIdsToGeneSet(differentialAbundance = res, geneSet = x)
})

geneset_ora <- map(names(gs), \(database) {
  map(unique(module_aptamers$module), \(moduleNumber) {
    fora(
      pathways = split(
        gs[[database]][["ncbi_gene"]],
        gs[[database]][["gs_name"]]
      ),
      genes = filter(module_aptamers, module == moduleNumber)$EntrezGeneID,
      universe = res$EntrezGeneID
    ) |>
      mutate(module = moduleNumber, .before = pathway)
  }) |>
    list_rbind() |>
    as_tibble()
}) |>
  set_names(names(gs))


# export
saveRDS(geneset_ora, "02_protein_modules/data/geneset_ora.wgcna_modules.rds")
