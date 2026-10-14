library(WGCNA)
library(tidyverse)
library(arrow)
source("functions.R")
theme_set(figTheme())
set.seed(1997)
# data import -------------------------------------------------------------
meta <- readRDS("00_setup/data/metadata.samples.postQC.rds")
counts <- read_parquet("00_setup/data/count_tidy.postQC.parquet")

# format counts -----------------------------------------------------------
counts <- counts %>%
  pivot_wider(
    id_cols = "apt_id",
    names_from = study_id,
    values_from = log2_RFU
  ) %>%
  column_to_rownames("apt_id") %>%
  as.matrix()
stopifnot(identical(colnames(counts), meta$study_id))
# transpose + regress covariates
counts <- resid(lm(
  t(counts) ~ poly(age_at_dx, 2) +
    sex +
    poly(iqr_log2_RFU, 2) +
    poly(gc_days, 2),
  data = meta
))
# counts %>% t() %>% as.data.frame() %>% rownames_to_column("sample_id") %>%
#   write_parquet("")
sample_clust <- hclust(dist(counts), method = "average")
saveRDS(sample_clust, "02_protein_modules/data/sample_hclust.rds")

# wgcna -------------------------------------------------------------------
# pick soft thresholding powers
powers <- c(c(1:10), seq(from = 10, to = 30, by = 1))
sft <- pickSoftThreshold(
  counts,
  powerVector = powers,
  verbose = 5,
  networkType = "signed",
  blockSize = 10000
)
sft_df <- sft$fitIndices %>%
  mutate(model_fit = -sign(slope) * SFT.R.sq)
power <- sort(sft_df[sft_df$model_fit > 0.85, ]$Power)[1]
sft <- list(sft = sft_df, power = power)
saveRDS(sft, "02_protein_modules/data/soft_thresholding.rds")


net <- blockwiseModules(
  datExpr = counts,
  power = power,
  TOMType = "signed",
  networkType = "signed",
  deepSplit = 3,
  maxBlockSize = 10000,
  numericLabels = TRUE,
  randomSeed = 1997,
  verbose = 3
)
saveRDS(net, "02_protein_modules/data/network.rds")

# connectivity -----------------------------------------------------------
# ie correlation between protein abundance and the module eigenprotein.
signed_connectivity <- signedKME(
  datExpr = counts,
  datME = net$MEs,
  outputColumnName = "ME"
) %>%
  rownames_to_column("apt_id") %>%
  pivot_longer(cols = !apt_id, names_to = "module", values_to = "pearson")


# an alternative method is below that also provides the p value
# Calculate all correlations at once
cor_matrix <- cor(counts, net$MEs, use = "complete.obs")

# collect p values
get_cor_pvals <- function(x, y) {
  n <- nrow(x)
  cor_mat <- cor(x, y, use = "complete.obs")
  # Convert correlations to t-statistics
  t_stats <- cor_mat * sqrt((n - 2) / (1 - cor_mat^2))
  # Get p-values
  p_vals <- 2 * pt(abs(t_stats), df = n - 2, lower.tail = FALSE)
  return(p_vals)
}
pval_matrix <- get_cor_pvals(counts, net$MEs)
# Convert to long format
correlation_results <- expand.grid(
  apt_id = colnames(counts),
  module = colnames(net$MEs)
)
correlation_results$pearson <- as.vector(cor_matrix)
correlation_results$pval <- as.vector(pval_matrix)
correlation_results

# prove equality
equal_test <- signed_connectivity %>%
  full_join(correlation_results, by = c("apt_id", "module"))
stopifnot(all.equal(equal_test$pearson.x, equal_test$pearson.y))

# export
write_parquet(
  correlation_results,
  "02_protein_modules/data/connectivity.parquet"
)
