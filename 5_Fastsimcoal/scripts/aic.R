# Model selection with the Akaike information criterion (AIC) for the best
# replicate of each population pair x demographic model.
#
# Run from this directory:  cd 5_Fastsimcoal/scripts && Rscript aic.R
#
# Input:  ../results/<rate>_rate/fastsimcoal2_bestlhoods_per_scenario_<rate>_mutation_rate.txt
#         (one row per pair x model: highest-likelihood replicate; `params` column =
#         number of parameters used for AIC)
# Output: ../results/<rate>_rate/model_selection_AIC_<rate>_mutation_rate.tsv
#         (pair, model, replicate, MaxEstLhood, params, AIC, deltaAIC, best)

for (rate in c("human", "mammal")) {
  res_dir <- file.path("..", "results", paste0(rate, "_rate"))
  fsc2_df <- read.delim(file.path(res_dir, paste0("fastsimcoal2_bestlhoods_per_scenario_", rate, "_mutation_rate.txt")))

  # fastsimcoal2 reports log10 likelihoods; convert to natural log
  fsc2_df$AIC <- 2*fsc2_df$params-2*(fsc2_df$MaxEstLhood/log10(exp(1)))

  fsc2_df$pair  <- sub("_m[1-4]_.*$", "", fsc2_df$X_group)
  fsc2_df$model <- sub("^.*_(m[1-4]_.*)$", "\\1", fsc2_df$X_group)
  fsc2_df$replicate <- basename(dirname(fsc2_df$X_source))

  fsc2_df$deltaAIC <- fsc2_df$AIC - ave(fsc2_df$AIC, fsc2_df$pair, FUN = min)
  fsc2_df$best <- fsc2_df$deltaAIC == 0

  out <- fsc2_df[order(fsc2_df$pair, fsc2_df$model),
                 c("pair", "model", "replicate", "MaxEstLhood", "params", "AIC", "deltaAIC", "best")]
  write.table(out, file.path(res_dir, paste0("model_selection_AIC_", rate, "_mutation_rate.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(out[out$best, c("pair", "model", "replicate")], row.names = FALSE)
}
