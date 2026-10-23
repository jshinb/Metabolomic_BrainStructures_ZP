# ------------------------------------------
# Sensitivity analysis: downsample ADO1 to SPS sample size
# ------------------------------------------

set.seed(2026)

# Example input objects from your script
# ADO1_dat: ADO1 analytic dataset
# SPS_dat: SPS analytic dataset
# trait: e.g. "Volume.WM.adjAgeSexICV" or "T1wSI.WM.adjAgeSex"
# metabolite_list: vector of key metabolites from ADO1 findings

# Use adult sample size as the target
n_target <- nrow(SPS_dat)

# function to run the same association model under a matched sample size
run_matched_n_test <- function(dat_large, dat_small, trait_name, metabolite_name, n_iter = 500) {
  out <- list()
  
  for (i in seq_len(n_iter)) {
    idx <- sample(seq_len(nrow(dat_large)), size = nrow(dat_small), replace = FALSE)
    dat <- dat_large[idx, ]
    
    dat[["y.rnt"]] <- inormal(dat[[trait_name]])
    dat[["m"]] <- inormal(dat[[metabolite_name]])
    
    fit <- lm(y.rnt ~ m + (age.c + age.c2) + Sex, data = dat)
    coef_tmp <- summary(fit)$coef["m", ]
    
    out[[i]] <- data.frame(
      iter = i,
      beta = as.numeric(coef_tmp["Estimate"]),
      se = as.numeric(coef_tmp["Std. Error"]),
      pvalue = as.numeric(coef_tmp["Pr(>|t|)"]),
      n = nrow(dat),
      sign = sign(as.numeric(coef_tmp["Estimate"]))
    )
  }
  
  out <- bind_rows(out)
  out
}

# Example: For a key metabolite known to be significant in ADO1
# Replace with actual ADO1 metabolite name
metab_name <- "PE(16:0/18:1)"  # example only
trait_name <- "T1wSI.WM.adjAgeSex"

ADO1_matched <- run_matched_n_test(
  dat_large = ADO1_dat,
  dat_small = SPS_dat,
  trait_name = trait_name,
  metabolite_name = metab_name,
  n_iter = 500
)

# Summaries
summary_ADO1 <- ADO1_matched %>%
  summarise(
    median_beta = median(beta, na.rm = TRUE),
    mean_p = mean(pvalue, na.rm = TRUE),
    prop_sig = mean(pvalue < 0.05, na.rm = TRUE),
    prop_same_direction = mean(sign == sign(median(beta, na.rm = TRUE)), na.rm = TRUE)
  )

summary_ADO1

# key multiple metabolites from ADO1: for example top metabolites by p-value
key_metabs <- c("PE(16:0/18:1)", "PC(16:0/16:1)", "SM(d18:1/16:0)")

# trait names of interest
traits <- c("Volume.WM.adjAgeSexICV", "T1wSI.WM.adjAgeSex", "MTR.WM.adjAgeSex")

results_list <- list()

for (trait_name in traits) {
  for (metab_name in key_metabs) {
    matched_res <- run_matched_n_test(
      dat_large = ADO1_dat,
      dat_small = SPS_dat,
      trait_name = trait_name,
      metabolite_name = metab_name,
      n_iter = 500
    )
    
    results_list[[length(results_list) + 1]] <- matched_res %>%
      summarise(
        trait = trait_name,
        metabolite = metab_name,
        median_beta = median(beta, na.rm = TRUE),
        prop_sig = mean(pvalue < 0.05, na.rm = TRUE),
        prop_same_direction = mean(sign == sign(median(beta, na.rm = TRUE)), na.rm = TRUE),
        mean_p = mean(pvalue, na.rm = TRUE)
      )
  }
}

matched_summary <- bind_rows(results_list)
matched_summary
