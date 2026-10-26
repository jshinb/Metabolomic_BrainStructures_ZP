# ------------------------------------------
# Sensitivity analysis: downsample ADO1 to SPS sample size
# ------------------------------------------

set.seed(2026)

# Example input objects from your script
# ADO1_dat: ADO1 analytic dataset
# SPS_dat: SPS analytic dataset
# trait: e.g. "Volume.WM.adjAgeSexICV" or "T1wSI.WM.adjAgeSex"
# metabolite_list: vector of key metabolites from ADO1 findings

# function to run the same association model under a matched sample size
run_matched_n_test <- function(dat_large, dat_small, trait_name, metabolite_name, n_iter = 500) {
  out <- list()
  #dat_small used for the size only
  for (i in seq_len(n_iter)) {
    analdat_large = dat_large %>% select(all_of(c('uniqueID',trait_name,metabolite_name,'age.c','age.c2','Sex')))
    analdat_small = dat_small %>% select(all_of(c('uniqueID',trait_name,metabolite_name,'age.c','age.c2','Sex')))
    match.size = nrow(na.omit(analdat_small))
    dat_large = na.omit(analdat_large)
    
    idx <- sample(seq_len(nrow(dat_large)), size = match.size, replace = FALSE)
    dat <- dat_large[idx, ]
    
    dat[["y.rnt"]] <- inormal(dat[[trait_name]])
    dat[["m"]] <- inormal(dat[[metabolite_name]])
    
    fit <- lm(y.rnt ~ m + (age.c + age.c2) * Sex, data = dat)
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
