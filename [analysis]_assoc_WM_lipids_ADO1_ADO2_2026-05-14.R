#* Lipidomic association analyses of white matter phenotypes for ADO2
#* Cleaned and reorganized version
#* Last updated: 2026-10-23

# 0. Setup ----
# NOTE: This file must exist in your project directory and contains the shared
# working directory and common settings used across scripts.
project_specific_file <- "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"

source(project_specific_file)
setwd(wd)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(psych)
  library(lme4)
  library(RNOmni)
})

select <- dplyr::select
describe <- psych::describe

ado2datadir <- "00_ADO2_data"

# 1. Helper functions ----
clean_column_prefix <- function(df, prefix) {
  names(df) <- str_remove(names(df), prefix)
  df
}

summarise_missing <- function(df, id_cols = c("uniqueID", "fam_id", "sub_id")) {
  x <- df %>% select(-all_of(id_cols))
  na.count <- apply(as.data.frame(x), 1, is.na)
  na.count <- apply(na.count, 2, sum)
  table(na.count)
}

harmonize_metabolite_ids <- function(x) {
  x %>%
    str_replace("[-]FA", "FA") %>%
    str_replace("PE[(]P[-]", "PE(P") %>%
    str_replace("PE[(]O[-]", "PE(O")
}

generate_volcano_plot <- function(df_plot,
                                  trait_plot,
                                  generation_plot,
                                  signif.level,
                                  point.size = 1.5,
                                  point.alpha = 0.5) {
  df_plot %>%
    filter(trait == trait_plot, generation == generation_plot) %>%
    ggplot(aes(x = beta, y = -log10(pvalue),
               shape = fdr.signif,
               fill = fdr.signif,
               color = fdr.signif)) +
    geom_point(alpha = point.alpha, size = point.size) +
    scale_shape_manual(values = c("p.fdr<0.05" = 21, "p.fdr>=0.05" = 19)) +
    scale_color_manual(values = c("p.fdr<0.05" = "#073B4C", "p.fdr>=0.05" = "grey")) +
    scale_fill_manual(values = c("p.fdr<0.05" = "#073B4C", "p.fdr>=0.05" = "grey")) +
    theme(
      panel.background = element_rect(fill = scales::alpha("#f8edeb", 0.2), size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey95"),
      panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey95"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 1)
    ) +
    geom_hline(yintercept = -log10(signif.level), linetype = 2, color = "darkred") +
    geom_vline(xintercept = 0, color = "darkgrey")
}

# 2. ADO2 data import and cleaning ----
# 2A. Covariate and lipid data
covariate_lipid_data <- fread(file.path(wd, ado2datadir,
                                        "jean_shin_dataanalysis_20231027_ado2_20260508_111322.txt"))

covariate_data <- covariate_lipid_data %>%
  select(1:4, contains("ADO2-FieldModule"))
names(covariate_data) <- str_remove(names(covariate_data), "ADO2-FieldModule.")
covariate_data %>% describe

table(covariate_data$visit_batch)

# Lipid data
lipid_data <- covariate_lipid_data %>% select(-contains("ADO2-FieldModule"))
lipid_data2 <- fread(file.path(wd, ado2datadir,
                               "jean_shin_dataanalysis_20231027_ado2_20260508_153940.txt"))

names(lipid_data) <- str_remove(names(lipid_data), "ADO2-Lipidomic_june2021_species_conc_1.")
names(lipid_data2) <- str_remove(names(lipid_data2), "ADO2-Lipidomic_june2021_species_conc_2.")

lipid_data <- full_join(lipid_data, lipid_data2)
lipid_data <- lipid_data %>% mutate(across(-c(1:4), as.numeric))

# Remove individuals with no lipid data
na.count.sample <- lipid_data %>% select(-c(1:4)) %>% data.frame()
na.count.sample <- apply(na.count.sample, 1, is.na)
na.count.sample <- apply(na.count.sample, 2, sum)
print(table(na.count.sample))
lipid_data <- lipid_data %>% filter(na.count.sample != 957)
rm(na.count.sample)

# Set negative values to NA
lipid_data <- lipid_data %>%
  mutate(across(-c(1:4), ~ if_else(. < 0, NA_real_, .)))

# Remove metabolites with excessive missingness
nobs <- lipid_data %>%
  select(-c(1:4)) %>%
  describe() %>%
  select(n) %>%
  mutate(prop.missing = 1 - n / nrow(lipid_data))

nobs.high <- nobs %>%
  filter(prop.missing > 0.2) %>%
  arrange(-prop.missing) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ID")

ID.high.missing <- nobs.high %>% pull(ID)
lipid_data <- lipid_data %>% select(-all_of(ID.high.missing))
rm(nobs, nobs.high, ID.high.missing)

# 2B. Brain data
volume_data <- fread(file.path(wd, ado2datadir, "TNT_ADO2_May2026_volumes.csv"))
WM_vol_vars <- c("WMVolume", "FrontalRightWMVolume", "ParietalRightWMVolume",
                 "TemporalRightWMVolume", "OccipitalRightWMVolume",
                 "FrontalLeftWMVolume", "ParietalLeftWMVolume", "TemporalLeftWMVolume",
                 "OccipitalLeftWMVolume")
WM_vol_vars <- str_remove(WM_vol_vars, "Volume")

volume_data <- volume_data %>%
  mutate(lobar.volume = volume_data %>%
           select(all_of(WM_vol_vars[-1])) %>%
           data.frame() %>%
           apply(1, sum))

T1wSI_data <- fread(file.path(wd, ado2datadir, "TNT_ADO2_May2026_intensities_t1w.csv"))
MTR_data <- fread(file.path(wd, ado2datadir, "TNT_ADO2_May2026_intensities_mtr.csv"))

# 2C. Create analytic dataset for ADO2
analdat <- lipid_data %>%
  select(uniqueID) %>%
  left_join(covariate_data) %>%
  left_join(volume_data %>%
              mutate(uniqueID = str_replace(subject_id, "fam", "FAM")) %>%
              select(uniqueID, lobar.volume.WM = lobar.volume, ICV)) %>%
  left_join(T1wSI_data %>%
              mutate(uniqueID = str_replace(subject_id, "fam", "FAM")) %>%
              select(uniqueID, T1wSI.WM = WM)) %>%
  left_join(MTR_data %>%
              mutate(uniqueID = str_replace(subject_id, "fam", "FAM")) %>%
              select(uniqueID, MTR.WM = WM))

analdat %>% describe %>% select(n)
analdat_ADO2 <- analdat
ID_in_ADO2 <- analdat$uniqueID

# 3. ADO1 data import and cleaning ----
# 3A. Brain metrics
load("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/data/sys_d_merge_2026-04-01.Rdata")

d_merge <- d_merge %>%
  rename(Volume.WM = Volume.MW) %>%
  select(-Volume.WM.adjAgeSexICV) %>%
  mutate(
    Volume.WM.adjAgeSexICV = resid(lmer(lobar.vol.WM.adjICV ~ (age.c + age.c2) * Sex + (1 | fam_id),
                                        data = d_merge, na.action = na.exclude)),
    T1wSI.WM.adjAgeSex = resid(lmer(T1wSI.WM ~ (age.c + age.c2) * Sex + (1 | fam_id),
                                    data = d_merge, na.action = na.exclude)),
    MTR.WM.adjAgeSex = resid(lmer(LobarWM_Z ~ (age.c + age.c2) * Sex + (1 | fam_id),
                                  data = d_merge, na.action = na.exclude))
  )

brain_ado <- d_merge %>%
  select(uniqueID,
         Volume.WM = lobar.vol.WM,
         Volume.WM.adjICV = lobar.vol.WM.adjICV,
         Volume.WM.adjAgeSexICV,
         T1wSI.WM,
         T1wSI.WM.adjAgeSex,
         MTR.WM = LobarWM_Z,
         MTR.WM.adjAgeSex,
         ICV = ICV_TNT2012,
         icv.c,
         icv.c2) %>%
  mutate(icv.c = icv.c[, 1],
         icv.c2 = icv.c^2)

# 3B. Load ADO1 lipid data and metadata
res_WM1 <- fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/TabS1_2026-03-17.txt")
res_WM1 <- res_WM1 %>% filter(metaboliteid != "GlycA")

metabodata_ado <- fread("data/metabodata_ados.txt")
names(metabodata_ado) <- str_remove(names(metabodata_ado), "_nmol/mL")
metabodata_ado <- metabodata_ado %>%
  select(all_of(c("uniqueID", res_WM1$metaboliteid)))

lipid_ado1 <- fread("~/Library/CloudStorage/OneDrive-Personal/00_Work/01_Grants/2022/NIH2022_metabolon_ApoE/dataFromJaveria/ADO1-Lipidomic_june2021_species_conc.txt")
names(lipid_ado1) <- str_remove(names(lipid_ado1), "ADO1-Lipidomic_june2021_species_conc.")
lipid_ado1 <- lipid_ado1 %>% mutate(across(-c(1:3), as.numeric))
lipid_ado1 <- lipid_ado1 %>% mutate(across(-c(1:3), ~ if_else(. < 0, NA_real_, .)))

# Remove low-quality or too-missing metabolites
lipid_ado1 <- lipid_ado1 %>% filter(!is.na(`CER(22:0)`))
rm_lipids <- lipid_ado1 %>%
  describe() %>%
  select(n) %>%
  mutate(nmiss = 995 - n,
         pct.miss = nmiss / 995) %>%
  filter(pct.miss > 0.2) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "metaboliteid")

lipid_ado1 <- lipid_ado1 %>% select(-all_of(rm_lipids %>% pull(metaboliteid)))
metabodata_ado <- lipid_ado1 %>% select(-sub_id, -fam_id)

# 3C. Covariate data for ADO1
covdata <- fread("data/covdata_ados.txt")

# 3D. Adjust metabolites for age and sex
mnames <- names(metabodata_ado)[-1]
metabo.adjAgeSex <- metabodata_ado %>% select(uniqueID)

for (mi in seq_along(mnames)) {
  mname <- mnames[mi]
  analdat_tmp <- covdata %>%
    right_join(metabodata_ado %>% select(uniqueID, all_of(mname)), by = "uniqueID")
  analdat_tmp[["x.rnt"]] <- inormal(analdat_tmp[[mname]])
  
  fit.metabo <- lmer(x.rnt ~ (age.c + age.c2) * Sex + (1 | fam_id),
                     data = analdat_tmp, na.action = na.exclude)
  
  metabo.adjAgeSex <- data.table(metabo.adjAgeSex, m.resid = resid(fit.metabo))
  names(metabo.adjAgeSex)[mi + 1] <- mname
}

# 4. ADO1 association testing ----
ynames <- c("Volume.WM.adjAgeSexICV", "T1wSI.WM.adjAgeSex", "MTR.WM.adjAgeSex")
mnames <- names(metabodata_ado)[-1]

res <- list()
for (yi in seq_along(ynames)) {
  yname <- ynames[yi]
  
  for (mi in seq_along(mnames)) {
    mname <- mnames[mi]
    
    dat <- covdata %>%
      left_join(brain_ado %>% select(uniqueID, all_of(yname)), by = "uniqueID") %>%
      left_join(metabo.adjAgeSex %>% select(uniqueID, all_of(mname)), by = "uniqueID")
    
    dat[["y.rnt"]] <- inormal(dat[[yname]])
    dat[["m"]] <- inormal(dat[[mname]])
    
    fit <- lm(y.rnt ~ m + (age.c + age.c2), data = dat)
    
    resi <- c(
      trait = yname,
      metaboliteid = mname,
      summary(fit)$coef["m", -3],
      n = nobs(fit)
    )
    res[[length(res) + 1]] <- resi
  }
}

res <- data.table::rbindlist(lapply(res, as.list))
res <- res %>%
  mutate(across(-c(1:2), as.numeric)) %>%
  rename(beta = Estimate, se = `Std. Error`, pvalue = `Pr(>|t|)`) %>%
  mutate(
    fdr.p2 = p.adjust(pvalue, method = "fdr"),
    fdr.signif = ifelse(fdr.p2 < 0.05, "p.fdr<0.05", "p.fdr>=0.05"),
    fdr.signif = factor(fdr.signif, levels = c("p.fdr<0.05", "p.fdr>=0.05"))
  )

res$trait <- str_remove(res$trait, ".adjAgeSexICV")
res$trait <- str_remove(res$trait, ".adjAgeSex")
res$trait <- factor(res$trait, levels = c("Volume.WM", "T1wSI.WM", "MTR.WM"))

# Write ADO1 results if needed
# write_tsv(res, "assoc_WM_lipids_ADO1_present_in_ADO2.out")

# 5. ADO2 vs ADO1 comparison ----
# 5A. Read ADO2 results
res_ADO1 <- res_WM1
rm(res_WM1)

res_ADO2 <- fread("ADO2_assoc_res_WM_vs_lipids_2026-05-08.txt")

res_ADO2$metaboliteid2 <- harmonize_metabolite_ids(res_ADO2$metaboliteid)
res_ADO2 <- res_ADO2 %>%
  mutate(metaboliteid = metaboliteid2) %>%
  select(-metaboliteid2)

# Harmonize ADO1 names to match ADO2
res1 <- res
res1$metaboliteid2 <- res1$metaboliteid
res1$metaboliteid2 <- harmonize_metabolite_ids(res1$metaboliteid2)

res1 <- res1 %>%
  filter(metaboliteid2 %in% res_ADO2$metaboliteid) %>%
  mutate(metaboliteid = metaboliteid2) %>%
  select(-metaboliteid2)

res1$trait <- str_remove(res1$trait, ".WM")
res2 <- res_ADO2 %>% filter(metaboliteid %in% res1$metaboliteid)
res12 <- rbind(
  res1 %>% mutate(group = "ados12"),
  res2 %>% mutate(group = "ados2")
)

res12 <- res12 %>%
  mutate(
    fdr.p2 = p.adjust(pvalue, method = "fdr"),
    fdr.signif = ifelse(fdr.p2 < 0.05, "p.fdr<0.05", "p.fdr>=0.05")
  )

# 5B. Significance threshold across combined dataset
signif.level12 <- res12 %>%
  filter(fdr.p2 < 0.05) %>%
  pull(pvalue) %>%
  max()

# 6. Volcano plots for ADO1 + ADO2 combined ----
plot.point.size <- 1.5
plot_traits <- c("Volume", "T1wSI", "MTR")

panel_ado12 <- list()
for (trait in plot_traits) {
  panel_ado12[[trait]] <- generate_volcano_plot(
    df_plot = res12 %>% filter(group == "ados12", trait == trait) %>% mutate(generation = "ados12"),
    trait_plot = trait,
    generation_plot = "ados12",
    signif.level = signif.level12,
    point.size = plot.point.size
  )
}

p_volcano_ado12 <- panel_ado12[["Volume"]] +
  panel_ado12[["T1wSI"]] +
  panel_ado12[["MTR"]] +
  plot_layout(ncol = 3, guides = "collect", axes = "collect") &
  coord_cartesian(
    xlim = max(abs(res1$beta), na.rm = TRUE) * 1.015 * c(-1, 1),
    ylim = c(0, max(-log10(res_ADO1 %>% select(contains("pvalue")) %>% select(contains("ados"))), na.rm = TRUE))
  ) &
  xlab(NULL) & ylab(NULL)

p_volcano_ado12
ggsave("results/[manuscript]_Fig1_volcalnoplots_ADO12_may2026_822_lipids_955_ADO1.png",
       plot = p_volcano_ado12,
       width = 10.99 * 2 * 1.45,
       height = 14.99 * 1.45 * 0.5,
       units = "cm",
       dpi = 300)

panel_ado2 <- list()
for (trait in plot_traits) {
  panel_ado2[[trait]] <- generate_volcano_plot(
    df_plot = res12 %>% filter(group == "ados2", trait == trait) %>% mutate(generation = "ados2"),
    trait_plot = trait,
    generation_plot = "ados2",
    signif.level = signif.level12,
    point.size = plot.point.size
  )
}

p_volcano_ado2 <- panel_ado2[["Volume"]] +
  panel_ado2[["T1wSI"]] +
  panel_ado2[["MTR"]] +
  plot_layout(ncol = 3, guides = "collect", axes = "collect") &
  coord_cartesian(
    xlim = max(abs(res1$beta), na.rm = TRUE) * 1.015 * c(-1, 1),
    ylim = c(0, max(-log10(res_ADO1 %>% select(contains("pvalue")) %>% select(contains("ados"))), na.rm = TRUE))
  ) &
  xlab(NULL) & ylab(NULL)

p_volcano_ado2
ggsave("results/[manuscript]_Fig1_volcalnoplots_ADO2_may2026.png",
       plot = p_volcano_ado2,
       width = 10.99 * 2 * 1.45,
       height = 14.99 * 1.45 * 0.5,
       units = "cm",
       dpi = 300)

# 7. Additional checks ----
res_ADO2_tmp <- res12 %>% filter(group == "ados2")
res_ADO2_tmp <- res_ADO2_tmp %>% mutate(fdr.p.ado2 = p.adjust(pvalue, method = "fdr"))
res_ADO2_tmp %>% arrange(pvalue)
