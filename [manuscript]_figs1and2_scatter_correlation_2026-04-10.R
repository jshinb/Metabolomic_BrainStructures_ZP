# Generate Table S1
# Created: Apr. 1, 2026
# Created by: Jean Shin
# Last modified: Apr. 1, 2026

# Inputs:
# - SPS data
# - SYS data
# - Age
# - Sex
# - BMI
# - % overweight/obesity
# - Type 2 diabetes
# - Hypertension
# - Elevated total cholesterol
# - Elevated TAGs
# - Elevated LDL-cholesterol
# - Glomerular filtration rate (GFR)
#
# Note: Some variables may need to be downloaded/derived separately.

# 0. Setup ----
library(tidyverse)
library(data.table)
library(readr)
library(ggpubr)
library(psych)
library(lmerTest)

rename <- dplyr::rename
select <- dplyr::select
describe <- psych::describe

inormal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# Working directory
wd <- "/Users/jshin/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP"
setwd(wd)

# 1. Load data ----
load("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/data/sys_d_merge_2026-04-01.Rdata")
d_merge <- d_merge %>% rename(Volume.WM = Volume.MW)

metaboinfo_ado <- fread("data/metabodata_info_ados.txt")
metaboinfo_adult <- fread("data/metabodata_info_adults.txt")

metaboinfo_ado <- metaboinfo_ado %>%
  mutate(metabolite_id = str_remove(metabolite_id, "_nmol/mL")) %>%
  filter(!str_detect(metabolite_id, "_mol%"))

metaboinfo_adult <- metaboinfo_adult %>%
  mutate(metabolite_id = str_remove(metabolite_id, "_nmol/mL")) %>%
  filter(!str_detect(metabolite_id, "_mol%"))

# 2. Metabolite metadata summaries ----
# Class counts in ADO
metaboinfo_ado %>%
  filter(platform == "CLP") %>%
  filter(str_detect(metabolite_id, "[:]")) %>%
  filter(!str_detect(metabolite_id, "\\[FA")) %>%
  pull(class) %>%
  table

# Class counts in adults
metaboinfo_adult %>%
  filter(platform == "CLP") %>%
  filter(str_detect(metabolite_id, "[:]")) %>%
  filter(!str_detect(metabolite_id, "\\[FA")) %>%
  pull(class) %>%
  table %>%
  sum

# Explore class-specific IDs
metaboinfo_ado %>% filter(class == "CE")
metaboinfo_ado %>% filter(str_detect(metabolite_id, "[:]"))
metaboinfo_ado %>% filter(class == "PI") %>% select(metabolite_id)

# Load ADO1 WM association results
res_WM1 <- fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/TabS1_2026-03-17.txt")
res_WM1 <- res_WM1 %>% filter(metaboliteid != "GlycA")

# 3. Filter metabolites by significance pattern ----
# ADO significance across any trait
ind_ado <- res_WM1$fdr.p2_WM.T1wSI_ados < 0.05
ind_ado <- ind_ado | (res_WM1$fdr.p2_WM.MTR_ados < 0.05)
ind_ado <- ind_ado | (res_WM1$fdr.p2_WM.Vol_ados < 0.05)
print(sum(ind_ado)) # 711 after removing GlycA

# ADO: all three traits significant
ind_ado1 <- res_WM1$fdr.p2_WM.T1wSI_ados < 0.05
ind_ado1 <- ind_ado1 & (res_WM1$fdr.p2_WM.MTR_ados < 0.05)
ind_ado1 <- ind_ado1 & (res_WM1$fdr.p2_WM.Vol_ados < 0.05)
print(sum(ind_ado1)) # 711 after removing GlycA

# ADO: relaxed threshold 0.1
ind_ado2 <- res_WM1$fdr.p2_WM.T1wSI_ados < 0.1
ind_ado2 <- ind_ado2 & (res_WM1$fdr.p2_WM.MTR_ados < 0.1)
ind_ado2 <- ind_ado2 & (res_WM1$fdr.p2_WM.Vol_ados < 0.1)
print(sum(ind_ado2)) # 564

# Adult significance across traits
ind_adult <- res_WM1$fdr.p2_WM.T1wSI_adults < 0.05
ind_adult <- ind_adult | (res_WM1$fdr.p2_WM.MTR_adults < 0.05)
ind_adult <- ind_adult | (res_WM1$fdr.p2_WM.Vol_adults < 0.05)
print(sum(ind_adult)) # 73

# Adult relaxed threshold 0.1
ind_adult2 <- res_WM1$fdr.p2_WM.T1wSI_adults < 0.1
ind_adult2 <- ind_adult2 & (res_WM1$fdr.p2_WM.MTR_adults < 0.1)
ind_adult2 <- ind_adult2 & (res_WM1$fdr.p2_WM.Vol_adults < 0.1)
print(sum(ind_adult2)) # 0

# Correlations among effect sizes for significant ADO results
res_WM1 %>%
  filter(ind_ado) %>%
  select(beta_WM.T1wSI_ados, beta_WM.MTR_ados, beta_WM.Vol_ados) %>%
  cor(use = "p")

res_WM1 %>%
  filter(ind_ado2) %>%
  select(beta_WM.T1wSI_ados, beta_WM.MTR_ados, beta_WM.Vol_ados) %>%
  cor(use = "p")

res_WM1 %>%
  filter(ind_adult) %>%
  select(beta_WM.T1wSI_adults, beta_WM.MTR_adults, beta_WM.Vol_adults) %>%
  cor(use = "p")

# 4. Metabolite class processing ----
lipid.class <- c("NeutL", "GPC", "SL", "NonL")
class2.levels <- c(
  "Neutral lipids",
  "Glycerophospholipids",
  "Sphingolipids",
  "Non-lipid metabolites"
)

res_WM1 <- res_WM1 %>% rename(metabolite_id = metaboliteid)

metaboinfo_ado %>%
  filter(metabolite_id %in% res_WM1$metabolite_id) %>%
  pull(class) %>%
  table %>%
  data.frame()

# Quick inspection of metabolite distributions
fread("data/metabodata_ados.txt") %>%
  select("L_VLDL_PL_pct") %>%
  describe()

fread("data/metabodata_ados.txt") %>%
  describe() %>%
  select(min) %>%
  arrange(min)

hist((fread("data/metabodata_ados.txt") %>% pull(`Acetoacetate`)))
hist(log(fread("data/metabodata_adults.txt") %>% pull(`TAG36:0FA12:0_nmol/mL`)))

# 5. Load brain phenotypes ----
# Volume: adults
Volume_data <- fread(file.path(wd, "data", "lobarVolData_WM_GM_adults_2025-05-26.tsv"))

# T1wSI: adults
T1wSI_data <- fread(file.path(wd, "data", "braindata_adults.txt"))

# MTR: adults
MTR_data <- fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/jean_shin_dataanalysis_completespssubjects_final_20251125_132248.txt")
names(MTR_data) <- str_remove(names(MTR_data), "archived.SPS-MRI-MTRpipe-MRI-Apr2014.")
MTR_data <- MTR_data %>%
  mutate(across(-c(1:3), as.numeric)) %>%
  mutate(across(c(wm, gm), ~ ifelse(.x < 0, NA, .x)))

# 6. ADO phenotype adjustments ----
base.size <- 12
plot.alpha <- 0.2

# Refit WM phenotypes to adjust for age/sex/ICV
d_merge <- as.data.frame(d_merge)
rownames(d_merge) <- d_merge$uniqueID

# For adult-like residualized age-sex adjustments
d_merge <- d_merge %>%
  mutate(
    Volume.WM.adjAgeSexICV = resid(
      lmer(lobar.vol.WM.adjICV ~ (age.c + age.c2) * Sex + (1 | fam_id),
           data = d_merge, na.action = na.exclude)
    ),
    T1wSI.WM.adjAgeSex = resid(
      lmer(T1wSI.WM ~ (age.c + age.c2) * Sex + (1 | fam_id),
           data = d_merge, na.action = na.exclude)
    ),
    MTR.WM.adjAgeSex = resid(
      lmer(LobarWM_Z ~ (age.c + age.c2) * Sex + (1 | fam_id),
           data = d_merge, na.action = na.exclude)
    )
  )

# Use sex-specific intercept adjustment for narrative plots
fit_ado_volume <- lmer(lobar.vol.WM ~ (age.c + age.c2 + icv.c + icv.c2) + (1 | fam_id),
                       data = d_merge %>% mutate(lobar.vol.WM = inormal(lobar.vol.WM)) %>% filter(Sex == "F"),
                       na.action = na.exclude)
mean.F <- summary(fit_ado_volume)$coef[1, 1]

fit_ado_volume <- lmer(lobar.vol.WM ~ (age.c + age.c2 + icv.c + icv.c2) + (1 | fam_id),
                       data = d_merge %>% mutate(lobar.vol.WM = inormal(lobar.vol.WM)) %>% filter(Sex == "M"),
                       na.action = na.exclude)
mean.M <- summary(fit_ado_volume)$coef[1, 1]

d_merge <- d_merge %>%
  mutate(
    Volume.WM.adjAgeICV = ifelse(Sex == "F",
                                Volume.WM.adjAgeSexICV + mean.F,
                                Volume.WM.adjAgeSexICV + mean.M)
  )

# T1wSI adjustment
fit_ado_t1wsi <- lmer(T1wSI.WM ~ (age.c + age.c2) + (1 | fam_id),
                     data = d_merge %>% filter(Sex == "F"),
                     na.action = na.exclude)
mean.F <- summary(fit_ado_t1wsi)$coef[1, 1]

fit_ado_t1wsi <- lmer(T1wSI.WM ~ (age.c + age.c2) + (1 | fam_id),
                     data = d_merge %>% filter(Sex == "M"),
                     na.action = na.exclude)
mean.M <- summary(fit_ado_t1wsi)$coef[1, 1]

d_merge <- d_merge %>%
  mutate(
    T1wSI.WM.adjAge = ifelse(Sex == "F",
                            T1wSI.WM.adjAgeSex + mean.F,
                            T1wSI.WM.adjAgeSex + mean.M)
  )

# MTR adjustment
d_merge <- d_merge %>% rename(MTR.WM = LobarWM_Z)

fit_ado_mtr <- lmer(MTR.WM ~ (age.c + age.c2) + (1 | fam_id),
                    data = d_merge %>% filter(Sex == "F"),
                    na.action = na.exclude)
mean.F <- summary(fit_ado_mtr)$coef[1, 1]

fit_ado_mtr <- lmer(MTR.WM ~ (age.c + age.c2) + (1 | fam_id),
                    data = d_merge %>% filter(Sex == "M"),
                    na.action = na.exclude)
mean.M <- summary(fit_ado_mtr)$coef[1, 1]

d_merge <- d_merge %>%
  mutate(
    MTR.WM.adjAge = ifelse(Sex == "F",
                          MTR.WM.adjAgeSex + mean.F,
                          MTR.WM.adjAgeSex + mean.M)
  )

# 7. Pairwise trait correlations ----
# Sex-combined data
wm_corr <- d_merge %>%
  select(Volume.WM.adjAgeSexICV, T1wSI.WM.adjAgeSex, MTR.WM.adjAgeSex) %>%
  cor(use = "p")

wm_corr

# Plot helper
make_pairwise_plot <- function(df, xvar, yvar, color = "Sex") {
  df %>%
    select(Sex, all_of(c(xvar, yvar))) %>%
    mutate(across(contains("WM"), inormal)) %>%
    ggplot(aes(x = .data[[xvar]], y = .data[[yvar]], color = Sex)) +
    geom_point(alpha = plot.alpha) +
    geom_smooth(method = "gam") +
    stat_cor()
}

p_1 <- make_pairwise_plot(d_merge, "Volume.WM.adjAgeSexICV", "T1wSI.WM.adjAgeSex")
p_2 <- make_pairwise_plot(d_merge, "Volume.WM.adjAgeSexICV", "MTR.WM.adjAgeSex")
p_3 <- make_pairwise_plot(d_merge, "T1wSI.WM.adjAgeSex", "MTR.WM.adjAgeSex")

p_sex_combined <- (p_1 + p_2 + p_3) +
  plot_layout(ncol = 3, guides = "collect") &
  geom_point(alpha = plot.alpha)

# 8. Adult phenotype preparation ----
analdat_adult <- Volume_data %>%
  left_join(fread("data/covdata_adults_v2.txt")) %>%
  rename(age = age_years) %>%
  mutate(sex = ifelse(Sex == 0, "M", "F")) %>%
  select(-Sex) %>%
  rename(Sex = sex) %>%
  mutate(age.c = scale(age)[, 1] / 2,
         age.c2 = age.c^2)

# WM volume in adults
fit.adult <- lm(lobar.vol.WM.adjICV ~ Sex * (age.c + age.c2 + on_lipid_low_med),
               data = analdat_adult, na.action = na.exclude)

analdat_adult <- analdat_adult %>%
  mutate(
    Volume.WM.adjICVSex = resid(fit.adult),
    Volume.WM.adjAgeSex = resid(fit.adult)
  )

# T1wSI in adults
analdat_adult <- analdat_adult %>% left_join(T1wSI_data)

fit.adult_t1wsi <- lm(NormWM ~ Sex, data = analdat_adult, na.action = na.exclude)
analdat_adult <- analdat_adult %>%
  mutate(T1wSI.WM.adjSex = resid(fit.adult_t1wsi),
         T1wSI.WM.adjAgeSex = resid(fit.adult_t1wsi))

# MTR in adults
analdat_adult <- analdat_adult %>% left_join(MTR_data) %>% rename(MTR.WM = wm)

fit.adult_mtr <- lm(MTR.WM ~ Sex, data = analdat_adult, na.action = na.exclude)
fit.adult_mtr2 <- lm(MTR.WM ~ Sex * (age.c + age.c2 + on_lipid_low_med),
                     data = analdat_adult, na.action = na.exclude)

analdat_adult <- analdat_adult %>%
  mutate(
    MTR.WM.adjSex = resid(fit.adult_mtr),
    MTR.WM.adjAgeSex = resid(fit.adult_mtr2)
  )

# 9. ADO phenotype preparation ----
# WM volume in ADO
fit.ado <- lmer(lobar.vol.WM.adjICV ~ Sex + (1 | fam_id),
                data = d_merge, na.action = na.exclude)
d_merge <- d_merge %>% mutate(Volume.WM.adjICVSex = resid(fit.ado))

# WM T1wSI in ADO
fit.ado_t1w <- lmer(T1wSI.WM ~ Sex + (1 | fam_id),
                    data = d_merge, na.action = na.exclude)
d_merge <- d_merge %>% mutate(T1wSI.WM.adjSex = resid(fit.ado_t1w))

# WM MTR in ADO
fit.ado_mtr <- lm(MTR.WM ~ Sex, data = d_merge, na.action = na.exclude)
d_merge <- d_merge %>% mutate(MTR.WM.adjSex = resid(fit.ado_mtr))

# 10. Figure S1: age vs white matter plots ----
p_Volume_age_ado <- d_merge %>%
  mutate(y = inormal(lobar.vol.WM)) %>%
  ggplot(aes(x = age, y = y, color = Sex)) +
  geom_point(alpha = plot.alpha) +
  geom_smooth(method = "gam", show.legend = FALSE) +
  stat_cor()

p_Volume_age_adult <- analdat_adult %>%
  mutate(y = inormal(lobar.vol.WM.adjICV)) %>%
  ggplot(aes(x = age, y = y, color = Sex)) +
  geom_point(alpha = plot.alpha) +
  geom_smooth(method = "gam", show.legend = FALSE) +
  stat_cor()

p_T1wSI_age_ado <- d_merge %>%
  mutate(y = inormal(T1wSI.WM.adjSex)) %>%
  ggplot(aes(x = age, y = y, color = Sex)) +
  geom_point(alpha = plot.alpha) +
  geom_smooth(method = "gam", show.legend = FALSE) +
  stat_cor()

p_T1wSI_age_adult <- analdat_adult %>%
  mutate(y = inormal(T1wSI.WM.adjSex)) %>%
  ggplot(aes(x = age, y = y, color = Sex)) +
  geom_point(alpha = plot.alpha) +
  geom_smooth(method = "gam", show.legend = FALSE) +
  stat_cor()

p_MTR_age_ado <- d_merge %>%
  mutate(y = inormal(MTR.WM.adjSex)) %>%
  ggplot(aes(x = age, y = y, color = Sex)) +
  geom_point(alpha = plot.alpha) +
  geom_smooth(method = "gam", show.legend = FALSE) +
  stat_cor()

p_MTR_age_adult <- analdat_adult %>%
  mutate(y = inormal(MTR.WM.adjSex)) %>%
  ggplot(aes(x = age, y = y, color = Sex)) +
  geom_point(alpha = plot.alpha) +
  geom_smooth(method = "gam", show.legend = FALSE) +
  stat_cor()

p_S1 <- (
  p_Volume_age_ado + p_Volume_age_adult +
  p_T1wSI_age_ado + p_T1wSI_age_adult +
  p_MTR_age_ado + p_MTR_age_adult
) +
  plot_layout(ncol = 2, byrow = TRUE, guides = "collect") &
  ggtitle(NULL) &
  ylab(NULL) &
  xlab(NULL) &
  theme_bw(base_size = base.size) &
  theme(
    panel.background = element_rect(fill = scales::alpha("#f8edeb", 0.2), size = 0.5, linetype = "solid"),
    panel.border = element_rect(colour = "grey80"),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey95"),
    panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey95")
  ) &
  coord_cartesian(ylim = c(-3.5, 4.25)) &
  stat_cor(label.y = 3.8)

p_S1

# Save figure
ggsave(
  "../Metabolomic_BrainStructures_ZP_clean/manuscript/figureS1.age_WM.png",
  plot = p_S1,
  height = 5.5 * 1.35,
  width = 4 * 1.35,
  units = "in",
  dpi = 300
)

# 11. Final checks ----
# Correlations in ADO after adjustments
d_merge %>%
  select(Volume.WM.adjAgeSexICV, T1wSI.WM.adjAgeSex, MTR.WM.adjAgeSex) %>%
  cor(use = "p")

# Stratified sex correlations
d_merge %>%
  filter(Sex == "F") %>%
  select(Volume.WM.adjAgeSexICV, T1wSI.WM.adjAgeSex, MTR.WM.adjAgeSex) %>%
  mutate(across(everything(), inormal)) %>%
  cor(use = "p")

d_merge %>%
  filter(Sex == "M") %>%
  select(Volume.WM.adjAgeSexICV, T1wSI.WM.adjAgeSex, MTR.WM.adjAgeSex) %>%
  mutate(across(everything(), inormal)) %>%
  cor(use = "p")
