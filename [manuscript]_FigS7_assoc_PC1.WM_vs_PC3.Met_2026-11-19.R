# ----------------------------@
# Supplementary Figure S2
# Created: 2026-10-27
# Modified: 2026-11-17
# ----------------------------@

# 0. Setup ----
wd <- "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP"
setwd(wd)

library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
select = dplyr::select
describe = psych::describe

source(file.path(
  wd,
  "99_scripts",
  "[function]_define_fit_function_2025-01-09.R"
))
d_pheno = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/results/used_for_manuscript/gwas_compVar1_unadjBMI_run_gwas_regenie_20260324_zp/pheno.input")
head(d_pheno)
covdataa_ado = fread('data/covdata_ados.txt')

d_merge = covdataa_ado %>% 
  mutate(IID = as.integer(str_split( uniqueID,"_",simplify = T  )[,2])) %>%
  left_join(d_pheno)
(d_merge %>% ggplot(aes(y=inormal(WM.PC1),x=BMI)) + geom_point() +
  geom_smooth(method = 'gam') + ggpubr::stat_cor()) +
(d_merge %>% ggplot(aes(y=inormal(PC3met),x=BMI)) + geom_point() +
  geom_smooth(method = 'gam') + ggpubr::stat_cor())

# BMI strongly associated with WM.PC1 but weaker with MetPC3
d_merge %>% ggplot(aes(y=WM.PC1,x=BMI)) + geom_point() +
  geom_smooth(method = 'gam') + ggpubr::stat_cor()

library(lmerTest)
## PC1.WM vs. PC3.Met ----
fit.lmer.PC1WM_PC3Met = lmer(inormal(WM.PC1)~inormal(PC3met) + (BMI + bmi2)*Sex + (1|fam_id), 
     data=d_merge %>% mutate(bmi2=BMI^2/4) )
qqnorm(resid(fit.lmer.PC1WM_PC3Met));abline(0,1)
plot(fit.lmer.PC1WM_PC3Met)

fit.lmer.PC1WM_PC3Met0 = lmer(inormal(WM.PC1)~inormal(PC3met) + Sex + (1|fam_id), 
                             data=d_merge %>% mutate(bmi2=BMI^2/4) %>% filter(!is.na(BMI)) )
qqnorm(resid(fit.lmer.PC1WM_PC3Met0));abline(0,1)
plot(fit.lmer.PC1WM_PC3Met0)
rbind(
summary(fit.lmer.PC1WM_PC3Met0)$coef%>%
  data.table %>% select(-contains("value"),-contains("df")) %>%
  slice(2),
summary(fit.lmer.PC1WM_PC3Met)$coef %>% 
  data.table %>% select(-contains("value"),-contains("df")) %>%
  slice(2)
) %>% bind_cols(adj = c("unadjBMI","adjBMI"))

## PC1.WM vs. PC1.Met ----
fit.lmer.PC1WM_PC1Met = lmer(
  inormal(WM.PC1)~inormal(PC1met) + (BMI + bmi2)*Sex + (1|fam_id), 
                             data=d_merge %>% mutate(bmi2=BMI^2/4) )
qqnorm(resid(fit.lmer.PC1WM_PC1Met));abline(0,1)
plot(fit.lmer.PC1WM_PC1Met)

fit.lmer.PC1WM_PC1Met0 = lmer(
  inormal(WM.PC1)~inormal(PC1met) + Sex + (1|fam_id), 
  data=d_merge %>% mutate(bmi2=BMI^2/4) %>% filter(!is.na(BMI)) )
qqnorm(resid(fit.lmer.PC1WM_PC3Met0));abline(0,1)
plot(fit.lmer.PC1WM_PC3Met0)

rbind(
  summary(fit.lmer.PC1WM_PC1Met0)$coef%>%
    data.table %>% select(-contains("value"),-contains("df")) %>%
    slice(2),
  summary(fit.lmer.PC1WM_PC1Met)$coef %>% 
    data.table %>% select(-contains("value"),-contains("df")) %>%
    slice(2)
) %>% bind_cols(adj = c("unadjBMI","adjBMI"))
