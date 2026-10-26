#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

# ADO1 data import and cleaning ----
## A. Brain metrics ----
load("data/sys_d_merge_2026-04-01.Rdata")

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

## B. Load ADO1 lipid data and metadata ----
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

## C. Covariate data for ADO1 ----
covdata <- fread("data/covdata_ados.txt")

## D. Adjust metabolites for age and sex ----
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

# SPS data import and cleaning ----
## A. brain metrics ----
MTR_data = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/jean_shin_dataanalysis_completespssubjects_final_20251125_132248.txt")
head(MTR_data) #n=403
names(MTR_data) = str_remove(names(MTR_data),"archived.SPS-MRI-MTRpipe-MRI-Apr2014.")
MTR_data = MTR_data %>% mutate_at(-c(1:3),as.numeric)
MTR_data = MTR_data %>% mutate(across(c(wm, gm), ~ ifelse(.x < 0, NA, .x)))#n=398

ICV_data= fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/jean_shin_dataanalysis_completespssubjects_final_20251003_185216.txt")
names(ICV_data) = str_remove(names(ICV_data),"archived.SPS-MRI-TNTpipe-Jan2016 --DO NOT USE.")
names(ICV_data) = str_remove(names(ICV_data),"archived.SPS-MRI-TNTpipe-MRI-Jan2014.")
names(ICV_data) = str_remove(names(ICV_data),"SPS-MRI-Freesurfer5.3.0-MRI-Aug2018.")
names(ICV_data) = str_remove(names(ICV_data),"SPS-MRI-Freesurfer5.3.0-MRI-Sept2013.")
names(ICV_data)[5] <- "QC.TNTpipe"
ICV_data = data.frame(ICV_data)
ICV_data[,-c(1:3)] = apply(ICV_data[,-c(1:3)],2,as.numeric)
ICV_data %>% describe %>% select(1:4)
ICV_data = ICV_data %>% rename(ICV_TNT2014=icv,ICV_TNT2016=ICV)
describe(ICV_data) %>% select(n)
#                          ICV_FS ICV_TNT2016 ICV_TNT2014 ICV_usedfor_FSAD_FSBA
# ICV_FS                1.0000000   0.7500408   0.7376012             1.0000000
# ICV_TNT2016           0.7500408   1.0000000   0.9914927             0.7445772
# ICV_TNT2014           0.7376012   0.9914927   1.0000000             0.7287898
# ICV_usedfor_FSAD_FSBA 1.0000000   0.7445772   0.7287898             1.0000000

brain_adult = fread(file.path(wd,"data","lobarVolData_WM_GM_adults_2025-05-26.tsv")) %>% 
  select(uniqueID,lobar.vol.WM,lobar.vol.WM.adjICV) %>%
  left_join(fread(file.path(wd,"data","braindata_adults.txt")) %>% 
              select(uniqueID,T1wSI.WM = NormWM,ICV_FS = ICV)) %>%
  left_join(MTR_data %>% select(uniqueID,wm) %>% rename(MTR.WM = wm)) %>%#n=398
  left_join(ICV_data %>% select(uniqueID,ICV_TNT2016,ICV_TNT2014,ICV_usedfor_FSAD_FSBA))

brain_adult = brain_adult %>% 
  select(uniqueID, Volume.WM = lobar.vol.WM,Volume.WM.adjICV = lobar.vol.WM.adjICV,
         T1wSI.WM, MTR.WM, ICV=ICV_TNT2014) %>% 
  mutate(icv.c = scale(ICV)[,1]/2) %>%
  mutate(icv.c2=icv.c^2)

# outliers
# get age- and sex-adjusted in adults
analdat = brain_adult %>% left_join(covdata_adult)
fit_vol = lm(Volume.WM.adjICV ~ Sex*(age.c+age.c2), data=analdat, na.action=na.exclude)
summary(fit_vol)
fit_T1wSI = lmer(T1wSI.WM ~ Sex*(age.c+age.c2)+(1|fam_id), data=analdat, na.action=na.exclude)
fit_MTR = lmer(MTR.WM ~ Sex*(age.c+age.c2)+(1|fam_id), data=analdat, na.action=na.exclude)
op = par(mfrow=c(2,2));plot(fit_vol);par(op)
op = par(mfrow=c(2,2));plot(fit_T1wSI);par(op)
op = par(mfrow=c(2,2));plot(fit_MTR);par(op)

brain_adult = brain_adult %>% 
  mutate(Volume.WM.adjAgeSexICV=resid(fit_vol),
         T1wSI.WM.adjAgeSex=resid(fit_T1wSI),
         MTR.WM.adjAgeSex=resid(fit_MTR))

### age- and sex-adjusted WM phenotypes in adults
p_adj = brain_adult %>% select(contains(".WM")) %>% 
  cor(use="p") %>% ggcorrplot::ggcorrplot(hc.order = T,lab=T)
p_adj = p_adj + ggtitle('untransformed (adjusted) phenotypes in adults')

### age- and sex-adjusted WM phenotypes and INT
p_adj_INT = brain_adult %>% select(contains(".WM")) %>% mutate_all(inormal) %>%
  cor(use="p") %>% ggcorrplot::ggcorrplot(hc.order = T,lab=T)
p_adj_INT = p_adj_INT + ggtitle('INT-transformed phenotypes in adults')

p_adj + p_adj_INT

## B. lipid and metabolite variables ----
metabodata_adult = fread('data/metabodata_adults.txt') 
names(metabodata_adult) = str_remove(names(metabodata_adult),"_nmol/mL")

metabodata_all = metabodata_ado %>% select(all_of(c('uniqueID',res_WM1$metaboliteid))) %>% 
  bind_rows(
    metabodata_adult %>% select(all_of(c('uniqueID',res_WM1$metaboliteid))) 
  )

## C. covariates ----
covdata_adult = fread('data/covdata_adults_v2.txt') %>%
  rename(age=age_years,sex=Sex) %>%
  mutate(age.c = scale(age)[,1]/2) %>%
  mutate(age.c2 = age.c^2, Sex=ifelse(sex==0,"M","F")) %>%
  mutate(fam_id = str_split(uniqueID,"_",simplify = T)[,1]) %>%
  mutate(group="adult") %>%
  select(all_of(c(names(covdata_ado)))) 

## D. merging ----
### get age- and sex-adjusted in adults
analdat = brain_adult %>% left_join(covdata_adult)
fit_vol = lm(Volume.WM.adjICV ~ Sex*(age.c+age.c2), data=analdat, na.action=na.exclude)
summary(fit_vol)
fit_T1wSI = lmer(T1wSI.WM ~ Sex*(age.c+age.c2)+(1|fam_id), data=analdat, na.action=na.exclude)
fit_MTR = lmer(MTR.WM ~ Sex*(age.c+age.c2)+(1|fam_id), data=analdat, na.action=na.exclude)
op = par(mfrow=c(2,2));plot(fit_vol);par(op)
op = par(mfrow=c(2,2));plot(fit_T1wSI);par(op)
op = par(mfrow=c(2,2));plot(fit_MTR);par(op)

brain_adult = brain_adult %>% 
  mutate(Volume.WM.adjAgeSexICV=resid(fit_vol),
         T1wSI.WM.adjAgeSex=resid(fit_T1wSI),
         MTR.WM.adjAgeSex=resid(fit_MTR))

### age- and sex-adjusted WM phenotypes in adults
p_adj = brain_adult %>% select(contains(".WM")) %>% 
  cor(use="p") %>% ggcorrplot::ggcorrplot(hc.order = T,lab=T)
p_adj = p_adj + ggtitle('untransformed (adjusted) phenotypes in adults')

### age- and sex-adjusted WM phenotypes and INT
p_adj_INT = brain_adult %>% select(contains(".WM")) %>% mutate_all(inormal) %>%
  cor(use="p") %>% ggcorrplot::ggcorrplot(hc.order = T,lab=T)
p_adj_INT = p_adj_INT + ggtitle('INT-transformed phenotypes in adults')

p_adj + p_adj_INT

## merge
### merge ado + adult data 
brain_all = brain_ado %>% 
  bind_rows(brain_adult %>% select(names(brain_ado)))
brain_all %>% describe
brain_all %>% left_join(covdata_all) %>% ggplot(aes(x=age, y=Volume.WM.adjAgeSexICV)) + geom_point()
brain_all %>% left_join(covdata_all) %>% 
  ggplot(aes(x=Volume.WM.adjAgeSexICV,fill=group,colour = group)) + 
  geom_histogram(alpha=0.5, position="identity")
describe(brain_all$Volume.WM.adjAgeSexICV)
brain_all %>% left_join(covdata_all) %>% 
  ggplot(aes(y=Volume.WM.adjAgeSexICV,x=age,fill=group,colour = group)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method="gam") + stat_cor() 

# Resampling ----
rfile='https://raw.githubusercontent.com/jshinb/Metabolomic_BrainStructures_ZP/refs/heads/main/%5Bmanuscript%5D_sample_size_matching_ADO1_SPS.R'
source(rfile)

# key metabolites from ADO1: for example top metabolites by p-value
# ADO1_dat = get_data(ado1) 
# SPS_dat = get_data (sps)

## extract the metabolites associated with all phenotypes in ADO1 ----
key_metabs <- [res_WM1$metaboliteid]
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
