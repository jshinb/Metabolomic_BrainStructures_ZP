#0. load libraries ----
rm(list=ls())
wd = "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP"
setwd(wd)

library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(psych)
library(patchwork)
library(ggrepel)
library(lme4)
library(coxme)
library(lmerTest)
library(outliers)
library(ggforestplot)
library(EnvStats)# outliers:
library(dplyr)
library(lubridate)
library(ggpp)
library(grid)
library(FactoMineR)
source(file.path(wd,'scripts','[function]_define_fit_function_2025-01-09.R'))
select = dplyr::select
describe = psych::describe


braindata = fread('data/braindata_ados.txt')

## covdata
covdata = fread(file.path(wd,"data","covdata_ados.txt"))

## metabolite and info data
metabodata_info = fread("data/metabodata_info_ados.txt")
head(metabodata_info)

## metabolite data
metabodata = fread(file.path(wd,"data","metabodata_ados.txt"))
metabodata = subset(
  metabodata,
  select=c('uniqueID',metabodata_info$metabolite_id)
)#
# checking nightigale:
d_NG = fread("~/Downloads/jean_shin_dataanalysis_completesyssubjects_rev1_20251203_111654.txt")
d_NG = d_NG %>% select(1:10)
names(d_NG) = str_remove(names(d_NG),"ADO1-Biomarkers-2016_pipeline.")
d_NG = data.frame(d_NG)
d_NG[,-c(1:3)] = apply(d_NG[,-c(1:3)],2,as.numeric)
d_NG[,-c(1:3)] %>% describe
# ados MTR data ----
source(file.path(wd,'scripts','[function]_define_fit_function_2025-01-09.R'))
MTR_data = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/jean_shin_dataanalysis_completesyssubjects_rev1_20251125_131931.txt")
names(MTR_data) = str_remove(names(MTR_data),'ADO1-MRI-MTRPipeAug2012.')
MTR_data = MTR_data %>% 
  mutate(LobarGM_Z = as.numeric(LobarGM_Z)) %>%
  mutate(LobarGM_Z = ifelse(LobarGM_Z==-666,NA,LobarGM_Z))
describe(MTR_data)#784

# adult brain ----
braindata = fread('data/braindata_adults.txt')
braindata %>% describe()#T1wSI (n=583)
source(file.path(wd,'scripts','[function]_define_fit_function_2025-01-09.R'))
d_T1wSI_2014 = fread("~/Downloads/jean_shin_dataanalysis_completespssubjects_20251203_115722.txt")
names(d_T1wSI_2014 )

WM_GM_volume = fread(file.path(wd,"data","lobarVolData_WM_GM_adults_2025-05-26.tsv"))
WM_GM_volume %>% describe
MTR_data = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/jean_shin_dataanalysis_completespssubjects_final_20251125_132248.txt")
head(MTR_data) #n=403
names(MTR_data) = str_remove(names(MTR_data),"archived.SPS-MRI-MTRpipe-MRI-Apr2014.")
MTR_data = MTR_data %>% mutate_at(-c(1:3),as.numeric)
MTR_data = MTR_data %>% mutate(across(c(wm, gm), ~ ifelse(.x < 0, NA, .x)))#n=398
MTR_data %>% describe

## adult volume data ----
braindata = fread(file.path(wd,"data","lobarVolData_WM_GM_adults_2025-05-26.tsv"))
describe(braindata )#403
# adult metabo and cov ----
metabodata = fread(file.path(wd,"data","metabodata_adults.txt"))
metabodata_info = fread(file.path(wd,"data","metabodata_info_adults.txt"))
metabodata %>% describe() %>% tail

## Nightingale
d_NG = fread("~/Downloads/jean_shin_dataanalysis_completespssubjects_final_20251203_112200.txt")
d_NG = d_NG %>% select(1:10)
names(d_NG) = str_remove(names(d_NG),"SPS-Biomarkers-2016_pipeline.")
d_NG = data.frame(d_NG)
d_NG[,-c(1:3)] = apply(d_NG[,-c(1:3)],2,as.numeric) 
d_NG %>% describe()

covdata = fread(file.path(wd,"data","covdata_adults.txt"))
covdata2 = fread(file.path(wd,"data","covdata_adults_v2.txt")) %>%
  mutate(fam_id = str_split(uniqueID,"_",simplify = T)[,1]) %>%
  mutate(sex = ifelse(Sex==0,"M","F"))
ynames=names(braindata)[-1]

#"~/Library/CloudStorage/OneDrive-SickKids/1Work/VF_PCA_CS/SPS QC'd Body and Brain Data.csv"

plot.df.volcano_main = fread("results/plot.df.volcano_main_2025-12-02.tsv")
plot.df.volcano_main %>% filter(generation=="ados") %>% pull(N) %>% range()
plot.df.volcano_main %>% filter(generation=="adults") %>% pull(N) %>% range()
                                                             