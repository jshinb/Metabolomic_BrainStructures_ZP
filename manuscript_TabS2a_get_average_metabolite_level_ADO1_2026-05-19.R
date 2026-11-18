#* Sex-specific lipidomic association analyses of white matter phenotypes in ADO1
#* Plan
#* 1. load brain data
#* 2. load metabolite data
#* 3. load covariate data
#* 4. run assoc analyses
#* 5. tabulate results
#* 6. get correlation coefficients among three types of esitmates

#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

#* 0. initialize ----
source(project_specific_file)
setwd(wd)

library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(patchwork)
library(tibble)#rownames_to_column
select=dplyr::select
describe=psych::describe

## 1. load brain data 
load("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/data/sys_d_merge_2026-04-01.Rdata")

d_merge = d_merge%>% 
  rename(Volume.WM = Volume.MW) 
d_merge = d_merge %>% select(-Volume.WM.adjAgeSexICV) %>% ## already INT)
  mutate( Volume.WM.adjAgeSexICV = # not inverse-normal transformed
            resid(lmer(lobar.vol.WM.adjICV~(age.c+age.c2)*Sex + (1|fam_id),
                       data = d_merge, na.action=na.exclude))) %>%
  mutate( T1wSI.WM.adjAgeSex = # not inverse-normal transformed
            resid(lmer(T1wSI.WM~(age.c+age.c2)*Sex + (1|fam_id),
                       data = d_merge, na.action=na.exclude))) %>%
  mutate(MTR.WM.adjAgeSex = # not inverse-normal transformed
           resid(lmer(LobarWM_Z~(age.c+age.c2)*Sex + (1|fam_id),
                      data = d_merge, na.action=na.exclude)))

d_merge %>% pull(Volume.WM.adjAgeSexICV) %>% hist
brain_ado = d_merge %>% 
  select(uniqueID,
         Volume.WM=lobar.vol.WM,
         Volume.WM.adjICV=lobar.vol.WM.adjICV,
         Volume.WM.adjAgeSexICV,## not INT
         T1wSI.WM,
         T1wSI.WM.adjAgeSex,
         MTR.WM = LobarWM_Z,
         MTR.WM.adjAgeSex,
         ICV=ICV_TNT2012,
         icv.c,
         icv.c2
  ) %>%
  mutate(icv.c=icv.c[,1]) %>%
  mutate(icv.c2=icv.c^2)

brain_ado %>% head

## 2. load metabodata
# to get the list of the tested metabolites (n=830)
res_WM1 = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/TabS1_2026-03-17.txt")
head(res_WM1 )
res_WM1 = res_WM1 %>% filter(metaboliteid!='GlycA')

metabodata_ado = fread('data/metabodata_ados.txt') 
names(metabodata_ado) = str_remove(names(metabodata_ado),"_nmol/mL")
metabodata_ado = metabodata_ado %>% select(all_of(c('uniqueID',res_WM1$metaboliteid)))
metabodata_ado %>% dim#830 metabolites

dim(brain_ado)
head(brain_ado)
ind = !is.na(brain_ado$Volume.WM.adjICV);print(sum(ind))
ind = ind | !is.na(brain_ado$T1wSI.WM);print(sum(ind))
ind = ind | !is.na(brain_ado$MTR.WM);print(sum(ind))#988
ID.wi.brain = brain_ado %>% filter(ind) %>% pull(uniqueID)#986

metabodata_wi_brain_ado = metabodata_ado %>% filter(uniqueID %in% ID.wi.brain)
metabodata_wi_brain_ado %>% dim#n=937

brain_wi_metabo = brain_ado %>% filter(uniqueID%in% metabodata_wi_brain_ado$uniqueID)
dim(brain_wi_metabo)#n=937

desc_metabodata_wi_brain_ado = metabodata_wi_brain_ado %>% select(-1) %>% describe %>%
  data.frame %>% tibble::rownames_to_column(var="Metabolite ID")
desc_metabodata_wi_brain_ado %>% head
desc_metabodata_wi_brain_ado %>% write_tsv('manuscript_tab_desc_metabodata_wi_brain.txt')

# adults ----
## brain data
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

ind = !is.na(brain_adult$Volume.WM.adjICV);sum( print(ind) )
ind = ind | !is.na(brain_adult$T1wSI.WM);sum( print(ind) )
ind = ind | !is.na(brain_adult$MTR.WM);sum( print(ind) )#603
ID.wi.brain_adult = brain_adult %>% filter(ind) %>% pull(uniqueID)

## metabodata_adults
metabodata_adult = fread('data/metabodata_adults.txt') 
names(metabodata_adult) = str_remove(names(metabodata_adult),"_nmol/mL")
metabodata_adult = metabodata_adult %>% select(all_of(c('uniqueID',res_WM1$metaboliteid)))
metabodata_adult[,1:5] %>% head; metabodata_adult %>% dim#830 metabolites
ID_incl = intersect(metabodata_adult$uniqueID,ID.wi.brain_adult);print(length(ID_incl))#586

metabodata_wi_brain_adult = metabodata_adult %>% filter(uniqueID %in% ID_incl)
metabodata_wi_brain_adult %>% dim# 586*831

desc_metabodata_wi_brain_adult = metabodata_wi_brain_adult %>% 
  select(-1) %>% 
  describe %>%
  data.frame %>% 
  tibble::rownames_to_column(var="Metabolite ID")
desc_metabodata_wi_brain_adult %>% head
desc_metabodata_wi_brain_adult %>% 
  write_tsv('manuscript_tab_desc_metabodata_wi_brain_adult.txt')

desc_metabodata_wi_brain_ado %>% select(`Metabolite ID`,mean,sd,median,mad,min,max) %>% 
  left_join(
    desc_metabodata_wi_brain_adult %>% select(`Metabolite ID`,mean,sd,median,mad,min,max), 
    join_by(`Metabolite ID`)
  ) %>% 
  write_tsv('manuscript_tab_desc_metabodata_wi_brain_ado_adult.txt')
