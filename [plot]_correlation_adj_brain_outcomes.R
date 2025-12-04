
rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
# outputs:
# "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/Correlation_brain_phenotypes_adjAgeSexICV_v3.pptx"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

# read in cov data ----
covdata = fread(file.path(wd,"data","covdata_ados.txt"))

#*[WM] ----
## read in result files: WM ----
load('results/ret_resid_unadjBMI_MTR_WM_ados_2025-04-08.Rdata')#mtr
WM.MTR_adjAge = ret_resid_unadjBMI$`LobarWM_Z MTR`$resid_out
rm(ret_resid_unadjBMI)
head(WM.MTR_adjAge)

load('results/ret_resid_unadjBMI_unadjICV_ados_2025-03-26.Rdata')#t1wSI
WM.T1wSI_adjAge = ret_resid_unadjBMI$NormWM$resid_out
head(WM.T1wSI_adjAge)
rm(ret_resid_unadjBMI)

load('results/ret_resid_unadjBMI_lobarVolWM_ados_2025-04-04.Rdata')#
WM.Vol_adjAge = ret_resid_unadjBMI$lobar.vol.WM.adjICV$resid_out
head(WM.Vol_adjAge)
rm(ret_resid_unadjBMI)

## merging:WM ----
braindata_adjAge = WM.T1wSI_adjAge %>% 
  left_join(WM.MTR_adjAge)%>% 
  left_join(WM.Vol_adjAge)

dim(braindata_adjAge)
head(braindata_adjAge)

## braindata_adjAge_wi_AgeSex(adj for sex): WM ----
braindata_adjAge_wi_AgeSex = covdata%>% dplyr::select(uniqueID,age,Sex) %>%
  left_join(braindata_adjAge) %>% mutate(age.c = scale(age)/2) %>%
  mutate(age.c2 = age.c^2)
names.braindata_adjAge_wi_AgeSex = names(braindata_adjAge_wi_AgeSex)
braindata_adjAge_wi_AgeSex = data.frame(braindata_adjAge_wi_AgeSex)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
names(braindata_adjAge_wi_AgeSex) = names.braindata_adjAge_wi_AgeSex
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(NormWM_fam~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),NormWM_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(inormal(`LobarWM_Z MTR_fam`)~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
nobs(fit1)#784
op=par(mfrow=c(2,2));plot(fit1,pch=20);par(op)
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),MTR_WM_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(lobar.vol.WM.adjICV_fam~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
nobs(fit1)#870
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),WMvol_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

p1=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=NormWM_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p2=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=MTR_WM_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p3=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=WMvol_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p1 + p2 + p3 + plot_layout(guides = "collect")

headTail(braindata_adjAge_wi_AgeSex)

## sex-combined correlation plot: WM ----
braindata_adjAge_wi_AgeSex_cleaned = braindata_adjAge_wi_AgeSex %>% 
  select(uniqueID,NormWM_adjAgeSex,WMvol_adjAgeSex,MTR_WM_adjAgeSex) %>%
  rename(nT1wSI_adjAgeSex = NormWM_adjAgeSex,
         VolWM_adjAgeSexICV = WMvol_adjAgeSex,
         MTR_adjAgeSex = MTR_WM_adjAgeSex)
braindata_adjAge_wi_AgeSex_cleaned = remove.outliers_grubbs(
  braindata_adjAge_wi_AgeSex_cleaned,
  varnames = names(braindata_adjAge_wi_AgeSex_cleaned)[-1])
braindata_adjAge_wi_AgeSex_cleaned$counts.NA#0,1,0
head(braindata_adjAge_wi_AgeSex_cleaned$clean_data)#64,

ggcorr_all = create_ggcorrplot (
  data = braindata_adjAge_wi_AgeSex_cleaned$clean_data ,
  sel_colnames = c("VolWM_adjAgeSexICV", "nT1wSI_adjAgeSex", "MTR_adjAgeSex"),
  new_colnames = c("WM-Vol","WM-SI","WM-MTR"),
  ID_colname = "uniqueID"
)

## female correlation plot ----
subgroup = 'F';sub_ids = covdata %>% filter(Sex==subgroup)%>% pull(uniqueID)
sel_colnames = c('lobar.vol.WM.adjICV_fam','NormWM_fam','LobarWM_Z MTR_fam')
braindata_adjAge_cleaned = subset(
  braindata_adjAge %>% filter(uniqueID %in% sub_ids),
  select = c('uniqueID',sel_colnames))
braindata_adjAge_cleaned = remove.outliers_grubbs(
  braindata_adjAge_cleaned, varnames = sel_colnames)
braindata_adjAge_cleaned $counts.NA %>% mutate(diff= after-before)
#                       var before after diff
# 1 lobar.vol.WM.adjICV_fam     81    83    2
# 2              NormWM_fam     33    33    0
# 3       LobarWM_Z MTR_fam    124   125    1

ggcorr_female = create_ggcorrplot (
  data = braindata_adjAge_cleaned$clean_data ,
  sel_colnames = c('lobar.vol.WM.adjICV_fam','NormWM_fam','LobarWM_Z MTR_fam'),
  new_colnames = c("WM-Vol","WM-SI","WM-MTR"),
  ID_colname = 'uniqueID'
)

## male correlation plot ----
subgroup = 'M';sub_ids = covdata %>% filter(Sex==subgroup)%>% pull(uniqueID)#496
sel_colnames = c('lobar.vol.WM.adjICV_fam','NormWM_fam','LobarWM_Z MTR_fam')
braindata_adjAge_cleaned = subset(
  braindata_adjAge %>% filter(uniqueID %in% sub_ids),
  select = c('uniqueID',sel_colnames))

braindata_adjAge_cleaned = remove.outliers_grubbs(
  braindata_adjAge_cleaned, varnames = sel_colnames)
braindata_adjAge_cleaned $counts.NA %>% mutate(diff= after-before)
#                       var before after diff
# 1 lobar.vol.WM.adjICV_fam     78    78    0
# 2              NormWM_fam     31    31    0
# 3       LobarWM_Z MTR_fam    121   124    3
ggcorr_male = create_ggcorrplot (
  data = braindata_adjAge_cleaned$clean_data ,
  sel_colnames = sel_colnames,
  new_colnames = c("WM-Vol","WM-SI","WM-MTR"),
  ID_colname = 'uniqueID'
)

## create png file for all three graphs ----
ggcorr_plot.WM =
  (ggcorr_all$p_ggcorr + theme(axis.text.y = element_text(size = 10)))+ 
  (ggcorr_female$p_ggcorr + theme(axis.text.y = element_blank())) + 
  (ggcorr_male$p_ggcorr + theme(axis.text.y = element_blank())) + 
  plot_layout(guides = 'collect') &
  theme(axis.text.x = element_text(size = 10))

ggsave(file.path(opt$output.dir,paste0("brain_corrplot_WM.jpg")), 
       bg = "white",
       plot = ggcorr_plot.WM, 
       width = 8.5*3*1.25, 
       height = 4.5*1*1.25, 
       units = "cm", 
       dpi = 300)

#*[GM] ----
## read in result files: GM ----
load('results/ret_resid_unadjBMI_LobarGM_Z_MTR_ados_2025-11-25.Rdata')#mtr
GM.MTR_adjAge = ret_resid_unadjBMI$LobarGM_Z$resid_out
rm(ret_resid_unadjBMI)
head(GM.MTR_adjAge)

load('results/ret_resid_unadjBMI_unadjICV_ados_2025-03-26.Rdata')#t1wSI
GM.T1wSI_adjAge = ret_resid_unadjBMI$NormGM$resid_out
head(GM.T1wSI_adjAge)
rm(ret_resid_unadjBMI)

load('results/ret_resid_unadj_lobarVolGM_ados_2025-05-26.Rdata')#
ret_resid_unadjBMI=ret_resid_unadj;rm(ret_resid_unadj)
names(ret_resid_unadjBMI)
GM.Vol_adjAge = ret_resid_unadjBMI$lobar.vol.GM.adjICV$resid_out
head(GM.Vol_adjAge)
rm(ret_resid_unadjBMI)

## braindata_adjAge_wi_AgeSex (adj for sex): GM ----
GM.braindata_adjAge = GM.Vol_adjAge %>% 
  left_join(GM.T1wSI_adjAge)%>% 
  left_join(GM.MTR_adjAge)

braindata_adjAge_wi_AgeSex = covdata %>% dplyr::select(uniqueID,age,Sex) %>%
  left_join(GM.braindata_adjAge) %>% mutate(age.c = scale(age)/2) %>%
  mutate(age.c2 = age.c^2)

names.braindata_adjAge_wi_AgeSex = names(braindata_adjAge_wi_AgeSex)
braindata_adjAge_wi_AgeSex = data.frame(braindata_adjAge_wi_AgeSex)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
names(braindata_adjAge_wi_AgeSex) = names.braindata_adjAge_wi_AgeSex
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(NormGM_fam~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),NormGM_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(inormal(`LobarGM_Z_fam`)~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
nobs(fit1)#769
op=par(mfrow=c(2,2));plot(fit1,pch=20);par(op)
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),MTR_GM_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(lobar.vol.GM.adjICV_fam~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
nobs(fit1)#870
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),GMvol_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

p1=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=NormGM_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p2=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=MTR_GM_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p3=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=GMvol_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p1 + p2 + p3 + plot_layout(guides = "collect")

headTail(braindata_adjAge_wi_AgeSex)

## sex-combined correlation plot:GM ----
braindata_adjAge_wi_AgeSex_cleaned = braindata_adjAge_wi_AgeSex %>% 
  select(uniqueID,NormGM_adjAgeSex,GMvol_adjAgeSex,MTR_GM_adjAgeSex) %>%
  rename(nT1wSI_adjAgeSex = NormGM_adjAgeSex,
         VolGM_adjAgeSexICV = GMvol_adjAgeSex,
         MTR_adjAgeSex = MTR_GM_adjAgeSex)
braindata_adjAge_wi_AgeSex_cleaned = remove.outliers_grubbs(
  braindata_adjAge_wi_AgeSex_cleaned,
  varnames = names(braindata_adjAge_wi_AgeSex_cleaned)[-1])
braindata_adjAge_wi_AgeSex_cleaned$counts.NA#1,2,0
head(braindata_adjAge_wi_AgeSex_cleaned$clean_data)#64,

ggcorr_all.GM = create_ggcorrplot (
  data = braindata_adjAge_wi_AgeSex_cleaned$clean_data ,
  sel_colnames = c("VolGM_adjAgeSexICV", "nT1wSI_adjAgeSex", "MTR_adjAgeSex"),
  new_colnames = c("GM-Vol","GM-T1wSI","GM-MTR"),
  ID_colname = "uniqueID"
)
ggsave(file.path(opt$output.dir,paste0("brain_corrplot_GM.jpg")), 
       bg = "white",
       plot = (ggcorr_all.GM$p_ggcorr+theme(axis.text.y = element_text(size = 10))), 
       width = 8.5*1*1.25, 
       height = 4.5*1*1.25, 
       units = "cm", 
       dpi = 300)
## scatter plot ----
head(braindata_adjAge_wi_AgeSex_cleaned$clean_data)
df = braindata_adjAge_wi_AgeSex_cleaned$clean_data %>%
  rename("GM-Vol"=VolGM_adjAgeSexICV,
         'GM-T1wSI'=nT1wSI_adjAgeSex,
         'GM-MTR'=MTR_adjAgeSex)
scatter1 = df %>% 
  ggplot(aes(x=`GM-Vol`,y=`GM-T1wSI`)) + 
  geom_point() + geom_smooth(method = 'gam') + stat_cor()
scatter2 = df %>% 
  ggplot(aes(x=`GM-Vol`,y=`GM-MTR`)) + 
  geom_point() + geom_smooth(method = 'gam') + stat_cor()
scatter3 = df %>% 
  ggplot(aes(x=`GM-T1wSI`,y=`GM-MTR`)) + 
  geom_point() + geom_smooth(method = 'gam') + stat_cor()

p = scatter1 + scatter2 + scatter3
ggsave(file.path(opt$output.dir,paste0("brain_scatterplot_brain_GM.jpg")), 
       bg = "white",
       plot = p,
       width = 8.5*3*1.25, 
       height = 4.5*1.5*1.25, 
       units = "cm", 
       dpi = 300)

#*[GM adults] ----
## read in result files: GM ----
load('results/ret_resid_unadjBMI_GM_WM_MTR_adults_wi_lipidLowMed_2025-11-25.Rdata')#mtr
GM.MTR_adjAge = ret_resid_unadjBMI$gm$resid_out
rm(ret_resid_unadjBMI)
head(GM.MTR_adjAge)

load('results/ret_resid_T1wSI_unadjBMI_adults_wi_lipidLowMed_2025-01-09.Rdata')#t1wSI
GM.T1wSI_adjAge = ret_resid_unadjBMI$NormGM$resid_out
head(GM.T1wSI_adjAge)
rm(ret_resid_unadjBMI)

load('results/')#volume
head(ret_resid_unadjBMI)
GM.Vol_adjAge = ret_resid_unadjBMI$lobar.vol.GM.adjICV$resid_out
head(GM.Vol_adjAge)
rm(ret_resid_unadjBMI)

## braindata_adjAge_wi_AgeSex (adj for sex): GM ----
GM.braindata_adjAge = GM.Vol_adjAge %>% 
  left_join(GM.T1wSI_adjAge)%>% 
  left_join(GM.MTR_adjAge)

braindata_adjAge_wi_AgeSex = covdata %>% dplyr::select(uniqueID,age,Sex) %>%
  left_join(GM.braindata_adjAge) %>% mutate(age.c = scale(age)/2) %>%
  mutate(age.c2 = age.c^2)

names.braindata_adjAge_wi_AgeSex = names(braindata_adjAge_wi_AgeSex)
braindata_adjAge_wi_AgeSex = data.frame(braindata_adjAge_wi_AgeSex)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
names(braindata_adjAge_wi_AgeSex) = names.braindata_adjAge_wi_AgeSex
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(NormGM_fam~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),NormGM_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(inormal(`LobarGM_Z_fam`)~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
nobs(fit1)#769
op=par(mfrow=c(2,2));plot(fit1,pch=20);par(op)
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),MTR_GM_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

fit1 = lm(lobar.vol.GM.adjICV_fam~(age.c + age.c2)* Sex, na.action = na.exclude, data=braindata_adjAge_wi_AgeSex)
nobs(fit1)#870
resid_fit1 = data.frame(uniqueID = names(resid(fit1)),GMvol_adjAgeSex=resid(fit1))
braindata_adjAge_wi_AgeSex = braindata_adjAge_wi_AgeSex %>% left_join(resid_fit1)
rownames(braindata_adjAge_wi_AgeSex) = braindata_adjAge_wi_AgeSex$uniqueID
head(braindata_adjAge_wi_AgeSex)

p1=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=NormGM_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p2=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=MTR_GM_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p3=braindata_adjAge_wi_AgeSex %>% ggplot(aes(x=age,y=GMvol_adjAgeSex,color=Sex,fill=Sex)) + geom_smooth(method="gam") + geom_point(alpha=0.3)
p1 + p2 + p3 + plot_layout(guides = "collect")

headTail(braindata_adjAge_wi_AgeSex)

## sex-combined correlation plot:GM ----
braindata_adjAge_wi_AgeSex_cleaned = braindata_adjAge_wi_AgeSex %>% 
  select(uniqueID,NormGM_adjAgeSex,GMvol_adjAgeSex,MTR_GM_adjAgeSex) %>%
  rename(nT1wSI_adjAgeSex = NormGM_adjAgeSex,
         VolGM_adjAgeSexICV = GMvol_adjAgeSex,
         MTR_adjAgeSex = MTR_GM_adjAgeSex)
braindata_adjAge_wi_AgeSex_cleaned = remove.outliers_grubbs(
  braindata_adjAge_wi_AgeSex_cleaned,
  varnames = names(braindata_adjAge_wi_AgeSex_cleaned)[-1])
braindata_adjAge_wi_AgeSex_cleaned$counts.NA#1,2,0
head(braindata_adjAge_wi_AgeSex_cleaned$clean_data)#64,

ggcorr_all.GM = create_ggcorrplot (
  data = braindata_adjAge_wi_AgeSex_cleaned$clean_data ,
  sel_colnames = c("VolGM_adjAgeSexICV", "nT1wSI_adjAgeSex", "MTR_adjAgeSex"),
  new_colnames = c("GM-Vol","GM-T1wSI","GM-MTR"),
  ID_colname = "uniqueID"
)
ggsave(file.path(opt$output.dir,paste0("brain_corrplot_GM.jpg")), 
       bg = "white",
       plot = (ggcorr_all.GM$p_ggcorr+theme(axis.text.y = element_text(size = 10))), 
       width = 8.5*1*1.25, 
       height = 4.5*1*1.25, 
       units = "cm", 
       dpi = 300)
## scatter plot ----
head(braindata_adjAge_wi_AgeSex_cleaned$clean_data)
df = braindata_adjAge_wi_AgeSex_cleaned$clean_data %>%
  rename("GM-Vol"=VolGM_adjAgeSexICV,
         'GM-T1wSI'=nT1wSI_adjAgeSex,
         'GM-MTR'=MTR_adjAgeSex)
scatter1 = df %>% 
  ggplot(aes(x=`GM-Vol`,y=`GM-T1wSI`)) + 
  geom_point() + geom_smooth(method = 'gam') + stat_cor()
scatter2 = df %>% 
  ggplot(aes(x=`GM-Vol`,y=`GM-MTR`)) + 
  geom_point() + geom_smooth(method = 'gam') + stat_cor()
scatter3 = df %>% 
  ggplot(aes(x=`GM-T1wSI`,y=`GM-MTR`)) + 
  geom_point() + geom_smooth(method = 'gam') + stat_cor()

p = scatter1 + scatter2 + scatter3
ggsave(file.path(opt$output.dir,paste0("brain_scatterplot_brain_GM.jpg")), 
       bg = "white",
       plot = p,
       width = 8.5*3*1.25, 
       height = 4.5*1.5*1.25, 
       units = "cm", 
       dpi = 300)

