
rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

# read in cov data ----
covdata = fread(file.path(wd,"data","covdata_ados.txt"))

# read in result files ----
load('results/ret_resid_unadjBMI_MTR_WM_ados_2025-04-08.Rdata')#mtr
mtrWM_adjAge = ret_resid_unadjBMI$`LobarWM_Z MTR`$resid_out
rm(ret_resid_unadjBMI)
head(mtrWM_adjAge)

load('results/ret_resid_unadjBMI_unadjICV_ados_2025-03-26.Rdata')#t1wSI
NormWM_adjAge = ret_resid_unadjBMI$NormWM$resid_out
head(NormWM_adjAge)
rm(ret_resid_unadjBMI)

load('results/ret_resid_unadjBMI_lobarVolWM_ados_2025-04-04.Rdata')#
WMvol_adjAge = ret_resid_unadjBMI$lobar.vol.WM.adjICV$resid_out
head(WMvol_adjAge)
rm(ret_resid_unadjBMI)

# merging ----
braindata_adjAge = NormWM_adjAge %>% 
  left_join(mtrWM_adjAge)%>% 
  left_join(WMvol_adjAge)

dim(braindata_adjAge)
head(braindata_adjAge)

# braindata_adjAge_wi_AgeSex: adj for sex ----
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

# braindata_adjAge_wi_AgeSex_cleaned (outliers have been removed) ----
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


# sex-combined correlation plot ----
sel_colnames = c("VolWM_adjAgeSexICV", "nT1wSI_adjAgeSex", "MTR_adjAgeSex")
df = braindata_adjAge_wi_AgeSex_cleaned$clean_data %>% select(all_of(sel_colnames));rm(sel_colnames)
ind = !is.na(df[[1]])
ind = ind | !is.na(df[[2]])
ind = ind | !is.na(df[[3]])
df = df[ind,]

## change the order of columns 
names(df) = c("WM-Vol", "WM-T1wSI", "WM-MTR")
cat(names(df),sep="\n")
subgroup = 'sex-combined'
vnames =  names(df)

suppressWarnings({
  # create correlation plot 
  heat = plot_heatmap(dat = df, axisNames = vnames)
  
  # annotate
  plot.title=paste0(opt$group_name, " (", subgroup, ", N = ",nrow(df),")")
  heat <- annotate_figure(
    heat, 
    top = text_grob(plot.title,
                    color = "black", 
                    face = "bold", 
                    size = 10)
  ) 
  
  # save plot
  ggsave(file.path(opt$output.dir,"brain_corrplot.jpg"), 
         bg = "white",
         plot = heat, 
         width = 8.5*1.35, 
         height = 4.5**1.35, 
         units = "cm", 
         dpi = 300)
})

# female correlation plot ----
subgroup = 'F'
sub_ids = covdata %>% filter(Sex==subgroup)%>% pull(uniqueID)
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

rownames.df = braindata_adjAge_cleaned$clean_data$uniqueID
df = braindata_adjAge_cleaned$clean_data %>% select(all_of(sel_colnames))
df = data.frame(df)
rownames(df) = rownames.df 

ind = !is.na(df[[1]])
ind = ind | !is.na(df[[2]])
ind = ind | !is.na(df[[3]])
df = df[ind,]

## change the order of columns 
names(df) = c("WM-Vol", "WM-T1wSI", "WM-MTR")
cat(names(df),sep="\n")
vnames =  names(df)

suppressWarnings({
  # create correlation plot 
  heat = plot_heatmap(dat = df, axisNames = vnames)
  
  # annotate
  plot.title=paste0(opt$group_name, " (", subgroup, ", N = ",nrow(df),")")
  heat <- annotate_figure(
    heat, 
    top = text_grob(plot.title,
                    color = "black", 
                    face = "bold", 
                    size = 10)
  ) 
  
  # save plot
  ggsave(file.path(opt$output.dir,paste0("brain_corrplot_",subgroup,".jpg")), 
         bg = "white",
         plot = heat, 
         width = 8.5*1.35, 
         height = 4.5**1.35, 
         units = "cm", 
         dpi = 300)
})

# male correlation plot ----
subgroup = 'M'
sub_ids = covdata %>% filter(Sex==subgroup)%>% pull(uniqueID)#496
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
rownames.df = braindata_adjAge_cleaned$clean_data$uniqueID
df = braindata_adjAge_cleaned$clean_data %>% select(all_of(sel_colnames))
df = data.frame(df)
rownames(df) = rownames.df 

ind = !is.na(df[[1]])
ind = ind | !is.na(df[[2]])
ind = ind | !is.na(df[[3]])
df = df[ind,]#476

## change the order of columns 
names(df) = c("WM-Vol", "WM-T1wSI", "WM-MTR")
cat(names(df),sep="\n")
vnames =  names(df)

suppressWarnings({
  # create correlation plot 
  heat = plot_heatmap(dat = df, axisNames = vnames)
  
  # annotate
  plot.title=paste0(opt$group_name, " (", subgroup, ", N = ",nrow(df),")")
  heat <- annotate_figure(
    heat, 
    top = text_grob(plot.title,
                    color = "black", 
                    face = "bold", 
                    size = 10)
  ) 
  
  # save plot
  ggsave(file.path(opt$output.dir,paste0("brain_corrplot_",subgroup,".jpg")), 
         bg = "white",
         plot = heat, 
         width = 8.5*1.35, 
         height = 4.5**1.35, 
         units = "cm", 
         dpi = 300)
})
