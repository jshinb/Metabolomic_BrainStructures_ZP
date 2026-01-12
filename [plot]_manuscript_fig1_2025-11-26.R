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

add_class2 = function(volcano_main_all){
  metabodata_info = fread("data/metabodata_info_ados.txt")
  
  if(sum(names(volcano_main_all)=='class')==1){
    volcano_main_all = volcano_main_all %>% 
      dplyr::select(-class) 
  }
  volcano_main_all = volcano_main_all %>% 
    mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>%
    left_join(
      metabodata_info %>%
        mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>% 
        dplyr::select(metabolite_id,class)
    )
  # class2.levels = c("CE", 
  #                   "TAG, DAG, or MAG", 
  #                   "PC, PE, or PI",
  #                   "LPC or LPE", 
  #                   "SM or ceramides",
  #                   "Nightingale")
  class2.levels = c("Cholesteryl esters",
                    "Acylglycerols",
                    "Phospholipids",
                    "Lysophospholipids",
                    "Sphingolipids",
                    "Inflammation/AA/metabolism")
  volcano_main_all[['class2']] <- class2.levels[1]
  volcano_main_all = volcano_main_all %>%
    mutate(class2 = ifelse(class %in% c('TAG','DAG','MAG'),class2.levels[2],class2)) %>%
    mutate(class2 = ifelse(class %in% c('PC','PE','PI'),class2.levels[3],class2)) %>%
    mutate(class2 = ifelse(class %in% c('LPC','LPE'),class2.levels[4],class2)) %>%
    mutate(class2 = ifelse(class %in% c('SM','CER',"DCER","HCER","LCER"),class2.levels[5],class2)) %>%
    mutate(class2 = ifelse(class %in% c('Amino acids','Glycolysis related metabolites','Inflammation','Ketone bodies'),
                           class2.levels[6],class2))
  table(volcano_main_all$class,volcano_main_all$class2,useNA = 'a')
  volcano_main_all = volcano_main_all %>%
    mutate(`lipid class` = factor(class2,levels=class2.levels)) 
  cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
  
  class2.labels = paste(names(table(volcano_main_all$`lipid class`))," (n=",
                        table(volcano_main_all$`lipid class`),")",sep='')
  ret=(list(volcano_main_all=volcano_main_all,class2.labels=class2.labels,cols=cols))
  ret
}
# load significant lipids ----
load('results/quad.lipids_NormWM_MTR_WMvol_assoc.Rdata')
siglipids = c(quad.lipids$quad1$metabolite_id,quad.lipids$quad2$metabolite_id,quad.lipids$quad4$metabolite_id)#616

# load volcano results ----
load("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
volcano_main_NormWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.Rdata")
volcano_main_lobarVolWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_LobarWM_Z MTR_sex-combined_ados_unadjBMI_2025-04-08.Rdata")
volcano_main_MTR_WM = volcano_main
rm(volcano_main)

#debug(add_class2)
plot.df.volcano_main_NormWM = add_class2(volcano_main_NormWM$data)
head(plot.df.volcano_main_NormWM$volcano_main_all)

plot.df.volcano_main_lobarVolWM = add_class2(volcano_main_lobarVolWM$data)
head(plot.df.volcano_main_lobarVolWM$volcano_main_all)

plot.df.volcano_main_MTR_WM = add_class2(volcano_main_MTR_WM$data)
volcano_main_NormWM$data = plot.df.volcano_main_NormWM$volcano_main_all

## top metabolites 
top_main_lobarVolWM = plot.df.volcano_main_lobarVolWM$volcano_main_all %>% 
  filter(beta>0) %>% arrange(pvalue) %>% slice(1:10)
top_main_lobarVolWM = top_main_lobarVolWM %>% 
  bind_rows(
    plot.df.volcano_main_lobarVolWM$volcano_main_all %>% 
      filter(beta<0) %>% arrange(pvalue) %>% slice(1:10)
  ) 

top_main_NormWM = plot.df.volcano_main_NormWM$volcano_main_all %>% 
  filter(beta>0) %>% arrange(pvalue) %>% slice(1:10)
top_main_NormWM = top_main_NormWM %>% bind_rows(
  plot.df.volcano_main_NormWM$volcano_main_all %>% 
    filter(beta<0) %>% arrange(pvalue) %>% slice(1:10)
) 

top_main_MTR_WM = plot.df.volcano_main_MTR_WM$volcano_main_all %>% 
  filter(beta>0) %>% arrange(pvalue) %>% slice(1:10)
top_main_MTR_WM = top_main_MTR_WM %>% 
  bind_rows(
    plot.df.volcano_main_MTR_WM$volcano_main_all %>% 
      filter(beta<0) %>% arrange(pvalue) %>% slice(1:10)
  ) 

metabolite_overlap = intersect(
  top_main_lobarVolWM$metabolite_id,
  top_main_NormWM$metabolite_id
)
metabolite_overlap = intersect(
  metabolite_overlap,
  top_main_MTR_WM$metabolite_id
)

top_main_lobarVolWM = top_main_lobarVolWM %>% 
  mutate(label.col = ifelse(metabolite_id %in% metabolite_overlap,"darkred","grey30"))
top_main_NormWM = top_main_NormWM %>% 
  mutate(label.col = ifelse(metabolite_id %in% metabolite_overlap,"darkred","grey30"))
top_main_MTR_WM = top_main_MTR_WM %>% 
  mutate(label.col = ifelse(metabolite_id %in% metabolite_overlap,"darkred","grey30"))


# load association results ----
## T1wSI volcano plot data frame -----
### ado data (found typo on Nov 14, 2025) ----
load("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
volcano_main_NormWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_NormGM_sex-combined_ados_adjBMI_unadjICV_2025-03-26.Rdata")
volcano_main_NormGM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.Rdata")
volcano_main_VolumeWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_ados_unadjBMI_2025-05-26.Rdata")
volcano_main_VolumeGM = volcano_main
rm(volcano_main)

### adult data ----
load("results/volcanoplot_NormWM_sex-combined_adults_unadjBMI_2025-01-10.Rdata")
volcano_main_NormWM_adults = volcano_main
rm(volcano_main)
table(volcano_main_NormWM_adults$data$N)

load("results/volcanoplot_NormGM_sex-combined_adults_unadjBMI_2025-01-10.Rdata")
volcano_main_NormGM_adults = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
volcano_main_VolumeWM_adults = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
volcano_main_VolumeGM_adults = volcano_main
rm(volcano_main)

dim(fread('data/metabodata_adults.txt'))#650
dim(fread('data/metabodata_ados.txt'))

table(volcano_main_VolumeWM_adults$data$N)
table(volcano_main_VolumeGM_adults$data$N)#267-361
table(volcano_main_VolumeWM$data$N)#669-826
table(volcano_main_VolumeGM$data$N)#669-826

## Volume: volcano plot data frame ----
## ados
load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.Rdata")
volcano_main_VolumeWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_ados_unadjBMI_2025-05-26.Rdata")
volcano_main_VolumeGM = volcano_main
rm(volcano_main)

## adults: 
load("results/volcanoplot_all_lobar.vol.WM_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
volcano_main_VolumeWM_adults = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
volcano_main_VolumeGM_adults = volcano_main
rm(volcano_main)

# brain-by-age ----
load('data/braindata_NormWM_2025-03-26.Rdata'); 
braindata_tmp  = braindata
load('data/braindata_MTR_WM_2025-04-08.Rdata'); 
braindata_tmp  = braindata_tmp %>% left_join(braindata, join_by(uniqueID))
describe(braindata_tmp)
load('data/braindata_WMvol_2025-04-04.Rdata'); 
braindata_tmp  = braindata_tmp %>% left_join(braindata, join_by(uniqueID))
describe(braindata_tmp)
covdata = fread('data/covdata_ados.txt')

# scatter plots: brain-by-brain ----
##lobar.vol.WM_vs_T1wSI
scatter_Volume_vs_T1wSI = braindata_tmp %>% ggplot(aes(x=lobar.vol.WM,y=NormWM)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method='gam') +
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  # xlab("WM volume") + ylab("WM nT1wSI") #+ ggpubr::stat_cor()
  xlab(NULL) + ylab(NULL) #+ ggpubr::stat_cor()

##lobar.vol.WM_vs_T1wSI
scatter_Volume_vs_MTR = braindata_tmp %>% 
  ggplot(aes(x=lobar.vol.WM,y=`LobarWM_Z MTR`)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method='gam') +
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  # xlab("Volume") + ylab("MTR") #+ ggpubr::stat_cor()
  xlab(NULL) + ylab(NULL) #+ ggpubr::stat_cor()
scatter_Volume_vs_MTR

##
scatter_T1wSI_vs_MTR = braindata_tmp %>% 
  ggplot(aes(x=NormWM,y=`LobarWM_Z MTR`)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method='gam') +
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  # xlab("T1wSI") + ylab("MTR") #+ ggpubr::stat_cor()
  xlab(NULL) + ylab(NULL) #+ ggpubr::stat_cor()

scatter_Volume_vs_T1wSI + scatter_Volume_vs_MTR + scatter_T1wSI_vs_MTR  +
  plot_layout(guides='collect') &
  theme(legend.position = 'none') &
  # coord_cartesian(ylim=c(-0.225,0.35),xlim=c(-0.225,0.225)) &
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=9),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) #&

# geom_abline(slope=1,intercept = 0, color='grey')
ggsave("results/raw_brainphenos_scatterplot_sex_combined.png",
       width=11.5,height=3.5,units='in',dpi=300)

# scatter plots: brain-by-age (raw) ----
## age vs. volume 
scatter_age_vs_Volume = braindata_tmp %>% 
  left_join(covdata) %>%
  ggplot(aes(x=age,y=lobar.vol.WM)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method='gam') +
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  xlab('age (years)') + ylab('WM Volume') #+ ggpubr::stat_cor()

#age_vs_T1wSI
scatter_age_vs_T1wSI = braindata_tmp %>% 
  left_join(covdata) %>%
  ggplot(aes(x=age,y=NormWM)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method='gam') +
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  # xlab("WM volume") + ylab("WM nT1wSI") #+ ggpubr::stat_cor()
  xlab('age (years)') + ylab('WM nT1wSI') #+ ggpubr::stat_cor()

##age_vs_MTR
scatter_age_vs_MTR = braindata_tmp %>% 
  left_join(covdata) %>%
  ggplot(aes(x=age,y=`LobarWM_Z MTR`)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method='gam') +
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  # xlab("Volume") + ylab("MTR") #+ ggpubr::stat_cor()
  xlab('age (years)') + ylab('WM MTR') #+ ggpubr::stat_cor()

##
scatter_age_vs_Volume + scatter_age_vs_T1wSI + scatter_age_vs_MTR  +
  plot_layout(guides='collect',axis_titles = "collect") &
  theme(legend.position = 'none') &
  # coord_cartesian(ylim=c(-0.225,0.35),xlim=c(-0.225,0.225)) &
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=9),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) #& ggpubr::stat_cor()#&

# geom_abline(slope=1,intercept = 0, color='grey')
ggsave("results/raw_age_vs_brainphenos_scatterplot_sex_combined.png",
       width=11.5,height=3.5,units='in',dpi=300)

## plot data frame ----
plot.df = braindata_tmp %>% left_join(covdata,join_by(uniqueID))
describe(plot.df)

plot.df.cleaned = remove.outliers_grubbs(
  plot.df,
  varnames=c('NormWM','NormGM','LobarWM_Z MTR','lobar.vol.WM.adjICV')
  )
plot.df.cleaned$counts.NA#0,0,1,7
point.size=1.25
## NormWM ----
sm=describe(plot.df.cleaned$cleaed_data$NormWM,na.rm=T)
ysep = (sm$max-sm$min)*0.02;sm$max = sm$max+4*ysep
y1 = sm$max+ysep;y0 = sm$max-3*ysep

NormWM_by_age = plot.df %>%
  ggplot(aes(x=age,y=NormWM,color=Sex)) + 
  coord_cartesian(ylim=c(NA,y1)) + 
  # ggpubr::stat_cor(label.x.npc = 0,label.y = c(y0,y1)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(values=c('darkred','darkblue')) + 
  geom_smooth(method='gam') + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  ylab("WM nT1wSI")

## MTR----
sm=describe(plot.df.cleaned$cleaed_data$`LobarWM_Z MTR`,na.rm=T)
ysep = (sm$max-sm$min)*0.02;sm$max = sm$max+4*ysep
y1 = sm$max+ysep;y0 = sm$max-3*ysep

MTR_WM_by_age =  plot.df.cleaned$cleaed_data %>%
  ggplot(aes(x=age,y=`LobarWM_Z MTR`,color=Sex)) + 
  coord_cartesian(ylim=c(NA,y1)) + 
  # ggpubr::stat_cor(label.x.npc = 0,label.y = c(y0,y1)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(values=c('darkred','darkblue')) + 
  geom_smooth(method='gam') + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1))  +
  ylab("WM MTR") 

## WMvol ----
yname='lobar.vol.WM.adjICV'
sm=describe((plot.df.cleaned$cleaed_data[[yname]]),na.rm=T)
ysep = (sm$max-sm$min)*0.02;sm$max = sm$max+4*ysep
y1 = sm$max+ysep;y0 = sm$max-3*ysep

volWM_by_age =  plot.df.cleaned$cleaed_data %>%
  ggplot(aes(x=age,y=(lobar.vol.WM.adjICV),color=Sex)) + 
  coord_cartesian(ylim=c(NA,y1)) + 
  # ggpubr::stat_cor(label.x.npc = 0,label.y = c(y0,y1)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(values=c('darkred','darkblue')) + 
  geom_smooth(method='gam') + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) +
  ylab('Volume')

## generate png file ----
NormWM_by_age+MTR_WM_by_age+volWM_by_age+plot_layout(axes='collect',guides='collect')
ggsave("results/brain_by_age_sex_combined_2025-05-21.png",
       width=11.5,height=3.5*1.15,units='in',dpi=300)

hist(plot.df.cleaned$cleaed_data$`LobarWM_Z MTR`)
hist(plot.df.cleaned$cleaed_data$NormWM)

# metabolite ----
testing.lipids = fread("testing_lipid_species_list.txt",header=F) %>%
  pull(V1)
length(testing.lipids)#984 -> 838

metabodata_info = fread(file.path(wd,"data","metabodata_info_ados.txt"))
head(metabodata_info)

metabodata = fread(file.path(wd,"data","metabodata_ados.txt"))
metabodata = subset(
  metabodata,
  select=c('uniqueID',metabodata_info$metabolite_id)
)#

names(metabodata) = str_remove(names(metabodata),"_nmol/mL")
names.metabodata = names(metabodata)
metabodata = data.frame(metabodata)
names(metabodata) = names.metabodata
rownames(metabodata) = metabodata$uniqueID
testing.lipids2 = volcano_main_NormWM$data$metabolite_id

pca_metabo_adults = FactoMineR::PCA(metabodata %>% select(all_of(testing.lipids2)))
d = metabodata  %>% 
  left_join(covdata %>% dplyr::select(uniqueID,age,age.c,age.c2,Sex), by="uniqueID")
rownames(d) = d$uniqueID
head(d[,1:3])
metabo.names = names.metabodata[-1]
#n=90 
1/90

headTail(setdiff(testing.lipids[1:19],testing.lipids2[1:19]))
sort(testing.lipids2[1:19])

metabodata_info_test = metabodata_info%>% 
  filter(str_remove(metabolite_id,"_nmol/mL") %in% testing.lipids2) 

metabodata_info_test[['class2']] <- "CE"
metabodata_info_test = metabodata_info_test %>%
  mutate(class2 = ifelse(class %in% c('TAG','DAG','MAG'),'TAG, DAG, or MAG',class2)) %>%
  mutate(class2 = ifelse(class %in% c('PC','PE','PI'),'PC, PE, or PI',class2)) %>%
  mutate(class2 = ifelse(class %in% c('LPC','LPE'),'LPC or LPE',class2)) %>%
  mutate(class2 = ifelse(class %in% c('SM','CER',"DCER","HCER","LCER"),'SM or ceramides',class2)) %>%
  mutate(class2 = ifelse(class %in% c('Amino acids','Glycolysis related metabolites','Inflammation','Ketone bodies'),'Nightingale',class2))
table(metabodata_info_test$class,metabodata_info_test$class2,useNA = 'a')
class2.levels = c("CE", 
                  "TAG, DAG, or MAG", 
                  "PC, PE, or PI",
                  "LPC or LPE", 
                  "SM or ceramides",
                  "Nightingale")
metabodata_info_test = metabodata_info_test %>%
  mutate(class2 = factor(class2,levels=class2.levels)) 
metabodata_info_test = metabodata_info_test %>% rename(lipid_class = class2)
metabodata_info_test = metabodata_info_test %>% relocate(lipid_class,.after = "class")
metabodata_info_test = metabodata_info_test %>% 
  mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))
head(metabodata_info_test)

metabodata_test = metabodata
names(metabodata_test) = str_remove(names(metabodata_test),"_nmol/mL")
metabodata_test = subset(metabodata_test,select=c('uniqueID',testing.lipids2))

mdat1 = subset(
  metabodata_test, 
  select=c('uniqueID',metabodata_info_test %>% filter(lipid_class=="CE") %>% pull(metabolite_id)))
mdat1.names = names(mdat1)[-1]
mdat1 = mdat1 %>% left_join(covdata,join_by(uniqueID))
table(mdat1$Sex)
mi=1
mdat1.cleaned = remove.outliers_grubbs(mdat1,varnames = mdat1.names) 
mdat1.cleaned$counts.NA
hist(mdat1$`CE(22:2)`)
hist(mdat1.cleaned$cleaed_data$`CE(22:2)`)

# mdat1.cleaned$cleaed_data %>% 
(mdat1 %>% 
  ggplot(aes(x=age,y=inormal(`CE(22:2)`),color=Sex)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method="gam")
)+

  (mdat1.cleaned$cleaed_data %>% 
  ggplot(aes(x=age,y=`CE(22:2)`,color=Sex)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(method="gam"))

# generate ado/adults volcano plots ----

##
plot.df.volcano_main_NormWM = add_class2(volcano_main_NormWM$data)
str(plot.df.volcano_main_NormWM$volcano_main_all)
plot.df.volcano_main_NormWM$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormWM$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_NormGM = add_class2(volcano_main_NormGM$data)
plot.df.volcano_main_NormGM$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormGM$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_NormGM$volcano_main_all)
table(plot.df.volcano_main_NormGM$volcano_main_all$class)

plot.df.volcano_main_NormWM_adults = add_class2(volcano_main_NormWM_adults$data)
plot.df.volcano_main_NormWM_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormWM_adults$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_NormWM_adults$volcano_main_all)
table(plot.df.volcano_main_NormWM_adults$volcano_main_all$`lipid class`)

plot.df.volcano_main_NormGM_adults = add_class2(volcano_main_NormGM_adults$data)
plot.df.volcano_main_NormGM_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormGM_adults$volcano_main_all$metabolite_id,"_nmol/mL")
table(plot.df.volcano_main_NormGM_adults$volcano_main_all$`lipid class`)

plot.df.volcano_main_WM.vol = add_class2(volcano_main_VolumeWM$data)
str(plot.df.volcano_main_WM.vol$volcano_main_all)
plot.df.volcano_main_WM.vol$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_WM.vol$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_GM.vol = add_class2(volcano_main_VolumeGM$data)
plot.df.volcano_main_GM.vol$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_GM.vol$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_GM.vol$volcano_main_all)
table(plot.df.volcano_main_GM.vol$volcano_main_all$class)

plot.df.volcano_main_WM.vol_adults = add_class2(volcano_main_VolumeWM_adults$data)
str(plot.df.volcano_main_WM.vol_adults$volcano_main_all)
plot.df.volcano_main_WM.vol_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_WM.vol_adults$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_GM.vol_adults = add_class2(volcano_main_VolumeGM_adults$data)
plot.df.volcano_main_GM.vol_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_GM.vol_adults$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_GM.vol_adults$volcano_main_all)
table(plot.df.volcano_main_GM.vol_adults$volcano_main_all$class)

p.volume = c(plot.df.volcano_main_WM.vol$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_GM.vol$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_WM.vol_adults$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_GM.vol_adults$volcano_main_all %>% pull(pvalue))
p.T1wSI = c(plot.df.volcano_main_NormWM$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_NormGM$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_NormWM_adults$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_NormGM_adults$volcano_main_all %>% pull(pvalue))

fdr.p.T1wSI = p.adjust(p.T1wSI,method='fdr')
max(p.T1wSI[fdr.p.T1wSI<0.05])
hist(p.T1wSI,breaks=100);abline(v=max(p.T1wSI[fdr.p.T1wSI<0.05]),col="red")#0.01103552

fdr.p.volume = p.adjust(p.volume,method='fdr')
max(p.volume[fdr.p.volume<0.05])
hist(p.volume,breaks=100);abline(v=max(p.volume[fdr.p.volume<0.05]),col="red")#0.007237435

p.all = c(p.volume,p.T1wSI)
fdr.p = p.adjust(p.all,method='fdr')
hist(p.all,breaks=100);abline(v=max(p.all[fdr.p<0.05]),col="red")#p = 0.009080251

## WM and GM volume - volcano plots (adolescents and adults) ----
point.size = 1.75
volcano_WM.vol = plot.df.volcano_main_WM.vol$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_GM.vol$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_GM.vol$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_GM.vol$cols,alpha=0.5))
  ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.vol

volcano_GM.vol = plot.df.volcano_main_GM.vol$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_GM.vol$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_GM.vol$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_GM.vol$cols,alpha=0.5))
  ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_GM.vol') 
volcano_GM.vol

volcano_WM.vol_adults = plot.df.volcano_main_WM.vol_adults$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%#830
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  # scale_color_manual(
  #   values=plot.df.volcano_main_NormWM_adults$cols,
  #   name="Lipid classes",
  #   labels = plot.df.volcano_main_NormWM_adults$class2.labels,
  #   guide = guide_legend(
  #     override.aes = list(color=plot.df.volcano_main_NormWM_adults$cols,alpha=0.5))
  # ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol_adults')
volcano_WM.vol_adults

volcano_GM.vol_adults = plot.df.volcano_main_GM.vol_adults$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_GM.vol_adults$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_GM.vol_adults$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_GM.vol_adults$cols,alpha=0.5))
  ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_GM.vol_adults')#
volcano_GM.vol_adults

volcano_WM.vol + volcano_GM.vol + volcano_WM.vol_adults + volcano_GM.vol_adults +
  plot_layout(guides='collect',axis_titles = "collect",nrow=2,byrow = T) &
  theme(legend.position = 'none') &
  coord_cartesian(ylim=c(0,12.5),xlim=c(-0.225,0.225)) 

ggsave("results/volcanoplot_Volume_sex_combined_unadjBMI_ados_adults.png",
       width=11.5,height=8,units='in',dpi=300)

## T1wSI volcano plots (adolescents and adults) ----
point.size = 1.75
volcano_NormWM = plot.df.volcano_main_NormWM$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_NormWM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_NormWM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_NormWM$cols,alpha=0.5))
  ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+  ggtitle("NormWM (adolescents)")
volcano_NormWM

volcano_NormGM = plot.df.volcano_main_NormGM$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_NormGM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_NormGM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_NormGM$cols,alpha=0.5))
  ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+  ggtitle("NormGM (adolescents)")
volcano_NormGM

volcano_NormWM_adults = plot.df.volcano_main_NormWM_adults$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%#830
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  # scale_color_manual(
  #   values=plot.df.volcano_main_NormWM_adults$cols,
  #   name="Lipid classes",
  #   labels = plot.df.volcano_main_NormWM_adults$class2.labels,
  #   guide = guide_legend(
  #     override.aes = list(color=plot.df.volcano_main_NormWM_adults$cols,alpha=0.5))
  # ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+  ggtitle("NormWM (adults)")
volcano_NormWM_adults

volcano_NormGM_adults = plot.df.volcano_main_NormGM_adults$volcano_main_all %>%
  filter(metabolite_id %in% testing.lipids2) %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_NormGM_adults$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_NormGM_adults$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_NormGM_adults$cols,alpha=0.5))
  ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  geom_hline(yintercept = -log10(max(p.all[fdr.p<0.05])), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+  ggtitle("NormGM (adults)")
volcano_NormGM_adults

volcano_NormWM+volcano_NormGM+volcano_NormWM_adults+volcano_NormGM_adults +
  plot_layout(guides='collect',axis_titles = "collect",nrow=2,byrow = T) &
  theme(legend.position = 'none') &
  coord_cartesian(ylim=c(0,12.5),xlim=c(-0.225,0.225)) 
ggsave("results/volcanoplot_sex_combined_unadjBMI_ados_adults.png",
       width=11.5,height=8,units='in',dpi=300)

## volcano 3 phenotypes ----
volcano_main_3brain = plot.df.volcano_main_lobarVolWM$volcano_main_all %>% 
  # filter(metabolite_id %in% siglipids) %>%
  select(metabolite_id,class,class2,`lipid class`,beta,pvalue) %>%
  rename(beta.Vol = beta, pvalue.Vol = pvalue) %>% 
  left_join(
    plot.df.volcano_main_NormWM$volcano_main_all %>% 
      # filter(metabolite_id %in% siglipids) %>%
      select(metabolite_id,class,class2,`lipid class`,beta,pvalue) %>%
      rename(beta.nT1wSI = beta, pvalue.nT1wSI = pvalue)
  ) %>%
  left_join(
    plot.df.volcano_main_MTR_WM$volcano_main_all %>% 
      # filter(metabolite_id %in% siglipids) %>%
      select(metabolite_id,class,class2,`lipid class`,beta,pvalue) %>%
      rename(beta.MTR = beta, pvalue.MTR = pvalue)
  ) %>%
  mutate(p.multi = pvalue.Vol*pvalue.nT1wSI*pvalue.MTR) %>%
  arrange(p.multi)
head(volcano_main_3brain)

# locuszoom plots ----
get_locuszoom_range = function(chr,b0,b1) {
  b0=b0-250000
  b1=b1+250000
  ret = c(b0,b1)/10^6
  ret = c(chr,ret)
  ret
}
## rs650612(DHCR24):
LZ_rs650612 = get_locuszoom_range(1,55119702,55377952)

## rs174566(FADS1):
LZ_rs174566 = get_locuszoom_range(11,61542006,61624181)

# [plot]_volcano-plots (adolescents) ----
rm(volcano_main,volcano_main_all,p,top_main)

## generate plot ----
point.size = 0.95
top_main <- volcano_main_3brain %>% 
  filter(beta.Vol<0, beta.nT1wSI>0, beta.MTR>0) %>% slice(1:10) %>% 
  bind_rows(
    volcano_main_3brain %>% filter(beta.Vol>0, beta.nT1wSI<0, beta.MTR<0) %>% slice(1:10)
  )

top_main_lobarVolWM <- plot.df.volcano_main_lobarVolWM$volcano_main_all %>% 
  filter(metabolite_id %in% top_main$metabolite_id) %>%
  mutate(label.col = ifelse(!metabolite_id %in% siglipids,"darkred","grey30"))
top_main_NormWM <- plot.df.volcano_main_NormWM$volcano_main_all %>% 
  filter(metabolite_id %in% top_main$metabolite_id) %>%
  mutate(label.col = ifelse(!metabolite_id %in% siglipids,"darkred","grey30"))
top_main_MTR_WM <- plot.df.volcano_main_MTR_WM$volcano_main_all %>% 
  filter(metabolite_id %in% top_main$metabolite_id) %>%
  mutate(label.col = ifelse(!metabolite_id %in% siglipids,"darkred","grey30"))

### Volume ----
volcano_lobarVolWM = plot.df.volcano_main_lobarVolWM$volcano_main_all%>%
  ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(
    values=plot.df.volcano_main_lobarVolWM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_lobarVolWM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_lobarVolWM$cols,alpha=0.5))
  ) + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  # geom_hline(yintercept = -log10(0.1), color="darkgrey") + 
  geom_text_repel(
    data = top_main_lobarVolWM,
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha(top_main_lobarVolWM$label.col,0.25),
    segment.size = 0.2,
    # segment.color = scales::alpha("lightgray",0.5),
    color = top_main_lobarVolWM$label.col,#'black',
    #fill = scales::alpha('white',0.5),
    size = 2.5) #theme_bw()

### NormWM ----
volcano_NormWM = plot.df.volcano_main_NormWM$volcano_main_all%>%
  ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(
    values=plot.df.volcano_main_NormWM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_NormWM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_NormWM$cols,alpha=0.5))
  ) +
  geom_vline(xintercept = 0, color="darkgrey") +
  # geom_hline(yintercept = -log10(0.1), color="darkgrey") + 
  geom_text_repel(
    data = top_main_NormWM,
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha(top_main_NormWM$label.col,0.25),
    segment.size = 0.2,
    # segment.color = scales::alpha("lightgray",0.5),
    color = top_main_NormWM$label.col,#'black',
    #fill = scales::alpha('white',0.5),
    size = 2.5) #theme_bw()

### MTR ----
volcano_MTR_WM = plot.df.volcano_main_MTR_WM$volcano_main_all %>%
  ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(
    values=plot.df.volcano_main_MTR_WM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_MTR_WM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_MTR_WM$cols,alpha=0.5))
  ) +
  geom_vline(xintercept = 0, color="darkgrey") +
  # geom_hline(yintercept = -log10(0.1), color="darkgrey") + 
  geom_text_repel(
    data = top_main_MTR_WM,
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha(top_main_MTR_WM$label.col,0.25),
    segment.size = 0.2,
    # segment.color = scales::alpha("lightgray",0.5),
    color = top_main_MTR_WM$label.col,#'black',
    #fill = scales::alpha('white',0.5),
    size = 2.5)

### put all three volcano plots ----
volcano_lobarVolWM + volcano_NormWM + volcano_MTR_WM +
  plot_layout(guides='collect',axis_titles = "collect") &
  theme(legend.position = 'none') &
  coord_cartesian(ylim=c(0,12.5),xlim=c(-0.225,0.225)) &
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=11),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) 

### generate png file ----
ggsave("results/volcanoplot_sex_combined_unadjBMI.png",
       width=11.5,height=3.5,units='in',dpi=300)

# [plot]_scatter plots: assoc_beta ----
point.size=0.95
## Volume vs T1wSI ----
top_main_lobarVolWM$label.col = 'gray30'
scatter_lobarVolWM_T1wSI = volcano_main_3brain %>%
  ggplot(aes(x=beta.Vol,y=beta.nT1wSI, color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) +
  scale_color_manual(
    values=plot.df.volcano_main_lobarVolWM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_lobarVolWM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_lobarVolWM$cols,alpha=0.5))
  ) + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_text_repel(
    data = volcano_main_3brain %>% filter(metabolite_id %in% top_main$metabolite_id),
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha(top_main_lobarVolWM$label.col,0.25),
    segment.size = 0.2,
    # segment.color = scales::alpha("lightgray",0.5),
    color = top_main_lobarVolWM$label.col,#'black',
    #fill = scales::alpha('white',0.5),
    size = 2.25)+ #theme_bw()
xlab(bquote(hat(beta) ~ "_WM Volume"))+
ylab(bquote(hat(beta) ~ "_WM nT1wSI"))#+  ggpubr::stat_cor()
scatter_lobarVolWM_T1wSI

## Volume vs. MTR ----
top_main_lobarVolWM$label.col = 'gray30'
scatter_lobarVolWM_MTR = volcano_main_3brain %>%
  # ggplot(aes(x=beta.Vol,y=beta.MTR)) + 
  ggplot(aes(x=beta.Vol,y=beta.MTR, color=`lipid class`)) +
  geom_point(shape=19,alpha=0.4,size=point.size) +
  scale_color_manual(
    values=plot.df.volcano_main_lobarVolWM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_lobarVolWM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_lobarVolWM$cols,alpha=0.5))
  ) + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_text_repel(
    data = volcano_main_3brain %>% filter(metabolite_id %in% top_main$metabolite_id),
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha(top_main_lobarVolWM$label.col,0.25),
    segment.size = 0.2,
    # segment.color = scales::alpha("lightgray",0.5),
    color = top_main_lobarVolWM$label.col,#'black',
    #fill = scales::alpha('white',0.5),
    size = 2.25)+ #theme_bw()
  xlab(bquote(hat(beta) ~ "_WM Volume"))+
  ylab(bquote(hat(beta) ~ "_WM MTR")) #+ ggpubr::stat_cor()
scatter_lobarVolWM_MTR

## T1wSI vs. MTR ----
scatter_T1wSI_MTR = volcano_main_3brain %>%
  ggplot(aes(x=beta.nT1wSI,y=beta.MTR, color=`lipid class`)) +
  # ggplot(aes(x=beta.nT1wSI,y=beta.MTR)) + ggpubr::stat_cor() + 
  geom_point(shape=19,alpha=0.4,size=point.size) + 
  scale_color_manual(
    values=plot.df.volcano_main_lobarVolWM$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_lobarVolWM$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_lobarVolWM$cols,alpha=0.5))
  ) + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_text_repel(
    data = volcano_main_3brain %>% filter(metabolite_id %in% top_main$metabolite_id),
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha(top_main_lobarVolWM$label.col,0.25),
    segment.size = 0.2,
    # segment.color = scales::alpha("lightgray",0.5),
    color = top_main_lobarVolWM$label.col,#'black',
    #fill = scales::alpha('white',0.5),
    size = 2.25)+ #theme_bw()
  xlab(bquote(hat(beta) ~ "_WM nT1wSI"))+
  ylab(bquote(hat(beta) ~ "_WM MTR"))

## generate scatter plots ----
scatter_lobarVolWM_T1wSI + scatter_lobarVolWM_MTR + scatter_T1wSI_MTR +
  plot_layout(guides='collect') &
  theme(legend.position = 'none') &
  coord_cartesian(ylim=c(-0.225,0.35),xlim=c(-0.225,0.225)) &
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=9),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) #&
  # geom_abline(slope=1,intercept = 0, color='grey')
ggsave(paste0("results/scatterplot_sex_combined_unadjBMI_",Sys.Date(),".png"),
       width=11.5,height=3.5,units='in',dpi=300)

cor(volcano_main_3brain %>%
      select(beta.Vol,beta.nT1wSI,beta.MTR))

# [analysis]: Assoc_cognition_vs_SNPs ----
replace.neg.wi.NA = function(x,negnum){
  ind = !is.na(x) & x == negnum
  x[ind] <- NA
  x
}

d_cognition = fread("~/Downloads/jean_shin_dataanalysis_completesyssubjects_rev1_20250528_121659.txt")
names(d_cognition) = str_remove(names(d_cognition),'ADO1-CognitionZscored.')
snp_dosages = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_dosages.csv")
snp_dosages_info = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_snps_info.csv")
names(snp_dosages)
head(snp_dosages_info)
head(snp_dosages)
geno_pcs = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/info_used_for_gwas.txt/sys_ado1_covariates_to_be_reordered.txt")
head(geno_pcs)

d_cognition = data.frame(d_cognition)
d_cognition[,-c(1:3)] = apply(d_cognition[,-c(1:3)],2,as.numeric)
d_cognition[,-c(1:3)] = apply(d_cognition[,-c(1:3)],2,replace.neg.wi.NA,negnum=-999)
d_cognition[,-c(1:3)] = apply(d_cognition[,-c(1:3)],2,replace.neg.wi.NA,negnum=-888)
describe(d_cognition) %>% select(min)

analdat = d_cognition %>% 
  left_join(snp_dosages %>% select(sub_id,rs650612,rs174566)) %>%
  left_join(geno_pcs %>% select(uniqueID,age,sex,genoPC1,genoPC2,genoPC3,genoPC4)) %>%
  mutate(prs_13 = -0.27*rs650612-0.26*rs174566,
         prs_unweighted = rs650612+rs174566)
head(analdat)
ynames = names(d_cognition)[-c(1:3)]
xnames = c("rs650612","rs174566","prs_13","prs_unweighted")
analdat.cleaned = remove.outliers_grubbs(analdat,varnames=ynames)
analdat.cleaned$counts.NA %>% mutate(diff = after-before)

# PCA of cognitive traits (2025-05-28) ----
library("factoextra")
head(analdat.cleaned$cleaed_data)
y1=which(names(analdat.cleaned$cleaed_data)=='RuAutoSpeed')
y2=which(names(analdat.cleaned$cleaed_data)=='TapNonDominantHandRep')

analdat.cleaned$cleaed_data = data.frame(analdat.cleaned$cleaed_data)
rownames(analdat.cleaned$cleaed_data) = analdat.cleaned$cleaed_data$uniqueID
pca_cog = PCA(na.omit(analdat.cleaned$cleaed_data[,y1:y2]),ncp=10)
# fviz_eig(pca_cog, addlabels = TRUE, ylim = c(0, 30))
# dim(pca_cog$ind$coord[,1:2])#792

# remove 8 individuals with missing scores for >=45 traits
table(apply(apply(analdat.cleaned$cleaed_data[,y1:y2],2,is.na),1,sum))
rm.IDs = analdat.cleaned$cleaed_data %>% slice(which(apply(apply(analdat.cleaned$cleaed_data[,y1:y2],2,is.na),1,sum)>=45)) %>% pull(uniqueID)

analdat.cleaned.data = analdat.cleaned$cleaed_data %>% 
  filter(!uniqueID %in% rm.IDs) 
analdat.cleaned.data = analdat.cleaned.data %>%
  filter(uniqueID %in% (analdat.cleaned.data %>% filter(!is.na(GFactorPC1))%>% pull(uniqueID) ))

analdat.cleaned.data = data.frame(analdat.cleaned.data)
rownames(analdat.cleaned.data) = analdat.cleaned.data$uniqueID
impute.dat = missMDA::imputePCA(analdat.cleaned.data[,y1:y2],ncp=10)
pca_cog_imputed = FactoMineR::PCA(impute.dat$completeObs,ncp=10)
fviz_eig(pca_cog_imputed, addlabels = TRUE, ylim = c(0, 30))

analdat.cleaned$cleaed_data = analdat.cleaned$cleaed_data %>% 
  left_join(data.frame(uniqueID = rownames(pca_cog_imputed$ind$coord[,1:3]), 
                       pca_cog_imputed$ind$coord[,1:3]))
analdat.cleaned$cleaed_data %>%
  ggplot(aes(x=GFactorPC1,y=Dim.1)) + geom_point() + 
  geom_abline(slope=1,intercept = 0)

## assoc_tests ----
ynames = c(ynames,paste("Dim",1:3,sep="."));
ynames = unique(ynames)
res_all=c()
for(xj in xnames){
  n =c(); assoc_res = c(); 
  for(yi in ynames){
    di = analdat
    di[['y']] = analdat.cleaned$cleaed_data[[yi]]
    di[['x']] = analdat.cleaned$cleaed_data[[xj]]
    fiti = lmer(inormal(y)~ x + age + sex + genoPC1 + genoPC2+ genoPC3+ genoPC4 + (1|fam_id),
                data=di,na.action = na.exclude)
    # lmi = lm(inormal(y)~ x + age + sex + genoPC1 + genoPC2+ genoPC3+ genoPC4 ,
    #            data=di,na.action = na.exclude)
    # plot(lmi)
    assoc_res = rbind(assoc_res,summary(fiti)$coef["x",])
    n = c(n,nobs(fiti))
  }
  res_all = rbind(res_all,data.table(yname=ynames,rsID = xj,assoc_res,n=n))
}
dim(res_all)
head(res_all)
write_tsv(res_all,"results/assoc_cognition_DHCR24_FADS1_snps.txt")

## forest plot ----
pos <- position_dodge(width=0.75)#
alpha.scale = 0.95
res_all = fread("results/assoc_cognition_DHCR24_FADS1_snps.txt") %>% 
  rename(beta=Estimate,se = `Std. Error`, pvalue=`Pr(>|t|)`) %>%
  mutate(L95CI = beta-1.96*se,U95CI = beta+1.96*se)
ylevel = res_all%>% filter(rsID=="prs_unweighted") %>% arrange(beta) %>% pull(yname)
res_all = res_all %>% mutate(yname = factor(yname,levels=rev(ylevel))) %>% arrange(yname)
res_all = res_all %>%
  mutate(signif = ifelse(pvalue<0.05,"p<0.05",'p>=0.05'))
res_all$rsID = str_replace(res_all$rsID,"rs174566","rs174566_G_FADS1") 
res_all$rsID = str_replace(res_all$rsID,"rs650612","rs650612_C_DHCR24") 

res_all %>% filter(yname!="Dim.1") %>%
  ggplot(aes(y=beta,x=yname ,ymin=L95CI,ymax=U95CI,group=rsID,col=rsID,
                       shape=signif,fill=signif)) +
  geom_point( position=pos,alpha=alpha.scale) + 
  scale_shape_manual(values=c("p<0.05"=19,"p>=0.05"=21)) + 
  scale_fill_manual(values=c("p<0.05"="black","p>=0.05"="white")) + 
  geom_path(aes(group=rsID), position=pos,alpha=alpha.scale) + 
  geom_errorbar(position=pos, width=0.25,linewidth=0.5,alpha=0.4) + 
  # facet_grid(cols=vars(rsID)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45,hjust=0,vjust=0.5))+
  xlab(NULL) + ylab("Estimate (95%CI)") +
  guides(col = guide_legend(reverse=T))#700*725/1440*510
ggsave('results/forest_cognition_traits_vs_DHCR24_FADS1_SNPs_2025-05-28.png',
       width=14.4,height=5.1,units='in',dpi=300)

## wide table ----------------------------------------------------
sort(ynames)
res_all = res_all %>% filter(yname!="Dim.1") %>% mutate(fdr.p = p.adjust(pvalue,method="fdr"))
res_all %>% arrange(pvalue)
res_all %>% filter(rsID == "prs_unweighted") %>% arrange(yname) %>% select(yname,beta,se,pvalue)

res_all_wide = res_all %>% 
  filter(rsID == "prs_unweighted") %>%
  select(yname,beta,se,pvalue) %>%
  rename(beta.combi = beta, se.combi=se, pvalue.combi=pvalue) %>%
  #
  left_join(
    res_all %>% 
      filter(rsID == "prs_13") %>% 
      select(yname,beta,se,pvalue) %>%
      rename(beta.prs_13 = beta, se.prs_13=se, pvalue.prs_13=pvalue)
  ) %>%
  left_join(
    res_all %>% 
      filter(rsID == 'rs650612_C_DHCR24') %>% 
      select(yname,beta,se,pvalue) %>%
      rename(beta.rs650612 = beta, se.rs650612=se, pvalue.rs650612=pvalue), join_by(yname)
  ) %>%
  left_join(
    res_all %>% 
      filter(rsID == 'rs174566_G_FADS1') %>% 
      select(yname,beta,se,pvalue) %>%
      rename(beta.rs174566 = beta, se.rs174566=se, pvalue.rs174566=pvalue), join_by(yname)
  )

write_tsv(res_all_wide,"results/assoc_cognition_DHCR24_FADS1_snps_wide.txt")
head(res_all %>% arrange(pvalue) %>% slice(1:10))

round((data.frame(pca_cog_imputed$var$coord) %>% filter(abs(Dim.2)>=0.4)) %>% arrange(Dim.2) %>% select(1:3),2)
## 
readr::write_tsv(data.frame(rownames(data.frame(pca_cog_imputed$var$coord)),pca_cog_imputed$var$coord),"~/Downloads/tab.txt")

# Associations between SNPs and process speed (AutomaticDetectionSpeed) ----
snp_dosages = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_dosages.csv")
snp_dosages_info = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_snps_info.csv")
names(snp_dosages)
head(snp_dosages_info)
head(snp_dosages)

d_process = fread("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Temporary/Manuscript_Eeva/Seladin1/Proc speed Jean.csv")
head(d_process)
names(d_process) = str_remove(names(d_process),'ADO1-E08f.')

geno_pcs = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/info_used_for_gwas.txt/sys_ado1_covariates_to_be_reordered.txt")
head(geno_pcs)

describe(d_process)

d_rs588709 = fread("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Temporary/Manuscript_Eeva/Seladin1/rs588709_mldoes.csv")
head(d_rs588709)
head(snp_dosages)
snp_dosages = snp_dosages%>% left_join(d_rs588709)
head(snp_dosages)
cor(snp_dosages[,-c(1:2)],use="p")

analdat = d_process %>% 
  left_join(snp_dosages %>% select(sub_id,rs650612,rs174566)) %>%
  left_join(d_rs588709) %>%
  left_join(geno_pcs %>% select(uniqueID,age,sex,genoPC1,genoPC2,genoPC3,genoPC4)) %>%
  mutate(prs_13 = -0.27*rs650612-0.26*rs174566,
         prs_unweighted = rs650612+rs174566)
head(analdat)
analdat = analdat %>% 
  filter(uniqueID %in% metabodata_adj$uniqueID) %>%
  mutate(y = inormal(AutomaticDetectionSpeed))

xnames = c('rs588709',"rs650612","rs174566","prs_13","prs_unweighted")
fit_res_all= c()
for(xi in 1:length(xnames)){
  df = data.frame(analdat)%>%
    bind_cols(x=analdat[[xnames[xi]]])
  head(df)
  fiti = lmer(y~ x + age * sex + genoPC1 + genoPC2+ genoPC3+ genoPC4 + (1|fam_id),
              data = df, na.action = na.exclude)
  fit_res_all = rbind(fit_res_all,summary(fiti)$coef['x',])
}
fit_res_all = data.table(rsID = xnames,fit_res_all) 

#check
