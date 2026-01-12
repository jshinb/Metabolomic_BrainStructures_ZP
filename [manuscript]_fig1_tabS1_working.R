# Generated output files ----
# c(
#   "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/SupplementaryTables_2025-11-27.xlsx",
#   "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/Manuscript_21-11-2025_js.docx",
#   "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/Assoc_lipids_MTR_2025-11-27.pptx"
# )
# modified: Dec.1, 2025 - 
#* remove GM results
#* updated the fdr-adjusted significance level (for 831*6 tests)
#* generated new volcano plots 
#* updated the number
#* 
#* To-Do: adults [sensitivity] vs adolescents [primary]
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
generate_volcano_plot = function(df_plot, df, trait_plot, generation_plot){
  
  p = df_plot %>%
    filter(trait == trait_plot, generation==generation_plot) %>%
    ggplot(aes(x=beta,y=-log10(pvalue),
               shape=fdr.signif,
               fill=fdr.signif,
               # color=fdr.signif)) + 
               color=`lipid class`)) + 
    geom_point(alpha=0.25,size=point.size,stroke = 2) +
    scale_shape_manual(values=c(19,21))+ 
    # scale_color_manual(values=c("#073B4C",'grey'))+
    scale_color_manual(
      values=df$cols,
      name="Lipid classes",
      labels = df$class2.labels,
      guide = guide_legend(
        override.aes = list(color=df$cols,alpha=0.5))
    ) + 
    scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
    geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
    geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
  p
}

# load significant lipids ----

## ado data (found typo on Nov 14, 2025) ----
load("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
table(volcano_main$data$N)#736-916
volcano_main_NormWM = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_NormGM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
table(volcano_main$data$N)#732 - 912
volcano_main_NormGM = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.Rdata")
table(volcano_main$data$N)#669-826
volcano_main_VolumeWM = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_ados_unadjBMI_2025-05-26.Rdata")
table(volcano_main$data$N)#669-826
volcano_main_VolumeGM = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_LobarWM_Z MTR_sex-combined_ados_unadjBMI_2025-04-08.Rdata")
table(volcano_main$data$N)#604-745
volcano_main_MTR_WM = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_LobarGM_Z_sex-combined_ados_unadjBMI_2025-11-25.Rdata")
table(volcano_main$data$N)#599-738
volcano_main_MTR_GM = volcano_main$data
rm(volcano_main)

## adult data ----
load("results/volcanoplot_all_NormWM_sex-combined_adults_unadjBMI_2025-01-09.Rdata")
table(volcano_main$data$N)#454-559
volcano_main_NormWM_adults = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_NormGM_sex-combined_adults_unadjBMI_2025-01-09.Rdata")
table(volcano_main$data$N)#454-559
volcano_main_NormGM_adults = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
table(volcano_main$data$N)#267-361
volcano_main_VolumeWM_adults = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
table(volcano_main$data$N)#267-361
volcano_main_VolumeGM_adults = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_wm_sex-combined_adults_unadjBMI_2025-11-25.Rdata")
table(volcano_main$data$N)#262-355
volcano_main_MTR_WM_adults = volcano_main$data
rm(volcano_main)

load("results/volcanoplot_all_gm_sex-combined_adults_unadjBMI_2025-11-25.Rdata")
table(volcano_main$data$N)#258-351
volcano_main_MTR_GM_adults = volcano_main$data
rm(volcano_main)

# volcano plot for all generation ----
## Volume ----
plot.df.volcano_main_WM.vol = add_class2(volcano_main_VolumeWM)
str(plot.df.volcano_main_WM.vol$volcano_main_all)
plot.df.volcano_main_WM.vol$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_WM.vol$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_GM.vol = add_class2(volcano_main_VolumeGM)
plot.df.volcano_main_GM.vol$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_GM.vol$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_GM.vol$volcano_main_all)
table(plot.df.volcano_main_GM.vol$volcano_main_all$class)

plot.df.volcano_main_WM.vol_adults = add_class2(volcano_main_VolumeWM_adults)
str(plot.df.volcano_main_WM.vol_adults$volcano_main_all)
plot.df.volcano_main_WM.vol_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_WM.vol_adults$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_GM.vol_adults = add_class2(volcano_main_VolumeGM_adults)
plot.df.volcano_main_GM.vol_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_GM.vol_adults$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_GM.vol_adults$volcano_main_all)
table(plot.df.volcano_main_GM.vol_adults$volcano_main_all$class)

## T1wSI ----
plot.df.volcano_main_NormWM = add_class2(volcano_main_NormWM)
str(plot.df.volcano_main_NormWM$volcano_main_all)
plot.df.volcano_main_NormWM$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormWM$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_NormGM = add_class2(volcano_main_NormGM)
plot.df.volcano_main_NormGM$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormGM$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_NormGM$volcano_main_all)
table(plot.df.volcano_main_NormGM$volcano_main_all$class)

plot.df.volcano_main_NormWM_adults = add_class2(volcano_main_NormWM_adults)
plot.df.volcano_main_NormWM_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormWM_adults$volcano_main_all$metabolite_id,"_nmol/mL")
str(plot.df.volcano_main_NormWM_adults$volcano_main_all)
table(plot.df.volcano_main_NormWM_adults$volcano_main_all$`lipid class`,useNA='a')

plot.df.volcano_main_NormGM_adults = add_class2(volcano_main_NormGM_adults)
plot.df.volcano_main_NormGM_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_NormGM_adults$volcano_main_all$metabolite_id,"_nmol/mL")
table(plot.df.volcano_main_NormGM_adults$volcano_main_all$`lipid class`,useNA='a')

## MTR ----
plot.df.volcano_main_WM.MTR = add_class2(volcano_main_MTR_WM)
str(plot.df.volcano_main_WM.MTR$volcano_main_all)
plot.df.volcano_main_WM.MTR$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_WM.MTR$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_GM.MTR = add_class2(volcano_main_MTR_GM)
str(plot.df.volcano_main_GM.MTR$volcano_main_all)
plot.df.volcano_main_GM.MTR$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_GM.MTR$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_WM.MTR_adults = add_class2(volcano_main_MTR_WM_adults)
str(plot.df.volcano_main_WM.MTR_adults$volcano_main_all)
plot.df.volcano_main_WM.MTR_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_WM.vol_adults$volcano_main_all$metabolite_id,"_nmol/mL")

plot.df.volcano_main_GM.MTR_adults = add_class2(volcano_main_MTR_GM_adults)
str(plot.df.volcano_main_GM.MTR_adults$volcano_main_all)
plot.df.volcano_main_GM.MTR_adults$volcano_main_all$metabolite_id = 
  str_remove(plot.df.volcano_main_GM.vol_adults$volcano_main_all$metabolite_id,"_nmol/mL")

## p-values ----
p.volume = c(plot.df.volcano_main_WM.vol$volcano_main_all %>% pull(pvalue),
             # plot.df.volcano_main_GM.vol$volcano_main_all %>% pull(pvalue),
             plot.df.volcano_main_WM.vol_adults$volcano_main_all %>% pull(pvalue))#,
             # plot.df.volcano_main_GM.vol_adults$volcano_main_all %>% pull(pvalue))

p.T1wSI = c(plot.df.volcano_main_NormWM$volcano_main_all %>% pull(pvalue),
            #plot.df.volcano_main_NormGM$volcano_main_all %>% pull(pvalue),
            plot.df.volcano_main_NormWM_adults$volcano_main_all %>% pull(pvalue))#,
            #plot.df.volcano_main_NormGM_adults$volcano_main_all %>% pull(pvalue))

p.MTR = c(plot.df.volcano_main_WM.MTR$volcano_main_all %>% pull(pvalue),
          #plot.df.volcano_main_GM.MTR$volcano_main_all %>% pull(pvalue),
          plot.df.volcano_main_WM.MTR_adults$volcano_main_all %>% pull(pvalue))#,
          #plot.df.volcano_main_GM.MTR_adults$volcano_main_all %>% pull(pvalue))

p.all = c(p.volume, p.T1wSI, p.MTR)
fdr.p.all = p.adjust(p.all, method='fdr')
signif.level = max(p.all[fdr.p.all<0.05])#0.01889204
hist(p.all,breaks=100);abline(v=max(p.all[fdr.p.all<0.05]),col="red")#0.01103552#0.01347401

# select lipids tested in both ados and adults ----
testing.lipids2 = intersect(
  plot.df.volcano_main_WM.vol$volcano_main_all$metabolite_id,
  plot.df.volcano_main_WM.vol_adults$volcano_main_all$metabolite_id
)

# augment all association results ----
plot.df.volcano_main = rbind(
  plot.df.volcano_main_WM.vol$volcano_main_all,
  plot.df.volcano_main_WM.vol_adults$volcano_main_all,
  #plot.df.volcano_main_GM.vol$volcano_main_all,
  #plot.df.volcano_main_GM.vol_adults$volcano_main_all,
  plot.df.volcano_main_NormWM$volcano_main_all,
  plot.df.volcano_main_NormWM_adults$volcano_main_all,
  #plot.df.volcano_main_NormGM$volcano_main_all,
  #plot.df.volcano_main_NormGM_adults$volcano_main_all,
  plot.df.volcano_main_WM.MTR$volcano_main_all,
  plot.df.volcano_main_WM.MTR_adults$volcano_main_all#,
  #plot.df.volcano_main_GM.MTR$volcano_main_all,
  #plot.df.volcano_main_GM.MTR_adults$volcano_main_all
)
table(plot.df.volcano_main %>% select(generation, trait))

plot.df.volcano_main = plot.df.volcano_main %>%
  filter(metabolite_id %in% testing.lipids2)
plot.df.volcano_main = plot.df.volcano_main %>% 
  mutate(trait = ifelse(trait == 'wm_fam',"MTR.WM_fam",trait))
plot.df.volcano_main = plot.df.volcano_main %>% 
  mutate(trait = ifelse(trait == 'gm_fam',"MTR.GM_fam",trait))

plot.df.volcano_main = plot.df.volcano_main %>% 
  mutate(trait = ifelse(trait == 'LobarWM_Z MTR_fam',"MTR.WM_fam",trait))
plot.df.volcano_main = plot.df.volcano_main %>% 
  mutate(trait = ifelse(trait == 'LobarGM_Z_fam',"MTR.GM_fam",trait))

table(plot.df.volcano_main %>% select(generation, trait))
write_tsv(plot.df.volcano_main,'results/plot.df.volcano_main_2025-12-02.tsv')

plot.df.volcano_main = plot.df.volcano_main %>% mutate(fdr.p = p.adjust(pvalue,method="BH"))
hist(plot.df.volcano_main %>% pull(pvalue),breaks=100); 
pvalue.max = max(plot.df.volcano_main %>% filter(fdr.p<0.05) %>% pull(pvalue))#0.01916851#0.01360136
abline(v=pvalue.max)#pvalue.max
dim(plot.df.volcano_main)[1]/6/1#831
plot.df.volcano_main = plot.df.volcano_main %>% mutate(fdr.signif = ifelse(fdr.p<0.05,"p.fdr<0.05","p.fdr>=0.05"))
# volcano plots ----
point.size = 1.75
## volcano plots: Volume ----
xlims = max(abs(plot.df.volcano_main$beta))*1.015*c(-1,1)
ylims = c(0,max(-log10(plot.df.volcano_main$pvalue)))

volcano_WM.vol = plot.df.volcano_main %>%
  filter(trait == 'lobar.vol.WM.adjICV_fam', generation=="ados") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             # color=fdr.signif)) + 
             color=`lipid class`)) + 
  geom_point(alpha=0.75,size=point.size,stroke = 0) +
  scale_shape_manual(values=c(19,21))+ 
  # scale_color_manual(values=c("#073B4C",'grey'))+
  scale_color_manual(
    values=plot.df.volcano_main_GM.vol_adults$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_GM.vol_adults$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_GM.vol_adults$cols,alpha=0.5))
  ) + 
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.vol
##
volcano_WM.vol_adults = plot.df.volcano_main %>%
  filter(trait == 'lobar.vol.WM.adjICV_fam', generation=="adults") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.vol_adults
## 
volcano_GM.vol = plot.df.volcano_main %>%
  filter(trait == 'lobar.vol.GM.adjICV_fam', generation=="ados") %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_WM.vol$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_WM.vol$class2.labels,
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
##
volcano_GM.vol_adults = plot.df.volcano_main %>%
  filter(trait == 'lobar.vol.GM.adjICV_fam', generation=="adults") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_GM.vol
volcano_GM.vol_adults

## volcano plots: T1wSI ----
volcano_WM.T1wSI = plot.df.volcano_main %>%
  filter(trait == 'NormWM_fam', generation=="ados") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.T1wSI
volcano_WM.T1wSI_adults = plot.df.volcano_main %>%
  filter(trait == 'NormWM_fam', generation=="adults") %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.T1wSI_adults

volcano_GM.T1wSI = plot.df.volcano_main %>%
  filter(trait == 'NormGM_fam', generation=="ados") %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_GM.T1wSI
volcano_GM.T1wSI_adults = plot.df.volcano_main %>%
  filter(trait == 'NormGM_fam', generation=="adults") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_GM.T1wSI_adults

## volcano plots: MTR ----
volcano_WM.MTR = plot.df.volcano_main %>%
  filter(trait == 'MTR.WM_fam', generation=="ados") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.MTR
volcano_WM.MTR_adults = plot.df.volcano_main %>%
  filter(trait == 'MTR.WM_fam', generation=="adults") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_WM.MTR_adults

volcano_GM.MTR = plot.df.volcano_main %>%
  filter(trait == 'MTR.GM_fam', generation=="ados") %>%
  ggplot(aes(x=beta,y=-log10(pvalue),
             shape=fdr.signif,
             fill=fdr.signif,
             color=fdr.signif)) + 
  geom_point(alpha=0.75,size=point.size) +
  scale_shape_manual(values=c(21,19))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  scale_fill_manual(values=c("#073B4C",'grey'))+ 
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_GM.MTR
volcano_GM.MTR_adults = plot.df.volcano_main %>%
  filter(trait == 'MTR.WM_fam', generation=="adults") %>%
  # ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
  ggplot(aes(x=beta,y=-log10(pvalue))) + 
  geom_point(shape=19,alpha=0.25,size=point.size,color="#073B4C") + 
  scale_color_manual(
    values=plot.df.volcano_main_WM.vol$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_WM.vol$class2.labels,
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
  geom_hline(yintercept = -log10(signif.level), linetype=2, color="darkred") + 
  geom_vline(xintercept = 0, color="darkgrey") #+ ggtitle('volcano_WM.vol') 
volcano_GM.MTR_adults

# for manuscript ----
## generate association results ----
df_wide <- plot.df.volcano_main %>% 
  select(metabolite_id,`lipid class`,trait,generation,
         beta,se,pvalue,fdr.p,N) %>%
  pivot_wider(
    names_from = c(trait,generation),
    values_from = c(beta,se,pvalue,fdr.p,N)
  )

## generate Supplementary table source ----
df_wide <- plot.df.volcano_main %>% 
  select(metabolite_id,`lipid class`,trait,generation,
         beta,se,pvalue,fdr.p) %>%
  pivot_wider(
    names_from = c(trait,generation),
    values_from = c(beta,se,pvalue,fdr.p)
  )

ordered_cols = c("metabolite_id", "lipid class");for (i in unique(plot.df.volcano_main$trait)){
  ordered_cols_i = paste(c('beta','se','pvalue','fdr.p'),i,sep='_')
  ordered_cols_i = paste(rep(ordered_cols_i,2),rep(c("ados","adults"),each=4),sep="_")
  ordered_cols = c(ordered_cols,ordered_cols_i)
  }

write_tsv(
  subset(df_wide,select = ordered_cols),
  file.path(wd,"results","TabS1_assoc_results_2025-11-27.tsv")
)
cat(ordered_cols,sep="\n",file = "results/TabS1_colnames.txt")

## generate plot ----
### volume ----
volcano_volume = (volcano_WM.vol+ggtitle('WM_ados')) + 
  (volcano_GM.vol+ggtitle('GM_ados')) + 
  (volcano_WM.vol_adults+ggtitle('WM_adults'))+ 
  (volcano_GM.vol_adults+ggtitle('GM_adults'))
volcano_volume = (volcano_WM.vol) + 
  (volcano_WM.vol_adults)

volcano_volume  + plot_layout(nrow = 2, axes= 'collect') & 
  coord_cartesian(xlim=xlims,ylim = ylims) &
  theme(legend.position = 'none') &
  labs(x = "Estimate (in SD units)", y = "-log10(p-value)") &
  geom_hline(yintercept = 0)
ggsave('results/[manuscript]_Fig1_volcalnoplots_volume_2025-12-01.png',
       width=10.99*1.45*0.75,height=14.99*1.45,units='cm',dpi=300)

### T1wSI ----
volcano_T1wSI = (volcano_WM.T1wSI+ggtitle('WM_ados')) + (volcano_GM.T1wSI+ggtitle('GM_ados')) + 
  (volcano_WM.T1wSI_adults+ggtitle('WM_adults'))+ (volcano_GM.T1wSI_adults+ggtitle('GM_adults'))
volcano_T1wSI
volcano_T1wSI = (volcano_WM.T1wSI) + (volcano_GM.T1wSI) + 
  (volcano_WM.T1wSI_adults)+ (volcano_GM.T1wSI_adults)

volcano_T1wSI  + plot_layout(nrow = 2, axes= 'collect') & 
  coord_cartesian(xlim=xlims,ylim = ylims) &
  labs(x = "Estimate (in SD units)", y = "-log10(p-value)") &
  geom_hline(yintercept = 0)
# ggsave('results/[manuscript]_Fig1_volcalnoplots_T1wSI_2025-11-27.png',
#        width=10.99*1.45,height=14.99*1.45,units='cm',dpi=300)

### MTR ----
volcano_MTR_wi_title = (volcano_WM.MTR+ggtitle('WM_ados')) + (volcano_GM.MTR+ggtitle('GM_ados')) + 
  (volcano_WM.MTR_adults+ggtitle('WM_adults'))+ (volcano_GM.MTR_adults+ggtitle('GM_adults'))
volcano_MTR_wi_title

volcano_MTR = (volcano_WM.MTR) + (volcano_GM.MTR) + 
  (volcano_WM.MTR_adults)+ (volcano_GM.MTR_adults)

volcano_MTR  + plot_layout(nrow = 2, axes= 'collect') & 
  coord_cartesian(xlim=xlims,ylim = ylims) &
  labs(x = "Estimate (in SD units)", y = "-log10(p-value)") &
  geom_hline(yintercept = 0)
# ggsave('results/[manuscript]_Fig1_volcalnoplots_MTR_2025-11-27.png',
#        width=10.99*1.45,height=14.99*1.45,units='cm',dpi=300)

# manuscript result section ----
p_WM = volcano_WM.T1wSI + volcano_WM.T1wSI_adults + 
  volcano_WM.MTR + volcano_WM.MTR_adults + 
  volcano_WM.vol + volcano_WM.vol_adults 
  
p_WM + plot_layout(nrow=2, byrow=F, guides = 'collect', axes = 'collect')& 
  coord_cartesian(xlim=xlims,ylim = ylims) &
  labs(x = "Estimate (in SD units)", y = "-log10(p-value)") &
  geom_hline(yintercept = 0)
ggsave('results/[manuscript]_Fig1_volcalnoplots_MTR_2025-12-01.png',
       width=10.99*2*1.45,height=14.99*1.45,units='cm',dpi=300)
range(plot.df.volcano_main %>% filter(generation=="ados") %>% pull(N))
range(plot.df.volcano_main %>% filter(generation=="adults") %>% pull(N))
## ados' results ----
### WM phenotypes ----
## (-ve) with volume, and (+ve) with T1wSI and MTR
df_wide %>% filter(beta_lobar.vol.WM.adjICV_fam_ados<0, fdr.p_lobar.vol.WM.adjICV_fam_ados<0.05) %>%dim#525
df_wide %>% #filter(beta_lobar.vol.WM.adjICV_fam_ados<0, fdr.p_lobar.vol.WM.adjICV_fam_ados<0.05) %>% 
  filter(beta_NormWM_fam_ados>0, fdr.p_NormWM_fam_ados<0.05) %>% dim#518/618
df_wide %>% 
  # filter(beta_lobar.vol.WM.adjICV_fam_ados<0, fdr.p_lobar.vol.WM.adjICV_fam_ados<0.05) %>% 
  filter(beta_MTR.WM_fam_ados>0, fdr.p_MTR.WM_fam_ados<0.05) %>% dim#517/632

## (+ve) with volume, and (-ve) with T1wSI and MTR
df_wide %>% 
  filter(beta_lobar.vol.WM.adjICV_fam_ados>0, fdr.p_lobar.vol.WM.adjICV_fam_ados<0.05) %>%dim#10
df_wide %>% 
  # filter(beta_lobar.vol.WM.adjICV_fam_ados>0, fdr.p_lobar.vol.WM.adjICV_fam_ados<0.05) %>% 
  filter(beta_NormWM_fam_ados<0, fdr.p_NormWM_fam_ados<0.05) %>% dim#10/18
df_wide %>% 
  # filter(beta_lobar.vol.WM.adjICV_fam_ados>0, fdr.p_lobar.vol.WM.adjICV_fam_ados<0.05) %>% 
  filter(beta_MTR.WM_fam_ados<0, fdr.p_MTR.WM_fam_ados<0.05) %>% dim#4
df_wide %>% 
  mutate(p_WM_ados = pvalue_lobar.vol.WM.adjICV_fam_ados*pvalue_NormWM_fam_ados*pvalue_MTR.WM_fam_ados) %>% 
  arrange(p_WM_ados) %>%
  select(metabolite_id,
         pvalue_lobar.vol.WM.adjICV_fam_ados,
         pvalue_NormWM_fam_ados,
         pvalue_MTR.WM_fam_ados)
# A tibble: 831 × 4
#   metabolite_id pvalue_lobar.vol.WM.adjICV_fam_ados pvalue_NormWM_fam_ados pvalue_MTR.WM_fam_ados
#   <chr>                                       <dbl>                  <dbl>                  <dbl>
# 1 TAG54:7FA16:1                       0.00000000138               2.89e-13               8.50e-13
# 2 TAG53:0FA16:0                       0.0000000531                6.80e-15               1.13e-10
# 3 TAG54:6FA16:0                       0.0000000352                3.59e-14               1.41e-10
# 4 TAG54:4FA20:4                       0.00000000494               1.74e-14               4.10e- 9
# 5 TAG54:6FA22:6                       0.00000000437               1.39e-14               1.08e- 8
# 6 TAG54:5FA20:4                       0.0000000356                8.30e-14               2.52e-10
# 7 TAG54:7FA22:6                       0.0000000131                4.25e-13               2.43e-10
# 8 TAG54:6FA16:1                       0.0000000153                3.68e-11               9.89e-12
# 9 TAG52:6FA20:5                       0.000000751                 7.30e-13               3.36e-11
#10 TAG58:8FA20:4                       0.00000000229               2.25e-12               3.73e- 9

### GM phenotypes ----
df_wide %>% filter(fdr.p_lobar.vol.GM.adjICV_fam_ados<0.05) %>%dim#11
df_wide %>% 
  filter(fdr.p_lobar.vol.GM.adjICV_fam_ados<0.05) %>% 
  arrange(pvalue_lobar.vol.GM.adjICV_fam_ados) %>%
  select(metabolite_id,
         beta_lobar.vol.GM.adjICV_fam_ados,
         se_lobar.vol.GM.adjICV_fam_ados,
         pvalue_lobar.vol.GM.adjICV_fam_ados,
         fdr.p_lobar.vol.GM.adjICV_fam_ados)
#    metabolite_id `lipid class`      beta_lobar.vol.GM.ad…¹ pvalue_lobar.vol.GM.…² fdr.p_lobar.vol.GM.a…³
# 1 PC(18:0/18:0) Phospholipids                      0.106                 0.00260                 0.0108
df_wide %>% 
  filter(fdr.p_NormGM_fam_ados<0.05) %>%dim#187
df_wide %>% 
  filter(fdr.p_NormGM_fam_ados<0.05) %>% 
  arrange(pvalue_NormGM_fam_ados) %>%
  select(metabolite_id,
         beta_NormGM_fam_ados,
         se_NormGM_fam_ados,
         pvalue_NormGM_fam_ados,
         fdr.p_NormGM_fam_ados)

df_wide %>% filter(fdr.p_MTR.GM_fam_ados<0.05) %>%dim#643
#    metabolite_id `lipid class`      beta_NormGM_fam_ados pvalue_NormGM_fam_ados fdr.p_NormGM_fam_ados
  # 1 MAG(16:1)     Acylglycerols                     0.118              0.0000786              0.000466
df_wide %>% 
  filter(fdr.p_MTR.GM_fam_ados<0.05) %>%dim#187
df_wide %>% 
  filter(fdr.p_MTR.GM_fam_ados<0.05) %>% 
  arrange(pvalue_MTR.GM_fam_ados) %>%
  select(metabolite_id,
         beta_MTR.GM_fam_ados,
         se_MTR.GM_fam_ados,
         pvalue_MTR.GM_fam_ados,
         fdr.p_MTR.GM_fam_ados)
#   metabolite_id  `lipid class` beta_MTR.GM_fam_ados pvalue_MTR.GM_fam_ados fdr.p_MTR.GM_fam_ados
# 1 TAG54:7FA16:1  Acylglycerols                0.240               6.10e-13              4.12e-10

## adults' results ----
### WM phenotypes ----
df_wide %>% filter(fdr.p_lobar.vol.WM.adjICV_fam_adults<0.05) %>%dim#6
df_wide %>% 
  filter(fdr.p_lobar.vol.WM.adjICV_fam_adults<0.05) %>% 
  arrange(pvalue_lobar.vol.WM.adjICV_fam_adults) %>%
  select(metabolite_id,
         beta_lobar.vol.WM.adjICV_fam_adults,
         se_lobar.vol.WM.adjICV_fam_adults,
         pvalue_lobar.vol.WM.adjICV_fam_adults,
         fdr.p_lobar.vol.WM.adjICV_fam_adults)
#   metabolite_id  `lipid class`      beta_lobar.vol.WM.ad…¹ pvalue_lobar.vol.WM.…² fdr.p_lobar.vol.WM.a…³
# 1 Glucose        Inflammation/AA/m…                 -0.214              0.0000416               0.000264
df_wide %>% filter(fdr.p_NormWM_fam_adults<0.05) %>%dim#23
data.table(df_wide) %>% 
  filter(fdr.p_NormWM_fam_adults<0.05) %>% 
  arrange(pvalue_NormWM_fam_adults) %>%
  select(metabolite_id,
         beta_NormWM_fam_adults,
         se_NormWM_fam_adults,
         pvalue_NormWM_fam_adults,
         fdr.p_NormWM_fam_adults) %>% slice(1)
#    metabolite_id  `lipid class`     beta_NormWM_fam_adults pvalue_NormWM_fam_ad…¹ fdr.p_NormWM_fam_adu…²
# 1 MAG(16:1)      Acylglycerols                      0.168              0.0000215               0.000145
df_wide %>% filter(fdr.p_MTR.WM_fam_adults<0.05) %>%dim#17
df_wide %>% filter(fdr.p_MTR.WM_fam_adults<0.05) %>%
  arrange(pvalue_MTR.WM_fam_adults) %>%
  select(metabolite_id,
         beta_MTR.WM_fam_adults,
         se_MTR.WM_fam_adults,
         pvalue_MTR.WM_fam_adults,
         fdr.p_MTR.WM_fam_adults) %>% slice(1)
#   metabolite_id `lipid class`      beta_MTR.WM_fam_adults pvalue_MTR.WM_fam_ad…¹ fdr.p_MTR.WM_fam_adu…²
# 1 PC(16:0/16:1) Phospholipids                       0.173               0.000758                0.00356

### GM phenotypes ----
df_wide %>% filter(fdr.p_lobar.vol.GM.adjICV_fam_adults<0.05) %>%dim#9
df_wide %>% 
  filter(fdr.p_lobar.vol.GM.adjICV_fam_adults<0.05) %>% 
  arrange(pvalue_lobar.vol.GM.adjICV_fam_adults) %>%
  select(metabolite_id,
         beta_lobar.vol.GM.adjICV_fam_adults,
         se_lobar.vol.GM.adjICV_fam_adults,
         pvalue_lobar.vol.GM.adjICV_fam_adults,
         fdr.p_lobar.vol.GM.adjICV_fam_adults)
#   metabolite_id  `lipid class`      beta_lobar.vol.GM.ad…¹ pvalue_lobar.vol.GM.…² fdr.p_lobar.vol.GM.a…³
# 1 Glucose        Inflammation/AA/m…                 -0.234              0.0000106              0.0000840
df_wide %>% filter(fdr.p_NormGM_fam_adults<0.05) %>%dim#23
df_wide %>% 
  filter(fdr.p_NormGM_fam_adults<0.05) %>% 
  arrange(pvalue_NormGM_fam_adults) %>%
  select(metabolite_id,
         beta_NormGM_fam_adults,
         se_NormGM_fam_adults,
         pvalue_NormGM_fam_adults,
         fdr.p_NormGM_fam_adults)
#   metabolite_id `lipid class` beta_NormGM_fam_adults pvalue_NormGM_fam_adults fdr.p_NormGM_fam_adults
# 1 TAG51:4FA15:0 Acylglycerols                 -0.115                  0.00487                  0.0205
df_wide %>% filter(fdr.p_MTR.GM_fam_adults<0.05) %>%dim#16
df_wide %>% filter(fdr.p_MTR.GM_fam_adults<0.05) %>%
  arrange(pvalue_MTR.GM_fam_adults) %>%
  select(metabolite_id,
         beta_MTR.GM_fam_adults,
         se_MTR.GM_fam_adults,
         pvalue_MTR.GM_fam_adults,
         fdr.p_NormGM_fam_adults)
#   metabolite_id  `lipid class`     beta_MTR.GM_fam_adults pvalue_MTR.GM_fam_ad…¹ fdr.p_MTR.GM_fam_adu…²
# 1 TAG44:0FA18:0  Acylglycerols                      0.191              0.0000705               0.000462

# checking ----
dim(fread('data/metabodata_adults.txt'))#650
dim(fread('data/metabodata_ados.txt'))

plot.df.volcano_main %>% filter(trait == 'MTR.WM_fam', generation=="ados") %>% 
  select(metabolite_id,beta) %>% rename(beta_MTR_WM = beta) %>%
  left_join(plot.df.volcano_main %>% filter(trait == 'MTR.GM_fam', generation=="ados") %>% 
              select(metabolite_id,beta) %>%  rename(beta_MTR_GM = beta), join_by(metabolite_id)) %>%
  ggplot(aes(x=beta_MTR_WM,y=beta_MTR_GM)) + geom_point() + ggpubr::stat_cor() + geom_abline(intercept = 0, slope=1)

plot.df.volcano_main %>% filter(trait == 'MTR.WM_fam', generation=="adults") %>% 
  select(metabolite_id,beta) %>% rename(beta_MTR_WM = beta) %>%
  left_join(plot.df.volcano_main %>% filter(trait == 'MTR.GM_fam', generation=="adults") %>% 
              select(metabolite_id,beta) %>% rename(beta_MTR_GM = beta), join_by(metabolite_id)) %>%
  ggplot(aes(x=beta_MTR_WM,y=beta_MTR_GM)) + geom_point() + ggpubr::stat_cor()#500*385

d_MTR_ADO1 = fread("~/Downloads/jean_shin_dataanalysis_completesyssubjects_rev1_20251126_231909.txt")
names(d_MTR_ADO1) = str_remove(names(d_MTR_ADO1),'ADO1-MRI-MTRPipeAug2012.')
d_MTR_ADO1 = data.frame(d_MTR_ADO1)
d_MTR_ADO1[,-c(1:3)] = apply(d_MTR_ADO1[,-c(1:3)],2,as.numeric)
d_MTR_ADO1 = d_MTR_ADO1 %>% mutate(across(-c(1:3), ~ ifelse(.x == -666, NA, .x)))
d_MTR_ADO1 %>% describe #n=784

d_MTR_ADO1 %>% ggplot(aes(x=LobarWM_Z,y=LobarGM_Z)) + geom_point() + 
  # ggpubr::stat_cor() + 
  geom_abline(intercept = 0, slope=1)
p_MTR = volcano_WM.MTR  + volcano_GM.MTR #950*500
p_MTR_adults = volcano_WM.MTR_adults  + volcano_GM.MTR_adults #950*500

(volcano_WM.MTR  + volcano_GM.MTR +volcano_WM.MTR_adults  + volcano_GM.MTR_adults ) + plot_layout(nrow=2)

plot.df.volcano_main_WM.vol_adults$volcano_main_all %>% select(metabolite_id,beta) %>%
  left_join(
    plot.df.volcano_main_GM.vol_adults$volcano_main_all %>% select(metabolite_id,beta),
    join_by(metabolite_id)
  ) %>%
  ggplot(aes(x=abs(beta.x),y=abs(beta.y))) + geom_point() + xlab('WM_volume_adults') + ylab('GM_volume_adults') + 
  geom_abline(slope = 1,intercept = 0)

# generate volcano plots (with class-indicating colors) ----
volcano_WM.vol = generate_volcano_plot(plot.df.volcano_main,df=plot.df.volcano_main_GM.vol_adults,'lobar.vol.WM.adjICV_fam','ados')
volcano_WM.vol_adults= generate_volcano_plot(plot.df.volcano_main,df=plot.df.volcano_main_GM.vol_adults,'lobar.vol.WM.adjICV_fam','adults')
(volcano_WM.vol+volcano_WM.vol_adults) + plot_layout(nrow=2,guides = 'collect',axes='collect') & coord_cartesian(xlim = xlims)

volcano_WM.MTR = generate_volcano_plot(plot.df.volcano_main,df=plot.df.volcano_main_GM.vol_adults,'MTR.WM_fam','ados')
volcano_WM.MTR_adults= generate_volcano_plot(plot.df.volcano_main,df=plot.df.volcano_main_GM.vol_adults,'MTR.WM_fam','adults')
(volcano_WM.MTR+volcano_WM.MTR_adults) + plot_layout(nrow=2,guides = 'collect',axes='collect') & coord_cartesian(xlim = xlims)
