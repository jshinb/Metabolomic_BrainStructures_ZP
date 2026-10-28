# CREATED: 2026-10-27 
# Modified: 2026-10-28
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
source(file.path(wd,'99_scripts','[function]_define_fit_function_2025-01-09.R'))
select = dplyr::select
describe = psych::describe

# 1. results data frame ----
metabodata_info = fread("data/metabodata_info_ados.txt")
metabodata_info %>% head

load("results/sample_size_matched_summary.Rdata",verbose = T)
traits = c("Volume.WM.adjAgeSexICV", "T1wSI.WM.adjAgeSex", "MTR.WM.adjAgeSex"
)
# 2. add class info to the result data frame ----
volcano_main_all = 
  matched_summary %>% rename(metabolite_id = metabolite) %>%
  left_join(
    metabodata_info %>%
      mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>% 
      dplyr::select(metabolite_id,class)
    ) %>%
  mutate(trait = factor(trait, levels = traits),
         generation = 'ado1_resample') 

class2.levels = c("Cholesteryl esters",
                  "Acylglycerols",
                  "Phospholipids",
                  "Lysophospholipids",
                  "Sphingolipids",
                  "Inflammation/AA/metabolism")
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

tab = table(volcano_main_all$`lipid class`)/length(unique(volcano_main_all$trait))
class2.labels = paste(names(table(volcano_main_all$`lipid class`))," (n=",
                      tab,")",sep='')
# 3. volcano plot ----
point.size = 1;plot.cols = cols
generation_plot = 'ado1_resample'
#signif.level = 0.05 #*[FOR NOW -- NEED TO THINK ABOUT THIS] 
load('results/signif_level_used_for_Fig1_apr2026.Rdata',verbose = T)#0.01889204
p = list();for (ti in 1:3) {
  trait_plot = traits[ti];

p[[ti]] = volcano_main_all %>% 
  mutate(prop.signif = ifelse(prop_sig>0.9, "%p<0.05 >90","%p<0.05 <=90")) %>%
  mutate(prop.signif = factor(prop.signif,levels=c("%p<0.05 >90","%p<0.05 <=90"))) %>%
  filter(trait == trait_plot, generation==generation_plot) %>%
  ggplot(aes(x=median_beta,y=-log10(mean_p),
             shape=prop.signif,
             fill=prop.signif,
             color=prop.signif)) + 
              # color=prop_sig)) + 
             #color=`lipid class`)) + 
  geom_point(alpha=0.5,size=point.size,stroke = 1) +
  scale_shape_manual(values=c(19,21))+ 
  scale_color_manual(values=c("#073B4C",'grey'))+
  # scale_color_manual(
  #   values=plot.cols,
  #   name="Lipid classes",
  #   labels = class2.labels,
  #   guide = guide_legend(
  #     override.aes = list(color=plot.cols,alpha=0.5))
  # ) + 
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
  geom_vline(xintercept = 0, color="darkgrey") + ggtitle(trait_plot) 
}
load("results/separate_Fig1_volcanoplots.Rdata",verbose = T)

p_resample_ados_adults = p[[1]] + p[[2]] + p[[3]] + 
  volcano_WM.vol_adults + volcano_WM.T1wSI_adults + volcano_WM.MTR_adults + 
  plot_layout(ncol=3, axes='collect', guides = 'collect') &
  coord_cartesian(ylim = c(0,13),xlim=c(-0.3,0.3));p_resample_ados_adults

ggsave(paste0('results/[manuscript]_sample_size_matched_volcalnoplots_ados_adults_',Sys.Date(),'.png'),
       plot=p_resample_ados_adults & ggtitle(NULL) &theme(legend.position = 'none')&xlab(NULL) &ylab(NULL),
       height=14.61,width=19.71,units='cm',dpi=300)
#

group = 'ados_resample';for (ti in 1:3) {
  trait_plot = traits[ti];
  
  plot.df = volcano_main_all %>% filter(trait == trait_plot)
  ind_sig = plot.df$mean_p < signif.level;sum(ind_sig)
  nmetabo <- c(sum(ind_sig),sum(!ind_sig))#830
  mylabels <- c("p<0.0189","p>=0.0189")
  df = data.frame(group = mylabels, value=nmetabo)
  
  # Basic piechart
  p_pi = ggplot(df, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    # geom_text(aes(label = group), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values=c("#9babbe","#dbd9da")) +
    theme(legend.position='none') + 
    theme_void();p_pi # remove background, grid, numeric labels
  pie_fname = paste0('[manuscript]_Fig1_pie_volcano_',trait_plot,'_',group,'_',Sys.Date(),'.png')
  ggsave(filename = pie_fname,
         plot=p_pi,
         width=6,height=5.5,units='in',dpi=300)
}
