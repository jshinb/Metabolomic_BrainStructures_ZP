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

# load lipid-brain association results (sex-combined) ----
load("results/volcanoplot_all_NormGM_sex-combined_adults_unadjBMI_2025-01-09.Rdata")
volcano_main_NormGM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.Rdata")
volcano_main_lobarVolGM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_gm_sex-combined_adults_unadjBMI_2025-11-25.Rdata")
volcano_main_MTR_GM = volcano_main
rm(volcano_main)

# add color codes to lipid classes ---
plot.df.volcano_main_lobarVolGM = add_class2(volcano_main_lobarVolGM$data)
head(plot.df.volcano_main_lobarVolGM$volcano_main_all)

plot.df.volcano_main_NormGM = add_class2(volcano_main_NormGM$data)
head(plot.df.volcano_main_NormGM$volcano_main_all)

plot.df.volcano_main_MTR_GM = add_class2(volcano_main_MTR_GM$data)
volcano_main_NormGM$data = plot.df.volcano_main_NormGM$volcano_main_all

# merge all three results
d_TabS1 = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/results/TabS1_assoc_results_2025-11-27.tsv")
lipids_incl = d_TabS1$metabolite_id
volcano_main_3brain = plot.df.volcano_main_lobarVolGM$volcano_main_all %>% 
  filter(metabolite_id %in% lipids_incl) %>%
  select(metabolite_id,class,class2,`lipid class`,beta,pvalue) %>%
  rename(beta.Vol = beta, pvalue.Vol = pvalue) %>% 
  left_join(
    plot.df.volcano_main_NormGM$volcano_main_all %>% 
      filter(metabolite_id %in% lipids_incl) %>%
      select(metabolite_id,class,class2,`lipid class`,beta,pvalue) %>%
      rename(beta.nT1wSI = beta, pvalue.nT1wSI = pvalue)
  ) %>%
  left_join(
    plot.df.volcano_main_MTR_GM$volcano_main_all %>% 
      filter(metabolite_id %in% lipids_incl) %>%
      select(metabolite_id,class,class2,`lipid class`,beta,pvalue) %>%
      rename(beta.MTR = beta, pvalue.MTR = pvalue)
  ) %>%
  mutate(p.multi = pvalue.Vol*pvalue.nT1wSI*pvalue.MTR) %>%
  arrange(p.multi)
head(volcano_main_3brain)

# correlation plot ----
new_colnames_beta = c("GM-Vol","GM-T1wSI","GM-MTR")
ggcorr_all = create_ggcorrplot (
  data = volcano_main_3brain ,
  sel_colnames = c("beta.Vol","beta.nT1wSI","beta.MTR"),
  new_colnames = new_colnames_beta,
  ID_colname = "metabolite_id"
)

ggsave(file.path(opt$output.dir,paste0("brain_corrplot_beta_GM_adults.jpg")), 
       bg = "white",
       # plot = ggcorr_plot_beta, 
       plot = (ggcorr_all$p_ggcorr + theme(axis.text.y = element_text(size = 10))),
       width = 8.5*1*1.25, 
       height = 4.5*1*1.25, 
       units = "cm", 
       dpi = 300)

# scatter plot ----
point.size=1.25
volcano_main_3brain <- add_class2(volcano_main_3brain)
head( volcano_main_3brain$volcano_main_all)
scatter1 =
  volcano_main_3brain$volcano_main_all %>% 
  ggplot(aes(x=beta.Vol,y=beta.nT1wSI,color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) +
  scale_color_manual(
    values=volcano_main_3brain$cols,
    name="Lipid classes",
    labels = volcano_main_3brain$class2.labels,
    guide = guide_legend(
      override.aes = list(color=volcano_main_3brain$cols,alpha=0.5))
  ) + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey") 

scatter2 = 
  volcano_main_3brain$volcano_main_all %>% 
  ggplot(aes(x=beta.Vol,y=beta.MTR,color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) +
  scale_color_manual(
    values=volcano_main_3brain$cols,
    name="Lipid classes",
    labels = volcano_main_3brain$class2.labels,
    guide = guide_legend(
      override.aes = list(color=volcano_main_3brain$cols,alpha=0.5))
  ) + 
  # geom_smooth(method = 'gam') + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey")

scatter3 = 
  volcano_main_3brain$volcano_main_all %>% 
  ggplot(aes(x=beta.nT1wSI,y=beta.MTR,color=`lipid class`)) + 
  geom_point(shape=19,alpha=0.4,size=point.size) +
  scale_color_manual(
    values=volcano_main_3brain$cols,
    name="Lipid classes",
    labels = volcano_main_3brain$class2.labels,
    guide = guide_legend(
      override.aes = list(color=volcano_main_3brain$cols,alpha=0.5))
  ) + 
  # geom_smooth(method = 'gam') + 
  geom_vline(xintercept = 0, color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey")

plot.lims = range(volcano_main_3brain$volcano_main_all %>% select(beta.Vol,beta.nT1wSI,beta.MTR))
p = (scatter1 + scatter2 + scatter3)+ plot_layout(guides = 'collect') & coord_cartesian(xlim=plot.lims,ylim=plot.lims)
p
ggsave(file.path(opt$output.dir,paste0("brain_scatterplot_beta_GM_adults.jpg")), 
       bg = "white",
       plot = p,
       width = 8.5*3*1.25, 
       height = 4.5*1.5*1.25, 
       units = "cm", 
       dpi = 300)

# [plot]: volcano plots: this part is not working currently (need to fix) ----
not.work <- function(){
  top_main_lobarVolGM = plot.df.volcano_main_lobarVolGM$volcano_main_all %>% 
    filter(beta>0) %>% arrange(pvalue) %>% slice(1:10)
  top_main_lobarVolGM = top_main_lobarVolGM %>% 
    bind_rows(
      plot.df.volcano_main_lobarVolGM$volcano_main_all %>% 
        filter(beta<0) %>% arrange(pvalue) %>% slice(1:10)
    ) 
  
  top_main_NormGM = plot.df.volcano_main_NormGM$volcano_main_all %>% 
    filter(beta>0) %>% arrange(pvalue) %>% slice(1:10)
  top_main_NormGM = top_main_NormGM %>% bind_rows(
    plot.df.volcano_main_NormGM$volcano_main_all %>% 
      filter(beta<0) %>% arrange(pvalue) %>% slice(1:10)
  ) 
  
  top_main_MTR_GM = plot.df.volcano_main_MTR_GM$volcano_main_all %>% 
    filter(beta>0) %>% arrange(pvalue) %>% slice(1:10)
  top_main_MTR_GM = top_main_MTR_GM %>% 
    bind_rows(
      plot.df.volcano_main_MTR_GM$volcano_main_all %>% 
        filter(beta<0) %>% arrange(pvalue) %>% slice(1:10)
    ) 
  
  metabolite_overlap = intersect(
    top_main_lobarVolGM$metabolite_id,
    top_main_NormGM$metabolite_id
  )
  metabolite_overlap = intersect(
    metabolite_overlap,
    top_main_MTR_GM$metabolite_id
  )
  
  top_main_lobarVolGM = top_main_lobarVolGM %>% 
    mutate(label.col = ifelse(metabolite_id %in% metabolite_overlap,"darkred","grey30"))
  top_main_NormGM = top_main_NormGM %>% 
    mutate(label.col = ifelse(metabolite_id %in% metabolite_overlap,"darkred","grey30"))
  top_main_MTR_GM = top_main_MTR_GM %>% 
    mutate(label.col = ifelse(metabolite_id %in% metabolite_overlap,"darkred","grey30"))
  
  point.size = 0.95
  top_main <- volcano_main_3brain %>% 
    filter(beta.Vol<0, beta.nT1wSI>0, beta.MTR>0) %>% slice(1:10) %>% 
    bind_rows(
      volcano_main_3brain %>% filter(beta.Vol>0, beta.nT1wSI<0, beta.MTR<0) %>% slice(1:10)
    )
  
  top_main_lobarVolGM <- plot.df.volcano_main_lobarVolGM$volcano_main_all %>% 
    filter(metabolite_id %in% top_main$metabolite_id) %>%
    mutate(label.col = ifelse(!metabolite_id %in% siglipids,"darkred","grey30"))
  top_main_NormGM <- plot.df.volcano_main_NormGM$volcano_main_all %>% 
    filter(metabolite_id %in% top_main$metabolite_id) %>%
    mutate(label.col = ifelse(!metabolite_id %in% siglipids,"darkred","grey30"))
  top_main_MTR_GM <- plot.df.volcano_main_MTR_GM$volcano_main_all %>% 
    filter(metabolite_id %in% top_main$metabolite_id) %>%
    mutate(label.col = ifelse(!metabolite_id %in% siglipids,"darkred","grey30"))
  
  ### Volume ----
  volcano_lobarVolGM = plot.df.volcano_main_lobarVolGM$volcano_main_all%>%
    ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
    geom_point(shape=19,alpha=0.4,size=point.size) + 
    scale_color_manual(
      values=plot.df.volcano_main_lobarVolGM$cols,
      name="Lipid classes",
      labels = plot.df.volcano_main_lobarVolGM$class2.labels,
      guide = guide_legend(
        override.aes = list(color=plot.df.volcano_main_lobarVolGM$cols,alpha=0.5))
    ) + 
    geom_vline(xintercept = 0, color="darkgrey") + 
    # geom_hline(yintercept = -log10(0.1), color="darkgrey") + 
    geom_text_repel(
      data = top_main_lobarVolGM,
      aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
      position=position_nudge_center(
        #x = 0,#, y = -log10(0.05),
        direction = "split"
      ),
      max.time = 1, max.iter = 1e5, max.overlaps = 25,
      box.padding = 0.3,
      min.segment.length = 0,
      segment.color = scales::alpha(top_main_lobarVolGM$label.col,0.25),
      segment.size = 0.2,
      # segment.color = scales::alpha("lightgray",0.5),
      color = top_main_lobarVolGM$label.col,#'black',
      #fill = scales::alpha('white',0.5),
      size = 2.5) #theme_bw()
  
  ### NormGM ----
  volcano_NormGM = plot.df.volcano_main_NormGM$volcano_main_all%>%
    ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
    geom_point(shape=19,alpha=0.4,size=point.size) + 
    scale_color_manual(
      values=plot.df.volcano_main_NormGM$cols,
      name="Lipid classes",
      labels = plot.df.volcano_main_NormGM$class2.labels,
      guide = guide_legend(
        override.aes = list(color=plot.df.volcano_main_NormGM$cols,alpha=0.5))
    ) +
    geom_vline(xintercept = 0, color="darkgrey") +
    # geom_hline(yintercept = -log10(0.1), color="darkgrey") + 
    geom_text_repel(
      data = top_main_NormGM,
      aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
      position=position_nudge_center(
        #x = 0,#, y = -log10(0.05),
        direction = "split"
      ),
      max.time = 1, max.iter = 1e5, max.overlaps = 25,
      box.padding = 0.3,
      min.segment.length = 0,
      segment.color = scales::alpha(top_main_NormGM$label.col,0.25),
      segment.size = 0.2,
      # segment.color = scales::alpha("lightgray",0.5),
      color = top_main_NormGM$label.col,#'black',
      #fill = scales::alpha('white',0.5),
      size = 2.5) #theme_bw()
  
  ### MTR ----
  volcano_MTR_GM = plot.df.volcano_main_MTR_GM$volcano_main_all %>%
    ggplot(aes(x=beta,y=-log10(pvalue), color=`lipid class`)) + 
    geom_point(shape=19,alpha=0.4,size=point.size) + 
    scale_color_manual(
      values=plot.df.volcano_main_MTR_GM$cols,
      name="Lipid classes",
      labels = plot.df.volcano_main_MTR_GM$class2.labels,
      guide = guide_legend(
        override.aes = list(color=plot.df.volcano_main_MTR_GM$cols,alpha=0.5))
    ) +
    geom_vline(xintercept = 0, color="darkgrey") +
    # geom_hline(yintercept = -log10(0.1), color="darkgrey") + 
    geom_text_repel(
      data = top_main_MTR_GM,
      aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
      position=position_nudge_center(
        #x = 0,#, y = -log10(0.05),
        direction = "split"
      ),
      max.time = 1, max.iter = 1e5, max.overlaps = 25,
      box.padding = 0.3,
      min.segment.length = 0,
      segment.color = scales::alpha(top_main_MTR_GM$label.col,0.25),
      segment.size = 0.2,
      # segment.color = scales::alpha("lightgray",0.5),
      color = top_main_MTR_GM$label.col,#'black',
      #fill = scales::alpha('white',0.5),
      size = 2.5)
}
