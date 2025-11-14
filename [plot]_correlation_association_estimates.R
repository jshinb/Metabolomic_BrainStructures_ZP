rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)


# load lipid-brain association results (sex-combined) ----
load("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
volcano_main_NormWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.Rdata")
volcano_main_lobarVolWM = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_LobarWM_Z MTR_sex-combined_ados_unadjBMI_2025-04-08.Rdata")
volcano_main_MTR_WM = volcano_main
rm(volcano_main)

# add color codings to lipid classes ---
plot.df.volcano_main_lobarVolWM = add_class2(volcano_main_lobarVolWM$data)
head(plot.df.volcano_main_lobarVolWM$volcano_main_all)

plot.df.volcano_main_NormWM = add_class2(volcano_main_NormWM$data)
head(plot.df.volcano_main_NormWM$volcano_main_all)

plot.df.volcano_main_MTR_WM = add_class2(volcano_main_MTR_WM$data)
volcano_main_NormWM$data = plot.df.volcano_main_NormWM$volcano_main_all

# merge all three results
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

# sex-specific ----
select = dplyr::select
load("results/volcanoplot_all_NormWM_F_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
volcano_main_NormWM_F = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_F_ados_unadjBMI_2025-04-04.Rdata")
volcano_main_lobarVolWM_F = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_LobarWM_Z MTR_F_ados_unadjBMI_2025-04-08.Rdata")
volcano_main_MTR_WM_F = volcano_main

load("results/volcanoplot_all_NormWM_M_ados_unadjBMI_unadjICV_2025-03-26.Rdata")
volcano_main_NormWM_M = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_lobar.vol.WM.adjICV_M_ados_unadjBMI_2025-04-04.Rdata")
volcano_main_lobarVolWM_M = volcano_main
rm(volcano_main)

load("results/volcanoplot_all_LobarWM_Z MTR_M_ados_unadjBMI_2025-04-08.Rdata")
volcano_main_MTR_WM_M = volcano_main

beta_3brain_F =  volcano_main_lobarVolWM_F$data %>% select(metabolite_id,beta) %>% rename(beta_Vol.WM=beta) %>%
  left_join(volcano_main_NormWM_F$data %>% select(metabolite_id,beta) %>% rename(beta_Norm.WM=beta), 
            join_by(metabolite_id)) %>%
  left_join(volcano_main_MTR_WM_F$data %>% select(metabolite_id,beta) %>% rename(beta_MTR.WM=beta), 
            join_by(metabolite_id))

beta_3brain_M =  volcano_main_lobarVolWM_M$data %>% select(metabolite_id,beta) %>% rename(beta_Vol.WM=beta) %>%
  left_join(volcano_main_NormWM_M$data %>% select(metabolite_id,beta) %>% rename(beta_Norm.WM=beta), 
            join_by(metabolite_id)) %>%
  left_join(volcano_main_MTR_WM_M$data %>% select(metabolite_id,beta) %>% rename(beta_MTR.WM=beta), 
            join_by(metabolite_id))

## correlation plot
new_colnames_beta = c("WM-Vol","WM-SI","WM-MTR")
ggcorr_all = create_ggcorrplot (
  data = volcano_main_3brain ,
  sel_colnames = c("beta.Vol","beta.nT1wSI","beta.MTR"),
  new_colnames = new_colnames_beta,
  ID_colname = "metabolite_id"
)

## female 
ggcorr_female = create_ggcorrplot (
  data = beta_3brain_F ,
  sel_colnames = c("beta_Vol.WM","beta_Norm.WM","beta_MTR.WM"),
  new_colnames = new_colnames_beta,
  ID_colname = "metabolite_id"
)

## male
ggcorr_male = create_ggcorrplot (
  data = beta_3brain_M ,
  sel_colnames = c("beta_Vol.WM","beta_Norm.WM","beta_MTR.WM"),
  new_colnames = new_colnames_beta,
  ID_colname = "metabolite_id"
)

ggcorr_plot_beta =
  (ggcorr_all$p_ggcorr + theme(axis.text.y = element_text(size = 10)))+ 
  (ggcorr_female$p_ggcorr + theme(axis.text.y = element_blank())) + 
  (ggcorr_male$p_ggcorr + theme(axis.text.y = element_blank())) + 
  plot_layout(guides = 'collect') &
  theme(axis.text.x = element_text(size = 10))

ggsave(file.path(opt$output.dir,paste0("brain_corrplot_beta.jpg")), 
       bg = "white",
       plot = ggcorr_plot_beta, 
       width = 8.5*3*1.25, 
       height = 4.5**1.25, 
       units = "cm", 
       dpi = 300)

# [plot]: volcano plots: this part is not working currently (need to fix) ----
not.work <- function(){
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
}
