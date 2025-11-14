rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)


#9. volcano-plots -------------------------------------------------------------
rm(volcano_main,volcano_main_all,p,top_main)

load("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_adjICV_2025-03-26.Rdata")
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
# generating volcano plots ----
point.size = 0.9

## volume ----
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

head(plot.df.volcano_main_lobarVolWM$volcano_main_all)
dim(plot.df.volcano_main_lobarVolWM$volcano_main_all)

# [plot]_scatter plots: assoc_beta ----
point.size=0.95
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

top_main <- volcano_main_3brain %>% 
  filter(beta.Vol<0, beta.nT1wSI>0, beta.MTR>0) %>% slice(1:10) %>% 
  bind_rows(
    volcano_main_3brain %>% filter(beta.Vol>0, beta.nT1wSI<0, beta.MTR<0) %>% slice(1:10)
  )

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

## generate scatter plots and create a png file ----
scatter_lobarVolWM_T1wSI + scatter_lobarVolWM_MTR + scatter_T1wSI_MTR +
  plot_layout(guides='collect') &
  theme(legend.position = 'none') &
  coord_cartesian(ylim=c(-0.225,0.35),xlim=c(-0.225,0.225)) &
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=9),
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) #&
# geom_abline(slope=1,intercept = 0, color='grey')
ggsave("results/scatterplot_sex_combined_unadjBMI.png",
       width=11.5,height=3.5,units='in',dpi=300)

cor(volcano_main_3brain %>%
      select(beta.Vol,beta.nT1wSI,beta.MTR))

## correlation plot ----
