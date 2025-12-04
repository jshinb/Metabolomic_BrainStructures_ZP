metaboinfo_adult = fread('data/metabodata_info_adults.txt')
metaboinfo_ado = fread('data/metabodata_info_ados.txt')
names(metaboinfo_ado)[-1] = paste(names(metaboinfo_ado)[-1],'ados',sep="_")
names(metaboinfo_adult)[-1] = paste(names(metaboinfo_adult)[-1],'adults',sep="_")

plot.df.volcano_main = plot.df.volcano_main
  
  
df = plot.df.volcano_main %>% right_join(
  metaboinfo_ado %>% 
    mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")), 
  join_by(metabolite_id)) %>%right_join(
    metaboinfo_adult %>% 
      mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")), 
    join_by(metabolite_id))
df %>% filter(generation == "adults",trait=="MTR.WM_fam") %>% arrange(N) %>% slice(1:10) %>%
  select(metabolite_id,N,trait,n_ados,missing_ados,n_adults,
         feature_missingness_ados,
         feature_missingness_adults 
  )

df %>% filter(generation == "ados",trait=="MTR.WM_fam") %>% arrange(N) %>% slice(1:10) %>%
  select(metabolite_id,N,trait,n_ados,missing_ados,n_adults,
         feature_missingness_ados,
         feature_missingness_adults 
  )
df_wide <- plot.df.volcano_main %>%
  select(metabolite_id,`lipid class`,trait,generation,
         beta,se,pvalue,fdr.p) %>%
  pivot_wider(
    names_from = c(trait,generation),
    values_from = c(beta,se,pvalue,fdr.p)
  )

df_wide %>% ggplot(aes(y=(beta_lobar.vol.WM.adjICV_fam_ados), 
                       x=(beta_lobar.vol.WM.adjICV_fam_adults),
                       color=`lipid class`,fill=`lipid class`)) + 
  geom_point(shape=19,alpha=0.75,size=point.size) + 
  scale_color_manual(
    values=plot.df.volcano_main_WM.vol$cols,
    name="Lipid classes",
    labels = plot.df.volcano_main_WM.vol$class2.labels,
    guide = guide_legend(
      override.aes = list(color=plot.df.volcano_main_GM.vol$cols,alpha=0.5))
  ) + 
  scale_fill_manual(
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
  geom_abline(slope = 1,intercept = 0,color="darkgrey") + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_vline(xintercept = 0, color="darkgrey") +
  coord_cartesian(xlim = xlims, ylim=xlims)

### FDR-correction-significant levels

