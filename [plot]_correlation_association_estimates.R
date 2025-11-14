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
load("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_adjICV_2025-03-26.Rdata")
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

## correlation plot ----
ggcorr_all = create_ggcorrplot (
  data = volcano_main_3brain ,
  sel_colnames = c("beta.Vol","beta.nT1wSI","beta.MTR"),
  new_colnames = c("beta-WM-Vol","beta-WM-SI","beta-WM-MTR"),
  ID_colname = "metabolite_id"
)

ggcorr_all
