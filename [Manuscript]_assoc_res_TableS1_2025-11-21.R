rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

# ados assoc results ----
volcano_main_lobarVolWM_ados = fread("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.tsv")
volcano_main_lobarVolGM_ados = fread("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_ados_unadjBMI_2025-05-26.tsv")
volcano_main_lobarVolWM_ados_wi_classes = add_class2(volcano_main_lobarVolWM_ados)
head(volcano_main_lobarVolWM_ados_wi_classes$volcano_main_all)

volcano_main_NormWM_ados = fread("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.tsv")
volcano_main_NormGM_ados = fread("results/volcanoplot_all_NormGM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.tsv")

volcano_main_MTR_WM_ados = fread("results/volcanoplot_all_LobarWM_Z MTR_sex-combined_ados_unadjBMI_2025-04-08.tsv")
# volcano_main_MTR_GM_ados #*[GM-MTR was not analysed]

assoc_res_ados = 
  volcano_main_lobarVolWM_ados_wi_classes$volcano_main_all %>% select(metabolite_id,`lipid class`,beta,se,pvalue) %>%
  rename("EstimateWM_Volume"=beta,SEWM_Volume=se,PWM_Volume=pvalue) %>%
  left_join(
    volcano_main_lobarVolGM_ados %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(EstimateGM_Volume=beta,SEGM_Volume=se,PGM_Volume=pvalue)
  )%>%
  left_join(
    volcano_main_NormWM_ados %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_WM_SI=beta,SE_WM_SI=se,P_WM_SI=pvalue)
  )%>%
  left_join(
    volcano_main_NormGM_ados %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_GM_SI=beta,SE_GM_SI=se,P_GM_SI=pvalue)
  )%>%
  left_join(
    volcano_main_MTR_WM_ados %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_WM_MTR=beta,SE_MTR_SI=se,P_MTR_SI=pvalue)
  )

write_tsv(assoc_res_ados,file.path("outputs","TableS1_assoc_ados_2025-11-21.txt"))

# adults assoc results ----
volcano_main_lobarVolWM_adults = fread("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.tsv")
volcano_main_lobarVolGM_adults = fread("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.tsv")

## need to add Nightingale metabolite (n=7)
volcano_main_NormWM_adults = fread("results/volcanoplot_NormWM_sex-combined_adults_unadjBMI_2025-01-10.tsv")
volcano_main_NormGM_adults = fread("results/volcanoplot_NormGM_sex-combined_adults_unadjBMI_2025-01-10.tsv")

# volcano_main_MTR_WM_adults #*[MTR was not analysed]
# volcano_main_MTR_GM_adults 

volcano_main_lobarVolWM_adults_wi_class = add_class2(volcano_main_lobarVolWM_adults,generation="adults")
dim(volcano_main_lobarVolWM_adults_wi_class$volcano_main_all)
table(volcano_main_lobarVolWM_adults_wi_class$volcano_main_all$class,useNA='a')

assoc_res_adults = volcano_main_lobarVolWM_adults_wi_class$volcano_main_all %>%
  select(metabolite_id,`lipid class`,beta,se,pvalue) %>%
  left_join(
    volcano_main_lobarVolGM_adults %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(EstimateGM_Volume=beta,SEGM_Volume=se,PGM_Volume=pvalue)
  )%>%left_join(
    volcano_main_NormWM_adults %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_WM_SI=beta,SE_WM_SI=se,P_WM_SI=pvalue)
  )%>%
  left_join(
    volcano_main_NormGM_adults %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_GM_SI=beta,SE_GM_SI=se,P_GM_SI=pvalue)
  )
head(assoc_res_adults)
write_tsv(assoc_res_adults,file.path("outputs","TableS1_assoc_adults_2025-11-21.txt"))
