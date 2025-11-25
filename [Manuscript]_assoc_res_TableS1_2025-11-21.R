rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)
signif.bonferroni = 0.05/831/6/2
# ados assoc results ----
volcano_main_lobarVolWM_ados = fread("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_ados_unadjBMI_2025-04-04.tsv")
volcano_main_lobarVolGM_ados = fread("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_ados_unadjBMI_2025-05-26.tsv")
volcano_main_lobarVolWM_ados_wi_classes = add_class2(volcano_main_lobarVolWM_ados)
head(volcano_main_lobarVolWM_ados_wi_classes$volcano_main_all)
volcano_main_lobarVolWM_ados %>% mutate(p.fdr = p.adjust(pvalue,method="BH")) %>% filter(beta<0,pvalue<0.05)#586

volcano_main_NormWM_ados = fread("results/volcanoplot_all_NormWM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.tsv")
volcano_main_NormGM_ados = fread("results/volcanoplot_all_NormGM_sex-combined_ados_unadjBMI_unadjICV_2025-03-26.tsv")
lipid1 = volcano_main_NormWM_ados %>% mutate(p.fdr = p.adjust(pvalue,method="BH")) %>% filter(beta>0,p.fdr <0.05) %>% pull(metabolite_id)
volcano_main_MTR_WM_ados = fread("results/volcanoplot_all_LobarWM_Z MTR_sex-combined_ados_unadjBMI_2025-04-08.tsv")
# volcano_main_MTR_GM_ados #*[GM-MTR was not analysed]

assoc_res_ados = 
  volcano_main_lobarVolWM_ados_wi_classes$volcano_main_all %>% select(metabolite_id,`lipid class`,beta,se,pvalue) %>%
  rename(Estimate_WM_Volume=beta,SE_WM_Volume=se,P_WM_Volume=pvalue) %>%
  left_join(
    volcano_main_lobarVolGM_ados %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_GM_Volume=beta,SE_GM_Volume=se,P_GM_Volume=pvalue)
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
      rename(Estimate_WM_MTR=beta,SE_WM_MTR=se,P_WM_MTR=pvalue)
  )

write_tsv(assoc_res_ados,file.path("outputs","TableS1_assoc_ados_2025-11-21.txt"))
assoc_res_ados %>% filter(P_WM_Volume<0.1,P_WM_SI<0.1,P_WM_MTR<0.1)#192
assoc_res_ados = assoc_res_ados %>% mutate(P_WM = P_WM_Volume*P_WM_SI*P_WM_MTR)
assoc_res_ados %>% 
  filter(Estimate_WM_Volume <0, Estimate_WM_SI>0, Estimate_WM_MTR >0) %>% 
  arrange(P_WM) %>% select(1,5,11,17,18) %>% rename(P_WM_overall=P_WM)
assoc_res_ados %>% filter(metabolite_id == "TAG54:4FA20:4") %>% select(P_WM_Volume,P_WM_SI,P_WM_MTR)
assoc_res_ados %>% filter(Estimate_WM_Volume>0,P_WM_Volume<0.05)#19

assoc_res_ados %>% filter(Estimate_WM_Volume<0,P_WM_Volume<0.05, Estimate_WM_SI>0,P_WM_SI<0.05)
assoc_res_ados %>% filter(Estimate_WM_Volume<0,P_WM_Volume<0.05, Estimate_WM_MTR>0,P_WM_MTR<0.05)

assoc_res_ados %>% filter(Estimate_WM_Volume>0,P_WM_Volume<0.05, Estimate_WM_SI<0,P_WM_SI<0.05) %>% arrange(P_WM)
assoc_res_ados %>% filter(Estimate_WM_Volume>0,P_WM_Volume<0.05, Estimate_WM_MTR<0,P_WM_MTR<0.05) %>% arrange(P_WM)

assoc_res_ados %>% filter(Estimate_WM_MTR<0,P_WM_SI<0.05)
assoc_res_ados %>% filter(Estimate_WM_Volume>0) %>% arrange(P_WM) %>% slice(1:3)

# adults assoc results ----
volcano_main_lobarVolWM_adults = fread("results/volcanoplot_all_lobar.vol.WM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.tsv")
volcano_main_lobarVolGM_adults = fread("results/volcanoplot_all_lobar.vol.GM.adjICV_sex-combined_adults_unadjBMI_2025-05-27.tsv")

## need to add Nightingale metabolite (n=7)
volcano_main_NormWM_adults = fread("results/volcanoplot_NormWM_sex-combined_adults_unadjBMI_2025-01-10.tsv")
volcano_main_NormGM_adults = fread("results/volcanoplot_NormGM_sex-combined_adults_unadjBMI_2025-01-10.tsv")
volcano_main_NormWM_adults = volcano_main_NormWM_adults %>%
  mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))

# volcano_main_MTR_WM_adults #*[MTR was not analysed]
# volcano_main_MTR_GM_adults 

volcano_main_lobarVolWM_adults_wi_class = add_class2(volcano_main_lobarVolWM_adults,generation="adults")
dim(volcano_main_lobarVolWM_adults_wi_class$volcano_main_all)
table(volcano_main_lobarVolWM_adults_wi_class$volcano_main_all$class,useNA='a')

assoc_res_adults = volcano_main_lobarVolWM_adults_wi_class$volcano_main_all %>%
  select(metabolite_id,`lipid class`,beta,se,pvalue) %>%
  rename(Estimate_WM_Volume=beta,SE_WM_Volume=se,P_WM_Volume=pvalue) %>%
left_join(
    volcano_main_lobarVolGM_adults %>% select(metabolite_id,beta,se,pvalue) %>%
      rename(Estimate_GM_Volume=beta,SE_GM_Volume=se,P_GM_Volume=pvalue)
  )%>%left_join(
    volcano_main_NormWM_adults %>% select(metabolite_id,beta,se,pvalue) %>%
      mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))%>%
      rename(Estimate_WM_SI=beta,SE_WM_SI=se,P_WM_SI=pvalue)
  )%>%
  left_join(
    volcano_main_NormGM_adults %>% select(metabolite_id,beta,se,pvalue) %>%
      mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))%>%
      rename(Estimate_GM_SI=beta,SE_GM_SI=se,P_GM_SI=pvalue)
  )
headTail(assoc_res_adults)
write_tsv(assoc_res_adults,file.path("outputs","TableS1_assoc_adults_2025-11-21.txt"))

# in adults but not in ados:
# DAG(16:0/16:0)
# DAG(16:0/18:0)
# MAG(14:0)
# MAG(16:0)
# MAG(18:0)
# MAG(20:0)
# MAG(22:5)
# PC(15:0/20:4)
# PC(17:0/20:3)
# PC(18:1/20:2)
# PC(18:1/22:5)
# PC(18:2/18:3)
# PC(18:2/20:5)
# PC(20:0/18:1)
# PC(20:0/20:4)
# PE(18:1/22:6)
# TAG53:4FA18:0

# in ados but not in adults:
# LCER(26:1)
# LPC(20:1)
# PE(18:2/18:2)
# PE(O16:0/22:4)
# PE(O16:0/18:1)
# PI(16:0/18:1)
# PI(16:0/18:2)

# metabo_info:
metaboinfo_ados = fread('data/metabodata_info_ados.txt')
metaboinfo_adults = fread('data/metabodata_info_adults.txt')

metaboinfo_ados = metaboinfo_ados %>%
  mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))
metaboinfo_adults = metaboinfo_adults %>%
  mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))

metaboinfo_adults %>% #metaboinfo_ados %>% 
  # filter(metabolite_id %in% assoc_res_ados$metabolite_id) %>%
  filter(metabolite_id %in% assoc_res_adults$metabolite_id) %>%
  filter(class=="Amino acids") %>%
  arrange(metabolite_id) %>%
  slice(1:10) %>%
  pull(metabolite_id)

only_in_adults = setdiff(assoc_res_adults$metabolite_id,assoc_res_ados$metabolite_id)
only_in_ados = setdiff(assoc_res_ados$metabolite_id,assoc_res_adults$metabolite_id)
length(only_in_adults)#17
length(only_in_ados)#7
in_both = intersect(assoc_res_adults$metabolite_id,assoc_res_ados$metabolite_id)#831
assoc_res_adults %>% filter(metabolite_id %in% only_in_adults) %>% select(c(1,which(str_detect(names(assoc_res_adults),"P_"))))
assoc_res_ados %>% filter(metabolite_id %in% only_in_ados) %>% select(c(1,which(str_detect(names(assoc_res_ados),"P_"))))
