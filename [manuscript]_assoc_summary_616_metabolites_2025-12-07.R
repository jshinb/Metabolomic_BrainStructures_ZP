rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

PC_loadings_616 = readxl::read_xlsx("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/outputs/MetabolitePC_loadings_616_2025-05-26_mod_nov14.xlsx")
PC_loadings_616 = data.table(PC_loadings_616)

plot.df.volcano_main = fread("results/plot.df.volcano_main_2025-12-02.tsv")
head(plot.df.volcano_main)
class1.levels = unlist(str_split('TAG,DAG,MAG,CE,PC,PE,PI,LPC,LPE,SM,CER,HCER,LCER,DCER,Other',','))
other.class = c("Amino acids","Glycolysis related metabolites","Inflammation","Ketone bodies")
plot.df.volcano_main = plot.df.volcano_main %>% 
  mutate(class1 = ifelse(class %in% other.class, "Other",class)) %>%
  mutate(class1 = factor(class1,levels=class1.levels))
assoc_res_part = (plot.df.volcano_main %>% filter(trait == "lobar.vol.WM.adjICV_fam", generation=="ados"))
tab_analyzed = as.data.frame.matrix(table(assoc_res_part %>% pull(class1)))
assoc_res_part = assoc_res_part %>%
  mutate(assoc.3.pheno = ifelse(metabolite_id %in% PC_loadings_616$metabolite_id,'wm.assoc','wm.unassoc'))

tab_analyzed = table(assoc_res_part %>% select(class1,assoc.3.pheno))
tab_analyzed = as.data.frame.matrix(tab_analyzed)
tab_analyzed = tab_analyzed%>% 
  mutate(all = wm.assoc + wm.unassoc) %>% 
  mutate(perc.associated = wm.assoc/all) %>%
  mutate(perc.class.assoc = wm.assoc/sum(wm.assoc),
         perc.class.all = all/sum(all))
beta.T1wSI <- beta.MTR <- beta.Vol <- c()
sd.T1wSI <- sd.MTR <- sd.Vol <- c()
for(i in class1.levels){
  beta.T1wSIi = plot.df.volcano_main %>% filter(generation=="ados",trait=="NormWM_fam",class1 == i) %>% select(beta) %>% describe %>% pull(mean)
  beta.MTRi = plot.df.volcano_main %>% filter(generation=="ados",trait=="MTR.WM_fam",class1 == i) %>% select(beta) %>% describe %>% pull(mean)
  beta.Voli = plot.df.volcano_main %>% filter(generation=="ados",trait=="lobar.vol.WM.adjICV_fam",class1 == i) %>% select(beta) %>% describe %>% pull(mean)
  sd.T1wSIi = plot.df.volcano_main %>% filter(generation=="ados",trait=="NormWM_fam",class1 == i) %>% select(beta) %>% describe %>% pull(sd)
  sd.MTRi = plot.df.volcano_main %>% filter(generation=="ados",trait=="MTR.WM_fam",class1 == i) %>% select(beta) %>% describe %>% pull(sd)
  sd.Voli = plot.df.volcano_main %>% filter(generation=="ados",trait=="lobar.vol.WM.adjICV_fam",class1 == i) %>% select(beta) %>% describe %>% pull(sd)
  beta.T1wSI <- c(beta.T1wSI ,beta.T1wSIi)
  sd.T1wSI <- c(sd.T1wSI ,sd.T1wSIi)
  beta.MTR <- c(beta.MTR ,beta.MTRi)
  sd.MTR <- c(sd.MTR ,sd.MTR)
  beta.Vol <- c(beta.Vol ,beta.Voli)
  sd.Vol <- c(sd.Vol ,sd.Voli)
  }
tab_analyzed = tab_analyzed %>% mutate(BetaT1wSI = beta.T1wSI, BetaMTR = beta.MTR, BetaVol = beta.Vol)
tab_analyzed = tab_analyzed %>% select(AnalyzedN = all, AssociatedN = wm.assoc, perc.Associated = perc.associated,
                                       BetaT1wSI,BetaMTR,BetaVol)
