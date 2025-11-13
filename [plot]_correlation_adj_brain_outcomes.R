#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#
#*****************************************************************************#

# 0. initialize ---------------------------------------------------------------
source(project_specific_file)
setwd(wd)

# 
load('results/ret_resid_unadjBMI_MTR_WM_ados_2025-04-08.Rdata')#mtr
mtrWM_adjAge = ret_resid_unadjBMI$`LobarWM_Z MTR`$resid_out
rm(ret_resid_unadjBMI)
head(mtrWM_adjAge)

load('results/ret_resid_unadjBMI_unadjICV_ados_2025-03-26.Rdata')#t1wSI
NormWM_adjAge = ret_resid_unadjBMI$NormWM$resid_out
head(NormWM_adjAge)
rm(ret_resid_unadjBMI)

load('results/ret_resid_unadjBMI_lobarVolWM_ados_2025-04-04.Rdata')#
WMvol_adjAge = ret_resid_unadjBMI$lobar.vol.WM.adjICV$resid_out
head(WMvol_adjAge)
rm(ret_resid_unadjBMI)

braindata_adjAge = NormWM_adjAge %>% 
  left_join(mtrWM_adjAge)%>% 
  left_join(WMvol_adjAge)

dim(braindata_adjAge)
head(braindata_adjAge)
