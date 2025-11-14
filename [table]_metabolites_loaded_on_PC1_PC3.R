rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

pc_res = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/results/pc_res_PC1_PC10_wi_classes.tsv")
head(pc_res)

ind = abs(pc_res$Dim.1) >0.4; print(sum(ind))#461
# ind = ind | abs(pc_res$Dim.2) >=0.4; print(sum(ind))#462
ind = ind | abs(pc_res$Dim.3) >0.4; print(sum(ind))#464

#     metabolite_id  class     Dim.1     Dim.2
#            <char> <char>     <num>     <num>
# 1: TAG60:12FA22:6    TAG 0.3700744 0.4416861

write_tsv( pc_res %>% filter(ind) %>% select(c(1:3,5) ),
           file.path(wd,"results","results_PC123_pc_res_2025apr04.tsv"))
#   metabolite_id        class     Dim.1       Dim.3
#          <char>       <char>     <num>       <num>
# 1:         GlycA Inflammation 0.6845942 -0.07396909
# 2:      CE(14:0)           CE 0.5146809  0.24409194
# 3:      CE(14:1)           CE 0.6785204  0.15994002

# Excel file sent via email: "pc_res_PC1_and_PC3_loading_0.4_2025-11-14.xlsx" (on Nov 14, 2025)