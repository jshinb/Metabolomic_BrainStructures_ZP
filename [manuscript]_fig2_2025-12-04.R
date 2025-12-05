rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

# figure 2: scatter plots of association estimates for WM phenotypes in adolescents ----
# A list of 14 distinct colors
# cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
# 
# my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
#                "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2")

my_colors <- c("#FFD166", "#118AB2", "#06D6A0", "#F0E442", "#0072B2", "#D55E00", "#EF476F",
               "#1F77B4", "#FF7F0E", "#073B4C", "#D62728", "#8f00ff", "#8C564B", "#E377C2")
# load association results:

plot.df.volcano_main = fread("results/plot.df.volcano_main_2025-12-02.tsv")
as.data.frame.matrix(table(plot.df.volcano_main %>% select(class,class2))) %>% dim

dput(rownames(as.data.frame.matrix(table(plot.df.volcano_main %>% select(class,class2)))))

class1.levels = unlist(str_split('TAG,DAG,MAG,CE,PC,PE,PI,LPC,LPE,SM,CER,HCER,LCER,DCER,Other',','))
other.class = c("Amino acids","Glycolysis related metabolites","Inflammation","Ketone bodies")

plot.df.volcano_main = plot.df.volcano_main %>%
  mutate(class1 = ifelse(class %in% other.class,"Other",class)) %>%
  mutate(class1 = factor(class1,levels = class1.levels))
table(plot.df.volcano_main %>% select(class1,class))
