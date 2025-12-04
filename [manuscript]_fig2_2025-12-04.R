rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

# figure 2 ----
# Example data (replace with your actual data)
# Let's create a sample data frame with 14 categories (A to N)
set.seed(123)
data <- data.frame(
  category = factor(sample(LETTERS[1:14], 100, replace = TRUE))
)

# Summarize the data to count 'n' for each category
summary_data <- data %>%
  count(category) %>%
  mutate(label_text = paste("n =", n)) # Create a label column

# A list of 14 distinct colors
# cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
# 
# my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
#                "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2")

my_colors <- c("#FFD166", "#118AB2", "#06D6A0", "#F0E442", "#0072B2", "#D55E00", "#EF476F",
               "#1F77B4", "#FF7F0E", "#073B4C", "#D62728", "#8f00ff", "#8C564B", "#E377C2")

ggplot(summary_data, aes(x = category, y = n, fill = category)) + # Map 'category' to fill for different colors
  geom_col() + # Use geom_col for pre-counted data
  geom_text(aes(label = label_text), vjust = -0.5, size = 3.5) + # Add count labels above bars
  scale_fill_manual(values = my_colors) + # Manually apply the 14 colors
  labs(
    title = "Bar Plot with 14 Colors and Counts",
    x = "Category",
    y = "Count (n)"
  ) +
  theme_minimal() + # Use a clean theme
  theme(legend.position = "none") # Hide the legend if it's redundant

add_class2 = function(volcano_main_all){
  metabodata_info = fread("data/metabodata_info_ados.txt")
  
  if(sum(names(volcano_main_all)=='class')==1){
    volcano_main_all = volcano_main_all %>% 
      dplyr::select(-class) 
  }
  volcano_main_all = volcano_main_all %>% 
    mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>%
    left_join(
      metabodata_info %>%
        mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>% 
        dplyr::select(metabolite_id,class)
    )
  # class2.levels = c("CE", 
  #                   "TAG, DAG, or MAG", 
  #                   "PC, PE, or PI",
  #                   "LPC or LPE", 
  #                   "SM or ceramides",
  #                   "Nightingale")
  class2.levels = c("Cholesteryl esters",
                    "Acylglycerols",
                    "Phospholipids",
                    "Lysophospholipids",
                    "Sphingolipids",
                    "Inflammation/AA/metabolism")
  volcano_main_all[['class2']] <- class2.levels[1]
  volcano_main_all = volcano_main_all %>%
    mutate(class2 = ifelse(class %in% c('TAG','DAG','MAG'),class2.levels[2],class2)) %>%
    mutate(class2 = ifelse(class %in% c('PC','PE','PI'),class2.levels[3],class2)) %>%
    mutate(class2 = ifelse(class %in% c('LPC','LPE'),class2.levels[4],class2)) %>%
    mutate(class2 = ifelse(class %in% c('SM','CER',"DCER","HCER","LCER"),class2.levels[5],class2)) %>%
    mutate(class2 = ifelse(class %in% c('Amino acids','Glycolysis related metabolites','Inflammation','Ketone bodies'),
                           class2.levels[6],class2))
  table(volcano_main_all$class,volcano_main_all$class2,useNA = 'a')
  volcano_main_all = volcano_main_all %>%
    mutate(`lipid class` = factor(class2,levels=class2.levels)) 
  class1_cols <- c("#FFD166", "#118AB2", "#06D6A0", "#F0E442", "#0072B2", "#D55E00", "#EF476F",
                               "#1F77B4", "#FF7F0E", "#073B4C", "#D62728", "#8f00ff", "#8C564B", "#E377C2")
  
  class2_cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
  class2.labels = paste(names(table(volcano_main_all$`lipid class`))," (n=",
                        table(volcano_main_all$`lipid class`),")",sep='')
  ret=(list(volcano_main_all=volcano_main_all,class2.labels=class2.labels,cols=cols))
  ret
}
