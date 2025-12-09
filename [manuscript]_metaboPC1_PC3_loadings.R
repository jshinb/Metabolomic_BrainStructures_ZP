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
# my_colors <- c("#FFD166", "#118AB2", "#06D6A0", "#F0E442", "#0072B2", "#D55E00", "#EF476F",
#                "#1F77B4", "#FF7F0E", "#073B4C", "#D62728", "#8f00ff", "#8C564B", "#E377C2",
#                "#C0C0C0")
cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
# load association results:

PC_loadings_616 = readxl::read_xlsx("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/outputs/MetabolitePC_loadings_616_2025-05-26_mod_nov14.xlsx")
PC_loadings_616 = data.table(PC_loadings_616)

#8. PCA-loadings for lipid-PC -------------------------------------------------
pca.lipids.plot.df = PC_loadings_616
head(pca.lipids.plot.df)

class1.levels = unlist(str_split('TAG,DAG,MAG,CE,PC,PE,PI,LPC,LPE,SM,CER,HCER,LCER,DCER,Other',','))
other.class = c("Amino acids","Glycolysis related metabolites","Inflammation","Ketone bodies")
pca.lipids.plot.df = pca.lipids.plot.df %>% 
  mutate(class1 = ifelse(class %in% other.class, "Other",class)) %>%
  mutate(class1 = factor(class1,levels=class1.levels))
class2.levels = c("Cholesteryl esters",
                  "Acylglycerols",
                  "Phospholipids",
                  "Lysophospholipids",
                  "Sphingolipids",
                  "Inflammation/AA/metabolism")
class1.labels = paste(names(table(pca.lipids.plot.df$class1))," (n=",
                      table(pca.lipids.plot.df$class1),")",sep='')
class2.labels = paste(names(table(pca.lipids.plot.df$`lipid class`))," (n=",
                      table(pca.lipids.plot.df$`lipid class`),")",sep='')
# boxplots ----
boxplot.ymax = -0.9
box.n.size = 4
cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
my_colors = c('#FFD166','#FFD166','#FFD166',
              '#EF476F',
              '#118AB2','#118AB2','#118AB2',
              '#073B4C','#073B4C',
              '#06D6A0','#06D6A0','#06D6A0','#06D6A0','#06D6A0',
              '#8f00ff')

names(my_colors) = class1.levels
names(class1.labels) = class1.levels
my_colors = my_colors[sort(unique(pca.lipids.plot.df$class1))]
class1.labels = class1.labels[sort(unique(pca.lipids.plot.df$class1))]

## PC1 ----
metaboPC1 = pca.lipids.plot.df  %>% 
  mutate(`lipid class` = factor(`lipid class`,levels = class2.levels)) %>%
  arrange(`lipid class`) %>% 
  ggplot(aes(x=`lipid class`,y=PC1,color=`lipid class`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=`lipid class`), width=0.15, size=1, shape=2, alpha=0.5)+
  geom_hline(yintercept = c(-0.4,0.4), color="grey75",linetype=2) + 
  scale_color_manual(values=cols,name="Lipid classes",
                     labels = class2.labels,
                     guide = guide_legend(override.aes = list(color=cols) ))+ # Edit legend title and label
  geom_hline(yintercept = 0, color = 'grey25') +
  xlab(NULL)+ylab(NULL)#+ ylab('loading') +   ggtitle('MetabolitePC1')

metaboPC1 + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=14),
        # axis.text.x = element_text(angle=-25,hjust=0,vjust=0.5)
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=13),
        legend.position = 'none',
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) + 
  coord_cartesian(ylim=c(-1,1))

## with 15 classes
metaboPC1 = 
  pca.lipids.plot.df  %>% 
  mutate(`lipid class` = factor(`lipid class`,levels = class2.levels)) %>%
  arrange(class1) %>% 
  ggplot(aes(x=class1,y=PC1,color=class1)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=class1), width=0.15, size=1, shape=2, alpha=0.5)+
  geom_hline(yintercept = c(-0.4,0.4), color="grey75",linetype=2) + 
  scale_color_manual(values=my_colors,name="Lipid classes",
                     labels = class1.labels,
                     guide = guide_legend(override.aes = list(color=my_colors) ))+ # Edit legend title and label
  geom_hline(yintercept = 0, color = 'grey25') +
  xlab(NULL)+ylab(NULL)#+ ylab('loading') +   ggtitle('MetabolitePC1')

plot_sample_sizes <- 
  pca.lipids.plot.df %>% #mutate(x=class1) %>%
  select(class1,PC1) %>%
  group_by(class1) %>%
  summarise(n = n())

# metaboPC1 = metaboPC1 + 
  # geom_text(data = plot_sample_sizes, 
  #           aes(x=class1, y=boxplot.ymax, label = paste0("(n=",n,")")),
  #           vjust = 1.5,  # Adjust vertical position of text
  #           size = box.n.size)+
  
##PC3 ----
metaboPC3 =
  pca.lipids.plot.df  %>%
  mutate(`lipid class` = factor(`lipid class`,levels = class2.levels)) %>%
  arrange(class1) %>%
  ggplot(aes(x=class1,y=PC3,color=class1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=class1), width=0.15, size=1, shape=2, alpha=0.5)+
  geom_hline(yintercept = c(-0.4,0.4), color="grey75",linetype=2) +
  scale_color_manual(values=my_colors,name="Lipid classes",
                     labels = class1.labels,
                     guide = guide_legend(override.aes = list(color=my_colors) ))+ # Edit legend title and label
  geom_hline(yintercept = 0, color = 'grey25') +
  xlab(NULL)+ylab(NULL)#+ ylab('loading') +   ggtitle('MetabolitePC3')
plot_sample_sizes <-
  pca.lipids.plot.df %>% #mutate(x=class1) %>%
  select(class1,PC3) %>%
  group_by(class1) %>%
  summarise(n = n())
# metaboPC3 = metaboPC3 +
  # geom_text(data = plot_sample_sizes,
  #           aes(x=class1, y=boxplot.ymax, label = paste0("(n=",n,")")),
  #           vjust = 1.5,  # Adjust vertical position of text
  #           size = box.n.size)+

# generating png file ----
(metaboPC1+metaboPC3 ) + plot_layout(ncol=2) & 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=14),
        axis.text.x = element_text(angle=30,hjust=1,vjust=1),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size=13),
        legend.position = 'none',
        panel.border = element_rect(colour = "grey85", fill=NA, linewidth=1)) &
  coord_cartesian(ylim=c(-0.6,1)) &  
  scale_x_discrete(labels= class1.labels)

ggsave('outputs/PC1_PC3_loadings.png',
       width=11*0.8,height = 5*0.8,units='in',dpi=300)
