# rm(list=ls())
#*****************************************************************************#
# updated: 2026-10-07 to address NatCommuns reviewers' comments
# load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file= "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
# REGENIE-input adjusted for BMI
#*****************************************************************************#
select=dplyr::select
describe = psych::describe

# initialize ----
if(file.exists(project_specific_file)){
  source(project_specific_file)
}else{
  stop("\nNo project_specific_file exists.\n")
}
setwd(wd)

load("results/pca_brain_adjBMI_final_imputePCA.Rdata",verbose = T)#pca_brain_final(for adjBMI)
#brainPCA_imputed = PCA(averaged_data%>% select(-1)):unadj
load('results/brainPCA_imputed_unadjBMI_MICE.Rdata', verbose = T)
brainPCA_imputed$eig
brainPCA_imputed$var$coord

# boxplot: WMPC1-loadings ----
base.size=10
box_WM.PC1 =
  brainPCA_imputed$var$coord %>% data.frame %>% 
  rownames_to_column(var = 'pheno') %>% 
  mutate(pheno = str_split(pheno,"[.]",simplify = T)[,1]) %>%
  select(pheno, WM.PC1 = Dim.1) %>%
  mutate(pheno = factor(pheno,levels=c('Volume','T1wSI','MTR')), 
         WM.PC1 = -(1)*WM.PC1) %>% #flip
  ggplot(aes(x=`pheno`, y=WM.PC1)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = base.size) + 
  geom_hline(yintercept = 0, color = 'grey75',linewidth=0.5,alpha=0.5) +
  xlab(NULL)+ylab(NULL) + #+ ylab('loading') +   ggtitle('MetabolitePC1')
  scale_y_continuous( minor_breaks = seq(-1, 1, 0.25) ) + 
  geom_point(size=2, shape=2, alpha=0.5, stroke=1)+
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = 'grey80'),
        panel.grid.major = 
          element_line(size = 0.5, linetype = 'solid',
                       colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(angle=-45,hjust=0,vjust=0.5),
        axis.text = element_text(size=10)) + ggtitle(NULL);box_WM.PC1 

# boxplot: metaboPC1-loadings ----
Essential_AA = 'His
Ile
Leu
Val
Phe'
Essential_AA = unlist(str_split(Essential_AA,"\n"))

Other_AA = 'Ala
Gln
Tyr
Gly'
Other_AA = unlist(str_split(Other_AA,"\n"))

Cell_energy='Glucose
Lactate
Pyruvate
Citrate
Glycerol
Acetate
Acetoacetate
bOHbutyrate'
Cell_energy = unlist(str_split(Cell_energy,"\n"))

# head(pca.adj.metabo$var$coord[,c(1,3)])
# pca.loadings = pca.adj.metabo$var$coord[,c(1,3)]
# pca.loadings = data.frame(metabolite_id = rownames(pca.loadings),pca.loadings)
pca.loadings_adjBMI = fread('results/pca.loadings_adjBMI_nov2026.txt') %>%
  rename(metaboPC1_adjBMI= Dim.1,metaboPC3_adjBMI= Dim.3)
pca.loadings_adjBMI %>% head
pca.loadings_unadjBMI = fread('assoc_WMphenos_with_lipidPC_loadings_525metabos_2026-03-24.txt')

head(pca.loadings_unadjBMI)
pca.loadings_all = pca.loadings_unadjBMI %>% left_join(pca.loadings_adjBMI)
p1 = pca.loadings_all %>% ggplot(aes(x=metaboPC1,y=metaboPC1_adjBMI)) + geom_point()
p2 = pca.loadings_all %>% ggplot(aes(x=metaboPC3,y=metaboPC3_adjBMI)) + geom_point()
p1+p2 + plot_layout(ncol=2) & geom_abline(slope = 1, intercept = 0) & stat_cor()

## extract only the metPC1 and metPC3 unadjBMI ----
pca.loadings = pca.loadings_all %>% 
  # take BMI-adjusted loadings
  select(metabolite_id,class,`lipid class`,metaboPC1=metaboPC1,metaboPC3=metaboPC3)
ind1 = pca.loadings$metabolite_id %in% Essential_AA
ind2 = pca.loadings$metabolite_id %in% Other_AA
ind3 = pca.loadings$metabolite_id %in% Cell_energy

pca.loadings$class[ind1] <- 'NonL'#"Essential_AA"
pca.loadings$class[ind2] <- 'NonL'#"Other_AA"
pca.loadings$class[ind3] <- 'NonL'#"Cell_energy"

names(pca.loadings) = str_replace(names(pca.loadings),'Dim.',"metaboPC")
lipid.class = c('NeutL',
                'GPC',
                'SL',
                'NonL')#neutral lipids
nclass = length(lipid.class)*1
cols = c('#FFD166','#06D6A0','#118AB2','#8f00ff')#,'#EF476F',
class2.levels = c("Neutral lipids",
                  "Glycerophospholipids",
                  "Sphingolipids",
                  "Non-lipid metabolites")
names(cols) = class2.levels
metabo_info = fread('data/metabodata_info_ados.txt')
table(metabo_info$class)

Cell_energy = unlist(str_split(Cell_energy,"\n"))
ind1 = res_WM1$metabolite_id %in% Essential_AA
ind2 = res_WM1$metabolite_id %in% Other_AA
ind3 = res_WM1$metabolite_id %in% Cell_energy

res_WM1$class[ind1] <- 'NonL'#"Essential_AA"
res_WM1$class[ind2] <- 'NonL'#"Other_AA"
res_WM1$class[ind3] <- 'NonL'#"Cell_energy"

class.levels ='TAG
DAG
MAG
CE
PC
PE
PI
LPC
LPE
SM
CER
HCER
LCER
DCER
NonL'
# Essential_AA
# Other_AA
# Cell_energy'
class.levels = unlist(str_split(class.levels,"\n"))
#
# range(pca.loadings%>% select(-1))
pca.loadings = pca.loadings %>% 
  mutate(class = factor(class, levels=class.levels),
         `lipid class` = factor(`lipid class`,levels = class2.levels))
pca.loadings %>% filter(!is.na(metaboPC1)) %>%
  select(class,`lipid class`) %>% table(useNA='a') %>% sum#525

top_loaded = pca.loadings %>%  filter(!is.na(metaboPC1)) %>%
  filter(class=="TAG") %>% arrange(metaboPC1) %>% slice(1:2) %>%
  bind_rows(pca.loadings %>%  filter(!is.na(metaboPC1)) %>%
              filter(class=="TAG") %>% arrange(-metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="DAG") %>% arrange(metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="DAG") %>% arrange(-metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="CE") %>% arrange(-metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="PC") %>% arrange(-metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="PE") %>% arrange(-metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="LPC") %>% arrange(metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="LPE") %>% arrange(metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="SM") %>% arrange(-metaboPC1) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC1)) %>%
              filter(class=="NonL") %>% arrange(metaboPC1) %>% slice(1:2)) %>%
  select(-metaboPC3)
top_loaded
col_ind = pca.loadings %>% filter(!is.na(metaboPC1)) %>% select(`lipid class`) %>% table(useNA='a') >0
metaboPC1 = pca.loadings %>% filter(!is.na(metaboPC1)) %>%
  mutate(class = factor(class,levels = class.levels)) %>%
  mutate(`lipid class` = factor(`lipid class`,levels = class2.levels)) %>%
  # mutate(`lipid class` = ifelse(`lipid class`=="NonL","Non-lipid metabolites",`lipid class`)) %>%
  # res_WM1 %>% 
  # right_join(pca.loadings, join_by(metabolite_id)) %>% 
  arrange(class,`lipid class`) %>% 
  ggplot(aes(x=class, y=metaboPC1, color=`lipid class`)) +
  theme_bw(base_size = base.size) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=`lipid class`), width=0.15, size=1, shape=2, alpha=0.5)+
  # geom_hline(yintercept = c(-0.4,0.4), color="red",linetype=3) + 
  scale_color_manual(values=cols[col_ind],name="Lipid classes",
                     labels = class2.labels,
                     guide = guide_legend(override.aes = list(color=cols[col_ind]) ))+ # Edit legend title and label
  geom_hline(yintercept = 0, color = 'grey75',linewidth=0.5,alpha=0.75) +
  xlab(NULL)+ylab(NULL) + #+ ylab('loading') +   ggtitle('MetabolitePC1')
  scale_y_continuous( minor_breaks = seq(-1, 1, 0.2) ) + 
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = 'grey80'),
        panel.grid.major = 
          element_line(size = 0.5, linetype = 'solid',
                       colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(angle=-45,hjust=0,vjust=0.5),
        axis.text = element_text(size=10)) + ggtitle(NULL) 

# merge boxplots ----
# metaboPC1 = metaboPC1 + geom_text_repel(data = top_loaded ,aes(label = metabolite_id,x = class),size = 2,color='black')
p = metaboPC1+ ylab('')+coord_cartesian(ylim=c(-0.3,1)) +
  (box_WM.PC1 + ylab('')+coord_cartesian(ylim=c(-1,1))) +
  plot_layout(ncol=2, widths = c(10,3)) &
  theme(legend.position='none',
        plot.margin = margin(r = 0.5, t=0.65,unit = "cm")
        # plot.margin = margin(t=0.33, b = 0.125, l=0.125, r = 1, unit = "cm")
  );p
ggsave("PCloading_WMPC1_metPC1_525metabos_unadjBMI_2026-11-13.png",
       plot=p,height=10,width=20,unit='cm',dpi = 300)#this 'ggsave' does not work

top_loaded %>% arrange(-metaboPC1)

## write dataset for this figure
pca.loadings %>% select(metabolite_id,class, `lipid class`, metaboPC1) %>%
  filter(!is.na(metaboPC1)) %>%
  mutate(class = factor(class,levels = class.levels)) %>%
  mutate(`lipid class` = factor(`lipid class`,levels = class2.levels)) %>%
  arrange(class,`lipid class`) %>%
  write_tsv('results/[manuscript]_dataset_PCloading_WMPC1_metPC1_525metabos_unadjBMI_2026-11-13.txt')

# get % variance explained by compVar1_unadjBMI ----
colnames(brainPCA_imputed$ind$coord) <- str_replace(colnames(brainPCA_imputed$ind$coord),"Dim.","wmPC")
colnames(brainPCA_imputed$ind$coord) <- paste(colnames(brainPCA_imputed$ind$coord),"adjBMI",sep="_")
##
PC_merge =  %>% data.frame
pca_PC1_na.omit = PCA(na.omit(PC_merge %>% select(PC1.WM,PC1.Met)),graph = F)

