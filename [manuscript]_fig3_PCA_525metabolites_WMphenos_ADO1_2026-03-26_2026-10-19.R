#*****************************************************************************#
# updated: 2026-10-07 to address NatCommuns reviewers' comments
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file= "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
# REGENIE-input adjusted for BMI
#*****************************************************************************#
select=dplyr::select
describe = psych::describe

# 0. initialize --describe()# 0. initialize ----
source(project_specific_file)
setwd(wd)
res_WM1 = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/TabS1_2026-03-17.txt")
head(res_WM1 )
res_WM1 = res_WM1 %>% filter(metaboliteid!='GlycA')

ind_ado = res_WM1$fdr.p2_WM.T1wSI_ados<0.05;print(sum(ind_ado))#651
ind_ado = ind_ado | res_WM1$fdr.p2_WM.MTR_ados <0.05;print(sum(ind_ado))#711# 
ind_ado = ind_ado | res_WM1$fdr.p2_WM.Vol_ados <0.05;print(sum(ind_ado))#711# after removing 'GlycA'

ind_ado1 = res_WM1$fdr.p2_WM.T1wSI_ados<0.05;print(sum(ind_ado1))#650
ind_ado1 = ind_ado1 & res_WM1$fdr.p2_WM.MTR_ados <0.05;print(sum(ind_ado1))#711# 
ind_ado1 = ind_ado1 & res_WM1$fdr.p2_WM.Vol_ados <0.05;print(sum(ind_ado1))#711# after removing 'GlycA'

ind_ado2 = res_WM1$fdr.p2_WM.T1wSI_ados<0.1;print(sum(ind_ado2))#673
ind_ado2 = ind_ado2 & res_WM1$fdr.p2_WM.MTR_ados <0.1;print(sum(ind_ado2))#618
ind_ado2 = ind_ado2 & res_WM1$fdr.p2_WM.Vol_ados <0.1;print(sum(ind_ado2))#564

ind_adult = res_WM1$fdr.p2_WM.T1wSI_adults<0.05;print(sum(ind_adult))#41
ind_adult = ind_adult | res_WM1$fdr.p2_WM.MTR_adults <0.05;print(sum(ind_adult))#63
ind_adult = ind_adult | res_WM1$fdr.p2_WM.Vol_adults <0.05;print(sum(ind_adult))#73

ind_adult2 = res_WM1$fdr.p2_WM.T1wSI_adults<0.1;print(sum(ind_adult2))#85
ind_adult2 = ind_adult2 & res_WM1$fdr.p2_WM.MTR_adults <0.1;print(sum(ind_adult2))#3
ind_adult2 = ind_adult2 & res_WM1$fdr.p2_WM.Vol_adults <0.1;print(sum(ind_adult2))#0

res_WM1 %>% filter(ind_ado) %>% select(beta_WM.T1wSI_ados,beta_WM.MTR_ados,beta_WM.Vol_ados) %>% cor(use='p')
res_WM1 %>% filter(ind_ado2) %>% select(beta_WM.T1wSI_ados,beta_WM.MTR_ados,beta_WM.Vol_ados) %>% cor(use='p')
res_WM1 %>% filter(ind_adult) %>% select(beta_WM.T1wSI_adults,beta_WM.MTR_adults,beta_WM.Vol_adults) %>% cor(use='p')

#
metabo_info = fread('data/metabodata_info_ados.txt')

# ----@
lipid.class = c('NeutL',
                'GPC',
                'SL',
                'NonL')#neutral lipids
nclass = length(lipid.class)*1
# cols = c('#E1AD01','#06D6A0','#118AB2','#8f00ff')#,'#EF476F',
cols <- c("#FFCE1B", "#06D6A0", "#118AB2", "#8f00ff")# barplot
class2.levels = c("Neutral lipids",
                  "Glycerophospholipids",
                  "Sphingolipids",
                  "Non-lipid metabolites")
metabo_info = fread('data/metabodata_info_ados.txt')
metabo_info = metabo_info %>% mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))
metabo_info = metabo_info%>% filter(metabolite_id %in% res_WM1$metaboliteid)#n=830
res_WM1 = res_WM1 %>% select(-`lipid class`)
res_WM1 = res_WM1 %>% rename(metabolite_id = metaboliteid) %>% left_join(metabo_info)
res_WM1 = res_WM1 %>%
  mutate(class2 = case_when(
    class %in% c('TAG','DAG','MAG','CE') ~ class2.levels[1],
    class %in% c('PC','PE','PI','LPC','LPE') ~ class2.levels[2],
    class %in% c('SM','CER',"DCER","HCER","LCER") ~ class2.levels[3],
    class %in% c('Amino acids','Glycolysis related metabolites',
                 'Inflammation',
                 'Ketone bodies') ~ class2.levels[4]
  )) %>% mutate(class2 = factor(class2, class2.levels)) %>%
  rename(`lipid class` = class2)
class2.labels = paste(names(table(res_WM1 %>% filter(ind_ado1) %>% pull(`lipid class`)))," (n=",
                      table(res_WM1 %>% filter(ind_ado1) %>% pull(`lipid class`)),")",sep='')
top_main = res_WM1 %>% 
  filter(
    beta_WM.T1wSI_ados<0,
    beta_WM.MTR_ados <0, 
    beta_WM.Vol_ados>0) %>%
  arrange(fdr.p2_WM.T1wSI_ados*fdr.p2_WM.MTR_ados*fdr.p2_WM.Vol_ados) %>%
  slice(1:10)

top_main = top_main %>% 
  bind_rows(res_WM1 %>% 
              filter(
                beta_WM.T1wSI_ados>0,
                beta_WM.MTR_ados >0, 
                beta_WM.Vol_ados<0) %>%
              arrange(fdr.p2_WM.T1wSI_ados*fdr.p2_WM.MTR_ados*fdr.p2_WM.Vol_ados) %>%
              slice(1:10))

p1 = res_WM1 %>%
  filter(ind_ado) %>% #n=711
  arrange(`lipid class`) %>% 
  ggplot(aes(x=beta_WM.Vol_ados,y=beta_WM.T1wSI_ados,
             color=`lipid class`)) + #
  xlab(bquote(hat(beta) ~ "_WM.Volume"))+
  ylab(bquote(hat(beta) ~ "_WM.T1wSI"))+
  # ylab(bquote(hat(beta) ~ "_MTR"))+
  scale_color_manual(values=cols,name="Lipid classes",
                     labels = class2.labels,
                     guide = guide_legend(override.aes = list(color=cols) )) + # Edit legend title and labels
  # stat_poly_line() + 
  # stat_poly_eq(use_label(c("eq", "R")), label.x = 1, label.y=1)+
  geom_point(alpha=0.6,shape=19) + 
  # geom_smooth(method = 'lm', formula=y~x) + 
  coord_cartesian(ylim=c(-0.35,0.35),xlim=c(-0.35,0.35)) +#600x550 
  # ggpubr::stat_cor(label.x.npc = 0.75) +
  # stat_summary(fun.data=mean_cl_normal) +
  # geom_abline(slope=-1, intercept=0, linetype=2) + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_vline(xintercept = 0, color="darkgrey") +
  scale_y_continuous( minor_breaks = seq(-1, 1, 0.1) ) + 
  geom_text_repel(
    data = top_main,
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha("lightgray",0.5),
    color = 'black',
    #fill = scales::alpha('white',0.5),
    size = 3) + #theme_bw()
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = 'grey80'),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

p1 + theme(legend.position = 'none')
ggsave('results/scatterplot_beta_WM.Vol_WM.T1wSI_sex-combined_ados_2026-03-20.png',
       plot = (p1 + theme(legend.position = 'none')), 
       width=6, height=5.5, units='in', dpi=300)
ggsave('results/scatterplot_beta_WM.Vol_WM.T1wSI_sex-combined_ados_2026-03-20_wi_legend.png',
       plot = (
         p1 + theme(
           legend.position = 'bottom', 
           legend.title = element_blank()
         )+
           guides(color = guide_legend(nrow = 2, byrow = TRUE))
       ), 
       width=12, height=5.5, units='in', dpi=300)

# 
p2 = res_WM1 %>%
  filter(ind_ado) %>% #n=525 n=711
  arrange(`lipid class`) %>% 
  ggplot(aes(x=beta_WM.Vol_ados,y=beta_WM.MTR_ados,
             color=`lipid class`)) + #
  xlab(bquote(hat(beta) ~ "_WM.Volume"))+
  ylab(bquote(hat(beta) ~ "_WM.MTR"))+
  # ylab(bquote(hat(beta) ~ "_MTR"))+
  scale_color_manual(values=cols,name="Lipid classes",
                     labels = class2.labels,
                     guide = guide_legend(override.aes = list(color=cols) )) + # Edit legend title and labels
  # stat_poly_line() + 
  # stat_poly_eq(use_label(c("eq", "R")), label.x = 1, label.y=1)+
  geom_point(alpha=0.5,shape=19) + 
  # geom_smooth(method = 'lm', formula=y~x) + 
  coord_cartesian(ylim=c(-0.35,0.35),xlim=c(-0.35,0.35)) +#600x550 
  # ggpubr::stat_cor(label.x.npc = 0.75) +
  # stat_summary(fun.data=mean_cl_normal) +
  # geom_abline(slope=-1, intercept=0, linetype=2) + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_vline(xintercept = 0, color="darkgrey") +
  geom_text_repel(
    data = top_main,
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha("lightgray",0.5),
    color = 'black',
    #fill = scales::alpha('white',0.5),
    size = 3) + #theme_bw()
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = 'grey80'),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

p2 + theme(legend.position = 'none')
ggsave('results/scatterplot_beta_WM.Vol_WM.MTR_sex-combined_ados_2026-03-20.png',
       plot = (p2 + theme(legend.position = 'none')), 
       width=6, height=5.5, units='in', dpi=300)

# T1wSI (x) vs. MTR (y)
p3 = res_WM1 %>%
  filter(ind_ado) %>% #n=711
  arrange(`lipid class`) %>% 
  ggplot(aes(x=beta_WM.T1wSI_ados,y=beta_WM.MTR_ados,
             color=`lipid class`)) + #
  xlab(bquote(hat(beta) ~ "_WM.T1wSI"))+
  ylab(bquote(hat(beta) ~ "_WM.MTR"))+
  # ylab(bquote(hat(beta) ~ "_MTR"))+
  scale_color_manual(values=cols,name="Lipid classes",
                     labels = class2.labels,
                     guide = guide_legend(override.aes = list(color=cols) )) + # Edit legend title and labels
  # stat_poly_line() + 
  # stat_poly_eq(use_label(c("eq", "R")), label.x = 1, label.y=1)+
  geom_point(alpha=0.5,shape=19) + 
  # geom_smooth(method = 'lm', formula=y~x) + 
  coord_cartesian(ylim=c(-0.35,0.35),xlim=c(-0.35,0.35)) +#600x550 
  # ggpubr::stat_cor(label.x.npc = 0.75) +
  # stat_summary(fun.data=mean_cl_normal) +
  # geom_abline(slope=-1, intercept=0, linetype=2) + 
  geom_hline(yintercept = 0, color="darkgrey") + 
  geom_vline(xintercept = 0, color="darkgrey") +
  geom_text_repel(
    data = top_main,
    aes(label = metabolite_id,hjust='outward'),# position = position_nudge_center(x=0.3,center_x = 0),
    position=position_nudge_center(
      #x = 0,#, y = -log10(0.05),
      direction = "split"
    ),
    max.time = 1, max.iter = 1e5, max.overlaps = 25,
    box.padding = 0.3,
    min.segment.length = 0,
    segment.color = scales::alpha("lightgray",0.5),
    color = 'black',
    #fill = scales::alpha('white',0.5),
    size = 3) + #theme_bw()
  theme(panel.background = element_rect(fill = scales::alpha("#f8edeb",0.2),
                                        #colour = scales::alpha("#f8edeb",0.2),
                                        size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = 'grey80'),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey95"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey95"),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

p3 + theme(legend.position = 'none')
ggsave('results/scatterplot_beta_WM.T1wSI_WM.MTR_sex-combined_ados_2026-03-20.png',
       plot = (p3 + theme(legend.position = 'none')), width=6, height=5.5, units='in', dpi=300)
res_WM1 %>% filter(metabolite_id == "CER(18:0)")


# clustering
library(pheatmap)
load('results/adj.metabo.Rdata',verbose = T)
mat = as.matrix(adj.metabo %>% select(all_of(res_WM1 %>% filter(ind_ado) %>% pull(metabolite_id))))
dim(mat)

library(factoextra)
corm = cor(mat,use="p")
fviz_nbclust(corm, hcut, method = "silhouette", hc_method = "complete",k.max = 20)
hclust.method="complete"
k=5
clustered_data <- pheatmap(corm,
                           main = paste0(k," clusters with age- and sex-adjusted metabolites"),
                           cluster_cols = T,
                           cutree_cols = k,
                           cutree_rows = k,
                           clustering_distance_cols = "euclidean",
                           clustering_distance_rows = "euclidean",
                           clustering_method = hclust.method,
                           fontsize = 3.5)
groups <- cutree(clustered_data$tree_col, k = k)
head(groups)
table(groups)

adj.metabo.df = data.frame(groups)
adj.metabo.df = data.frame(metabolite_id = rownames(adj.metabo.df),cluster = adj.metabo.df)
head(adj.metabo.df)

## use the imputed metabo generated for the 711 metabolites associated with at least 
## one WM phenotype at fdr-p<0.05, adjusted for multiple comparisons among
## 3 WM phenotypes in ADO and PAR groups (mar2026)
library(mice)
load('data/imputed_adj.metabo.Rdata')
cluster.length=length(imputed_metabolites)
for(i in 1:cluster.length){
  tmpi = mat[,adj.metabo.df %>% filter(groups==i) %>% pull(metabolite_id)]
  dim(tmpi)#41
  names.tmpi = colnames(tmpi)
  tmpi = data.frame(tmpi)
  imputed_data.metabo[[i]] <-
    complete(imputed_metabolites[[i]], action = "long") %>%
    group_by(.id) %>%
    summarise(across(everything(), mean)) %>%
    select(-.imp) # Remove the imputation index column
  names( imputed_data.metabo[[i]])[1] = "uniqueID"
  names( imputed_data.metabo[[i]])[-1] = names.tmpi
  if(i ==1 ){
    adj.metabo.imputed = imputed_data.metabo[[i]]
  }else{
    adj.metabo.imputed = adj.metabo.imputed %>%
      left_join(imputed_data.metabo[[i]])
  }
  cat(paste0('\n***********cluster',i," ends here.\n"))
}
## PCA of imputed metabo ----
adj.metabo.imputed.pca.input = as.matrix(adj.metabo.imputed %>% select(-1) )
rownames(adj.metabo.imputed.pca.input) = adj.metabo.imputed %>% pull(1)

library(FactoMineR)
pca.adj.metabo = FactoMineR::PCA(adj.metabo.imputed.pca.input,ncp = 20)
plot(pca.adj.metabo,choix = 'var',axes = c(1,2))
plot(pca.adj.metabo,choix = 'var',axes = c(3,4))

rm(ind);ind = res_WM1$fdr.p2_WM.Vol_ados<0.05;print(sum(ind))#548
ind = ind & res_WM1$fdr.p2_WM.T1wSI_ados<0.05;print(sum(ind))#542
ind = ind & res_WM1$fdr.p2_WM.MTR_ados<0.05;print(sum(ind))#525

res_WM1 %>% filter(ind) %>% select(contains("unique"),contains("ados"))
sel_metabo = res_WM1 %>% filter(ind) %>% pull(metabolite_id)

adj.metabo.imputed.pca.input = as.matrix(adj.metabo.imputed %>% select(-1) %>% select(all_of(sel_metabo)) )
rownames(adj.metabo.imputed.pca.input) = adj.metabo.imputed %>% pull(1)
dim(adj.metabo.imputed.pca.input)#final: 974*525
head(adj.metabo.imputed.pca.input[,1:5])

library(FactoMineR)
pca.adj.metabo = FactoMineR::PCA(adj.metabo.imputed.pca.input,ncp = 20)
plot(pca.adj.metabo,choix = 'var',axes = c(1,2))
plot(pca.adj.metabo,choix = 'var',axes = c(3,4))
head(pca.adj.metabo$ind$coord)
## 
dt_metabo = 
  rownames_to_column(
    data.frame(pca.adj.metabo$eig) %>% slice(1:3) %>% select(2),var ='Dim')
#550*150
graphics.off();
gridExtra::grid.table(dt_metabo %>% mutate(across(-1,~round(., digits=2))))

PC2.Met.scores = pca.adj.metabo$ind$coord %>% data.frame %>%
  rownames_to_column(var='uniqueID') %>% select(uniqueID,PC2.Met= Dim.2)

## PCA of brain data ----
# load covdata for adjusting for sex
covdata = fread(file.path(wd,"data","covdata_ados.txt"))

## WM: age- and/or ICV-adjusted
load('results/ret_resid_unadjBMI_MTR_WM_ados_2025-04-08.Rdata')#mtr
mtrWM_adjAge = ret_resid_unadjBMI$`LobarWM_Z MTR`$resid_out
rm(ret_resid_unadjBMI)
head(mtrWM_adjAge)
mtrWM_adjAge %>% describe#n=784;n=773

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
  left_join(WMvol_adjAge) %>%
  filter(uniqueID %in% adj.metabo$uniqueID) %>%
  select(-contains('gkin')) %>%
  rename(
    Volume = lobar.vol.WM.adjICV_fam,
    T1wSI = NormWM_fam,
    MTR = `LobarWM_Z MTR_fam`) %>%
  relocate(Volume,.after="uniqueID") %>%
  left_join(covdata,join_by(uniqueID))
adj.metabo = adj.metabo %>% filter(uniqueID %in% braindata_adjAge$uniqueID)

# adjust brain for sex ---
brain.adjAgeSex = braindata_adjAge %>% select(1)
fit = lm(Volume ~ Sex, data=braindata_adjAge, na.action = na.exclude)
identical(names(resid(fit)),brain.adjAgeSex %>% pull(1))
brain.adjAgeSex = data.frame(brain.adjAgeSex, Volume = resid(fit))

fit = lm(T1wSI ~ Sex, data=braindata_adjAge, na.action = na.exclude)
identical(names(resid(fit)),brain.adjAgeSex %>% pull(1))
brain.adjAgeSex = data.frame(brain.adjAgeSex, T1wSI=resid(fit))

fit = lm(MTR ~ Sex, data=braindata_adjAge, na.action = na.exclude)
identical(names(resid(fit)),brain.adjAgeSex %>% pull(1))
brain.adjAgeSex = data.frame(brain.adjAgeSex, MTR = resid(fit))
describe(brain.adjAgeSex)

#join with metabo
brain.adjAgeSex = brain.adjAgeSex %>% left_join(adj.metabo)

# impute missiong data
library(mice)
brain_missing = md.pattern(brain.adjAgeSex%>%select(Volume,T1wSI,MTR))
p <- md.pairs(brain.adjAgeSex%>%select(Volume,T1wSI,MTR))
p$rm/(p$rm+p$rr)

ind = !is.na(brain.adjAgeSex%>%pull(Volume));print(sum(ind))
ind = ind | !is.na(brain.adjAgeSex%>% pull(T1wSI));print(sum(ind))
ind = ind | !is.na(brain.adjAgeSex%>% pull(MTR));print(sum(ind)) #n=937: at least one brain (37 with no brain phenotypes among the ones with metabolites)
brain.adjAgeSex = brain.adjAgeSex %>% filter(ind)
rm(ind)

# 1. Generate the imputations
# m=50 creates five different imputed datasets
brain.adjAgeSex = data.frame(brain.adjAgeSex)
rownames(brain.adjAgeSex) = brain.adjAgeSex$uniqueID

## female
brain.adjAgeSex.F = brain.adjAgeSex %>% 
  left_join(covdata,join_by(uniqueID)) %>%
  filter(Sex=="F") %>% select(uniqueID, Volume, T1wSI, MTR)#493
brain.adjAgeSex.F = data.frame(brain.adjAgeSex.F)
rownames(brain.adjAgeSex.F) = brain.adjAgeSex.F$uniqueID

imputed_data50 <- mice(brain.adjAgeSex.F %>% select(-1), m=50, method='pmm', seed=500)
fit50 <- with(imputed_data50, lm(Volume ~ T1wSI+MTR))
pooled_results50 <- pool(fit50)
summary(pooled_results50, conf.int = TRUE, custom.stats = c("fmi", "lambda"))
pooled_results50$pooled[, c("term", "fmi", "lambda")]
pool.r.squared(fit50)
stripplot(imputed_data50, pch = 20, cex = 1.2)

# 1. Stack all 50 versions
all_imputations <- complete(imputed_data50, action = "long")

# 2. Group by the original row ID and calculate the mean for each column
averaged_data.F <- all_imputations %>%
  group_by(.id) %>%
  summarise(across(everything(), mean)) %>%
  select(-.imp) # Remove the imputation index column

# Now 'averaged_data' is a single dataframe with no missing values
averaged_data.F %>% describe()
averaged_data.F = averaged_data.F %>% rename(uniqueID=.id)
averaged_data.F %>% describe

## male
brain.adjAgeSex.M = brain.adjAgeSex %>% left_join(covdata,join_by(uniqueID)) %>%
  filter(Sex=="M") %>% select(uniqueID, Volume, T1wSI, MTR)#493
brain.adjAgeSex.M = data.frame(brain.adjAgeSex.M)
rownames(brain.adjAgeSex.M) = brain.adjAgeSex.M$uniqueID

imputed_data50 <- mice(brain.adjAgeSex.M %>% select(-1), m=50, method='pmm', seed=500)
fit50 <- with(imputed_data50, lm(Volume ~ T1wSI+MTR))
pooled_results50 <- pool(fit50)
summary(pooled_results50, conf.int = TRUE, custom.stats = c("fmi", "lambda"))
pooled_results50$pooled[, c("term", "fmi", "lambda")]
pool.r.squared(fit50)
stripplot(imputed_data50, pch = 20, cex = 1.2)

# 1. Stack all 50 versions
all_imputations <- complete(imputed_data50, action = "long")

# 2. Group by the original row ID and calculate the mean for each column
averaged_data.M <- all_imputations %>%
  group_by(.id) %>%
  summarise(across(everything(), mean)) %>%
  select(-.imp) # Remove the imputation index column

# Now 'averaged_data' is a single dataframe with no missing values
averaged_data.M %>% describe()
averaged_data.M = averaged_data.M %>% rename(uniqueID=.id)
averaged_data = rbind(averaged_data.F,averaged_data.M)
averaged_data %>%head
averaged_data = data.frame(averaged_data)
rownames(averaged_data) = averaged_data %>% pull(1)

## PCA of brain phenotypes ----
require(FactoMineR)#
na.omit(brain.adjAgeSex %>% select(uniqueID, Volume, T1wSI, MTR)) %>% describe#n=680
brainPCA_naomit = PCA(na.omit(brain.adjAgeSex %>% select(uniqueID, Volume, T1wSI, MTR))%>% 
                        select(-1));plot(brainPCA_naomit,choix = 'var')
brainPCA_imputed = PCA(averaged_data%>% select(-1))
brainPCA_naomit$var$coord
graphics.off()
#             Dim.1      Dim.2       Dim.3
# Volume -0.8664433  0.2780914  0.41465795
# T1wSI   0.8908548 -0.1483761  0.42937418
# MTR     0.4098620  0.9103846 -0.05668365
brainPCA_imputed$var$coord
#             Dim.1      Dim.2       Dim.3
# Volume -0.8771081  0.2919792  0.38135216
# T1wSI   0.9032718 -0.1638853  0.39653720
# MTR     0.4544694  0.8892353 -0.05213502
describe(averaged_data)
describe(brain.adjAgeSex)
describe(na.omit(brain.adjAgeSex))
op = par(mfcol=c(3,3))
#
hist(brain.adjAgeSex[["Volume"]])
hist(averaged_data[["Volume"]])
hist(na.omit(brain.adjAgeSex)[["Volume"]])
#
hist(brain.adjAgeSex[['T1wSI']])
hist(averaged_data[['T1wSI']])
hist(na.omit(brain.adjAgeSex)[['T1wSI']])
#
hist(brain.adjAgeSex[['MTR']])
hist(averaged_data[['MTR']])
hist(na.omit(brain.adjAgeSex)[['MTR']])
par(op)

# change the signs of PC1 
head(brainPCA_imputed$ind$coord)

# composite variable ----
WM.PC1 = brainPCA_imputed$ind$coord[,1,drop=F]
colnames(WM.PC1) = 'WM.PC1'
WM.PC1 = data.frame(uniqueID = rownames(WM.PC1),-WM.PC1)#change sign
WM.PC1 %>% dim #n=937
metabo.PC = pca.adj.metabo$ind$coord[,1:3]
pca.adj.metabo$eig %>% head
colnames(metabo.PC) = str_replace(colnames(metabo.PC),"Dim.","metaboPC")
metabo.PC = data.frame(uniqueID = rownames(metabo.PC),metabo.PC)
metabo.PC %>% head
#        eigenvalue percentage of variance cumulative percentage of variance
# comp 1  281.84800              53.685334                          53.68533
# comp 2   69.56195              13.249896                          66.93523
# comp 3   25.93634               4.940256                          71.87549

PCall = WM.PC1 %>% filter(uniqueID %in% metabo.PC$uniqueID) %>% 
  left_join(metabo.PC) %>%
  left_join(adj.metabo %>% 
              select(uniqueID,`TAG52:4FA18:0`,`Gly`,`PC(16:0/20:4)`,`TAG54:6FA18:2`))
PCall %>% left_join(averaged_data)
head(PCall)

ggcorrplot::ggcorrplot(
  cor(PCall %>% left_join(averaged_data) %>% select(-1), use='p'),
  hc.order = T, lab = T)
pheatmap(cor(PCall %>% left_join(averaged_data) %>% select(-1), use='p'))

pa=PCall %>% left_join(averaged_data) %>% ggplot(aes(x=inormal(metaboPC1),y=inormal(WM.PC1))) + 
  geom_point() + geom_smooth(method="gam") + 
  stat_cor()
pb=PCall %>% left_join(averaged_data) %>% ggplot(aes(x=inormal(metaboPC3),y=inormal(WM.PC1))) + 
  geom_point() + geom_smooth(method="gam") + 
  stat_cor()
pc=PCall %>% left_join(averaged_data) %>% ggplot(aes(x=metaboPC1,y=WM.PC1)) + 
  geom_point() + geom_smooth(method="gam") + 
  stat_cor()
pd=PCall %>% left_join(averaged_data) %>% ggplot(aes(x=metaboPC3,y=WM.PC1)) + 
  geom_point() + geom_smooth(method="gam") + 
  stat_cor()
pa + pb + pc + pd + plot_layout(ncol=2)

#composite variable
comp1 = PCA(PCall %>% mutate(WM.PC1=inormal(WM.PC1),metaboPC1=inormal(metaboPC1)) %>% 
              select(WM.PC1,metaboPC1))#500*450
comp1$var$coord

comp2 = PCA(PCall %>% mutate(WM.PC1=inormal(WM.PC1),metaboPC3=inormal(metaboPC2)) %>% 
              select(WM.PC1,metaboPC2))#500*450
comp3 = PCA(PCall %>% mutate(WM.PC1=inormal(WM.PC1),metaboPC3=inormal(metaboPC3)) %>% 
              select(WM.PC1,metaboPC3))#500*450
comp3$var$coord

hist(comp1$ind$coord[,1])
hist(comp2$ind$coord[,1])
hist(comp3$ind$coord[,1])

# [figure4]_boxplots ----
## boxplot: WM.PC1
brainPCA_naomit$var$coord
brainPCA_naomit$eig#PC1 explains 57%
brainPCA_imputed$eig#PC1 explains 60%
plot(brainPCA_imputed, choix = 'var', axes = c(1,2) )
dt = rownames_to_column(data.frame(brainPCA_imputed$var$coord),var="WM measure") %>%
  mutate(WM.PC1_direction_adjusted = (-1)*Dim.1) %>%
  select(1,5,2,3,4)
#550*150
graphics.off();
gridExtra::grid.table(dt %>% mutate(across(-1,~round(., digits=2))))
dt_brain_pca = rownames_to_column(data.frame(brainPCA_imputed$eig),var="Dim")
graphics.off();
gridExtra::grid.table(dt_brain_pca %>% select(1,3) %>% mutate(across(-1,~round(., digits=2))))

base.size=10
box_WM.PC1 =
  brainPCA_imputed$var$coord %>% data.frame %>% 
  rownames_to_column(var = 'pheno') %>% 
  select(pheno, WM.PC1 = Dim.1) %>%
  mutate(pheno = factor(pheno,levels=c('Volume','T1wSI','MTR')), 
         WM.PC1 = -(1)*WM.PC1) %>% #flip
  ggplot(aes(x=`pheno`, y=WM.PC1)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = base.size) + 
  geom_hline(yintercept = 0, color = 'grey75',linewidth=0.5,alpha=0.5) +
  xlab(NULL)+ylab(NULL) + #+ ylab('loading') +   ggtitle('MetabolitePC1')
  scale_y_continuous( minor_breaks = seq(-1, 1, 0.1) ) + 
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
        axis.text = element_text(size=10)) + ggtitle(NULL)

#boxplots: metaboPC-loadings ----
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
pca.loadings = fread('assoc_WMphenos_with_lipidPC_loadings_525metabos_2026-03-24.txt')
pca.loadings = pca.loadings%>% select(metabolite_id,class,`lipid class`,metaboPC1,metaboPC3)
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
pca.loadings %>% filter(!is.na(metaboPC3)) %>%
  select(class,`lipid class`) %>% table(useNA='a') %>% sum#525

top_loaded = pca.loadings %>%  filter(!is.na(metaboPC3)) %>%
  filter(class=="TAG") %>% arrange(metaboPC3) %>% slice(1:2) %>%
  bind_rows(pca.loadings %>%  filter(!is.na(metaboPC3)) %>%
              filter(class=="TAG") %>% arrange(-metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="DAG") %>% arrange(metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="DAG") %>% arrange(-metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="CE") %>% arrange(-metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="PC") %>% arrange(-metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="PE") %>% arrange(-metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="LPC") %>% arrange(metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="LPE") %>% arrange(metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="SM") %>% arrange(-metaboPC3) %>% slice(1:2)) %>%
  bind_rows(pca.loadings %>% filter(!is.na(metaboPC3)) %>%
              filter(class=="NonL") %>% arrange(metaboPC3) %>% slice(1:2)) %>%
  select(-metaboPC1)
top_loaded

metaboPC3 = pca.loadings %>% filter(!is.na(metaboPC3)) %>%
  mutate(class = factor(class,levels = class.levels)) %>%
  mutate(`lipid class` = factor(`lipid class`,levels = class2.levels)) %>%
  # mutate(`lipid class` = ifelse(`lipid class`=="NonL","Non-lipid metabolites",`lipid class`)) %>%
  # res_WM1 %>% 
  # right_join(pca.loadings, join_by(metabolite_id)) %>% 
  arrange(class,`lipid class`) %>% 
  ggplot(aes(x=class, y=metaboPC3, color=`lipid class`)) +
  theme_bw(base_size = base.size) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=`lipid class`), width=0.15, size=1, shape=2, alpha=0.5)+
  # geom_hline(yintercept = c(-0.4,0.4), color="red",linetype=3) + 
  scale_color_manual(values=cols,name="Lipid classes",
                     labels = class2.labels,
                     guide = guide_legend(override.aes = list(color=cols) ))+ # Edit legend title and label
  geom_hline(yintercept = 0, color = 'grey75',linewidth=0.5,alpha=0.5) +
  xlab(NULL)+ylab(NULL) + #+ ylab('loading') +   ggtitle('MetabolitePC1')
  scale_y_continuous( minor_breaks = seq(-1, 1, 0.1) ) + 
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

# metaboPC3 = metaboPC3 + geom_text_repel(data = top_loaded ,aes(label = metabolite_id,x = class),size = 2,color='black')
p = metaboPC3+ ylab('PC loading')+
  box_WM.PC1 + ylab('PC loading')+
  plot_layout(ncol=2, widths = c(10,3)) &
  theme(legend.position='none',
        plot.margin = margin(r = 0.5, t=0.65,unit = "cm")
        # plot.margin = margin(t=0.33, b = 0.125, l=0.125, r = 1, unit = "cm")
  )
p
ggsave("PCloading_WMPC1_metPC_525metabos_2026-03-26.png",
       plot=p,height=10,width=20,unit='cm',dpi = 300)

to_write = res_WM1 %>% 
  mutate(class = factor(class,levels = class.levels)) %>%
  left_join(pca.loadings %>% select(-class,-`lipid class`), 
            join_by(metabolite_id)) %>%
  arrange(class,`lipid class`)
# write_tsv(to_write, 'assoc_WMphenos_with_lipidPC_loadings_525metabos_2026-03-24.txt')
write_tsv(to_write, 'assoc_WMphenos_with_lipidPC_loadings_525metabos_2026-10-19.txt')

# create gwas inputs ----
inputdir="GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos_2026-10-19"
dir.create(inputdir)

## create composite variables:
WM.PC1 = brainPCA_imputed$ind$coord[,1,drop=F]
WM.PC1 = data.frame(uniqueID = rownames(WM.PC1),-WM.PC1) %>% #flip the sign
  rename(WM.PC1=Dim.1)
WM.PC1 %>% describe#n=937

PCmet = pca.adj.metabo$ind$coord[,c(1:3)]
PCmet = data.frame(uniqueID = rownames(PCmet),PCmet) %>% 
  rename(PC1met = Dim.1,PC2met = Dim.2,PC3met = Dim.3)#n=974
PCmet %>% describe

compVar_input = WM.PC1 %>% inner_join(PCmet) %>% data.frame
compVar_input %>% head
rownames(compVar_input) = compVar_input$uniqueID
compVar_input = as.matrix(compVar_input %>% select(-1))
compVar_pca1 = PCA(compVar_input[,c(1,2)] )
compVar_pca3 = PCA(compVar_input[,c(1,3)])
## 
pheno = data.frame(uniqueID = rownames(compVar_pca1$ind$coord[,1,drop=F]),
                   compVar1 = compVar_pca1$ind$coord[,1]) %>%
  left_join(
    data.frame(uniqueID = rownames(compVar_pca3$ind$coord[,1,drop=F]),
               compVar3 = compVar_pca3$ind$coord[,1])
  ) %>%
  full_join(brain.adjAgeSex %>% select(uniqueID,Volume,T1wSI,MTR)) %>%
  full_join(adj.metabo %>% 
              select(uniqueID,`TAG52:4FA18:0`,`Gly`,`PC(16:0/20:4)`,`TAG54:6FA18:2`))

pheno %>% describe
pheno %>% ggplot(aes(x=compVar1,y=compVar3)) + geom_point()
# remove one outlier 
pheno$compVar3[which.min(pheno$compVar3)] <- NA
pheno %>% ggplot(aes(x=compVar1,y=compVar3)) + geom_point()

names(averaged_data)[-1]=paste(names(averaged_data)[-1],"imputed",sep=".")
ggcorrplot::ggcorrplot(pheno %>% 
                         left_join(averaged_data)  %>%
                         left_join(PCmet) %>% 
                         left_join(WM.PC1) %>% 
                         select(-1) %>% cor(use='p'),
                       hc.order = T)
phenodata = pheno %>% 
  full_join(averaged_data)  %>%
  left_join(PCmet) %>% 
  left_join(WM.PC1) %>% 
  relocate(uniqueID)

phenodata %>% describe
phenodata = data.frame(phenodata)
ind = !names(phenodata) %in% c('FID','IID','uniqueID')
sum(ind)
for (i in 2:ncol(phenodata)){
  yi = names(phenodata)[i]
  print(yi)
  hist(phenodata[[yi]],main=yi)
  phenodata[[yi]] = inormal(phenodata[[yi]])
  hist(phenodata[[yi]],main=yi)
}

hist(phenodata$T1wSI)
## correlation between WM.PC1 vs. PC.Met 
cor.test(phenodata$WM.PC1,phenodata$PC1met)#p=7.1e-13 (R=-0.23)
cor.test(phenodata$WM.PC1,phenodata$PC2met)#p=0.09 (R=-0.06)
cor.test(phenodata$WM.PC1,phenodata$PC3met)#p=8.0e-07 (R=-0.16)
cor.test(phenodata$WM.PC1,phenodata$PC1met,method = 's')#p=7.1e-13 (R=-0.23)
cor.test(phenodata$WM.PC1,phenodata$PC2met,method = 's')#p=0.09 (R=-0.06)
cor.test(phenodata$WM.PC1,phenodata$PC3met,method = 's')#p=8.0e-07 (R=-0.16)

phenodata%>% #select(WM.PC1,PC2met) %>% 
  ggplot(aes(x=inormal(WM.PC1),y=inormal(PC2met))) +
  geom_point() + 
  geom_smooth() + stat_cor() #non-linear relationship with PC2met
pca.adj.metabo$var$coord %>% data.frame %>% select(1:3) %>% arrange(-Dim.2) %>% tail(20)

# merge datasets 
geno_pcs = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/info_used_for_gwas.txt/sys_ado1_covariates_to_be_reordered.txt")

merge_data = phenodata %>%
  left_join(covdata %>% select(-age) %>% 
              mutate(IID=str_split(uniqueID,"_",simplify = T)[,2]) %>%
              mutate(IID=as.integer(IID))) %>%
  left_join(geno_pcs) %>%
  mutate(fam_id = str_split(uniqueID,"_",simplify = T)[,1],
         sub_id = str_split(uniqueID,"_",simplify = T)[,2])
merge_data %>% describe

## pheno
pheno = merge_data %>% 
  select(FID=IID,IID=IID,all_of(names(phenodata)[-1])) %>%
  mutate(compVar1 = inormal(compVar1),
         compVar3 = inormal(compVar3))
head(pheno)
write_tsv(pheno,file.path(inputdir,"pheno.input"))

## cov
cov = merge_data %>% select(FID=IID,IID=IID,age.c,age.c2,Sex, starts_with("geno"))
cov %>% head
cov = cov %>% mutate(Sex = as.numeric(factor(Sex, levels = c("F","M")))) #f=1, m=2
write_tsv(cov,file.path(inputdir,"cov.input"))

## fam
fam = merge_data %>% select(FID,IID)
write_tsv(fam,file.path(inputdir,"fam.input"))

## phenonames
cat(names(pheno%>% select(-IID,-FID)),sep='\n',
    file = file.path(inputdir,'phenoNames.input'))

## move the input files to SciNet: '/scratch/jshinb/run_gwas_regenie_20260324_zp'
load('results/adj.metabo.Rdata')
cov_input_regenie = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos/cov.input")
snp_dosages = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_dosages.csv")
snp_dosages = snp_dosages %>% left_join(
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/compVar_WM.PC1_Met.PC3_20260421_dosages.csv"))
snp_dosages_info = rbind(
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_snps_info.csv") %>% select(-1),
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/compVar_WM.PC1_Met.PC3_20260421_snps_info.csv")
  ) %>% arrange(row.id)

compvar = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos/pheno.input") %>% 
  left_join(
    snp_dosages %>%
      rename(IID = sub_id) %>% 
      select(IID,rs174566,rs174546) %>%
      rename(rs174566_G=rs174566,rs174546_T=rs174546)
  ) %>% 
  left_join(
    adj.metabo %>% select(uniqueID,`TAG50:3FA18:2`,`TAG50:4FA20:4`,
                          `TAG54:7FA16:1`,`LPC(18:2)`) %>% 
      left_join(covdata %>% select(uniqueID,fam_id,BMI)) %>%
              mutate(IID = as.integer(str_split(uniqueID,"_",simplify = T)[,2]))) 
  
# compvar = compvar %>% 
#   left_join(
#     PC2.Met.scores %>% mutate(IID = as.integer(str_split(uniqueID,"_",simplify = T)[,2]))
#     )
compvar %>% select(rs174546_T,
                   rs174566_G,
                   compVar1,compVar2,PC1.WM=WM.PC1.rnt,
                   PC1.Met = Met.PC1, 
                   # PC2.Met = PC2.Met, 
                   PC3.Met = Met.PC3, 
                   `TAG50:3FA18:2`,
                   `TAG50:4FA20:4`,
                   `TAG54:7FA16:1`,
                   `LPC(18:2)`,
                   Volume.WM.adjAgeSexICV= Volume.WM.imp.adjAgeSexICV,#Volume.WM.adjAgeSexICV.rnt,
                   T1wSI.WM.adjAgeSex=T1wSI.WM.imp.adjAgeSex,#T1wSI.WM.adjAgeSex.rnt,
                   MTR.WM.adjAgeSex=MTR.WM.imp.adjAgeSex,BMI#MTR.WM.adjAgeSex.rnt,
) %>%
  mutate(compVar2 = -compVar2) %>% #flipped
  cor(use='p') %>% ggcorrplot(lab = T)

compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age)) %>% 
  filter(!is.na(rs174546_T)) %>%
  # mutate(rs174546_T = factor(rs174546_T)) %>%
  
 # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
ggplot(aes(x=rs174546_T,y=inormal(`TAG54:7FA16:1`))) +
  # geom_boxplot() + 
  geom_point() + geom_smooth(method="lm") + stat_cor()

compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age,fam_id)) %>% 
  filter(!is.na(rs174546_T)) %>%
  mutate(rs174546_T = factor(rs174546_T)) %>%
  
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(`TAG54:7FA16:1`))) +
  geom_boxplot()

compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age,fam_id)) %>% 
  filter(!is.na(rs174546_T)) %>%
  mutate(rs174546_T = factor(rs174546_T)) %>%
  
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(`LPC(18:2)`))) +
  geom_boxplot()
graphics.off()

analdat = compvar%>% left_join(cov_input_regenie) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age,fam_id))
names(analdat)

fit = lmer(inormal(`TAG54:7FA16:1`)~as.factor(rs174546_T)+genoPC1+genoPC2+genoPC3+genoPC4+(1|fam_id),
           data=analdat)

fit.recessive = lmer(inormal(`LPC(18:2)`)~as.factor(rs174546_T==2)+genoPC1+genoPC2+genoPC3+genoPC4+(1|fam_id),
           data=analdat)
fit.additive = lmer(inormal(`LPC(18:2)`)~rs174546_T+genoPC1+genoPC2+genoPC3+genoPC4+(1|fam_id),
           data=analdat)
summary(fit.additive)
summary(fit.recessive)#p=0.034

analdat %>%
  filter(!is.na(rs174546_T)) %>%
  mutate(rs174546_T = factor(rs174546_T)) %>%
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(MTR.WM.adjAgeSex))) +
  geom_boxplot()

fit.additive = lm(inormal(T1wSI.WM.imp.adjAgeSex.rnt)~rs174546_T,
                     data=analdat)
summary(fit.additive)
fit.recessive = lm(inormal(MTR.WM.imp.adjAgeSex)~as.factor(rs174546_T==2),
                     data=analdat)
summary(fit.recessive)#p=0.16

fit.recessive = lm(inormal(WM.PC1)~as.factor(rs174546_T==2),
                   data=analdat)
summary(fit.recessive)#p=0.114 (additive); p=0.234 (recessive)

# create gwas inputs: adjBMI (2026-10-16) ----
## previous input datasets ----
pheno_2026mar26_file = "GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos/pheno.input"
pheno_2026mar26 = fread(pheno_2026mar26_file)

inputdir="GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos_adjBMI"
dir.create(inputdir)

## adjust for BMI ----
pheno_2026mar26$compVar2[which.min(pheno_2026mar26$compVar2)] <- NA
pheno_2026mar26$compVar2[which.max(pheno_2026mar26$compVar2)] <- NA
hist(pheno_2026mar26$compVar2)

data.fit.compVar3_BMI = pheno_2026mar26 %>% select(IID,contains('compVar')) %>%
  left_join(
    covdata %>% select(uniqueID,BMI,fam_id,Sex,age,age.c,age.c2) %>%
      mutate(IID = str_split(uniqueID,"_",simplify = T)[,2]) %>%
      mutate(IID = as.integer(IID)) %>%
      mutate(bmi2 = (BMI/2)^2,bmi3 = (BMI/2)^3))
data.fit.compVar3_BMI.names = names(data.fit.compVar3_BMI)
data.fit.compVar3_BMI = data.frame(data.fit.compVar3_BMI)
rownames(data.fit.compVar3_BMI) = data.fit.compVar3_BMI$uniqueID
names(data.fit.compVar3_BMI) = data.fit.compVar3_BMI.names

### linear BMI
fit.compVar3_BMI.linear = lmerTest::lmer(
  compVar2~(BMI+age.c+age.c2)*Sex + (1|fam_id), 
  data=data.fit.compVar3_BMI,
  na.action=na.exclude)
# add adjBMI.linear
data.fit.compVar3_BMI = data.fit.compVar3_BMI %>% 
  mutate(compVar3.adjBMI.linear = resid(fit.compVar3_BMI.linear) )

### cubic BMI (too complicated??)
fit.compVar3_BMI = lmerTest::lmer(
  compVar2~(age.c+age.c2 + BMI+bmi2+bmi3)*Sex + (1|fam_id), 
  data=data.fit.compVar3_BMI,
  na.action=na.exclude)
plot(fit.compVar3_BMI)

# what to do: non-constant variation in residuals (fitted vs, residulas)
qqnorm(resid(fit.compVar3_BMI),main="reisd(fit.compVar3_BMI)");qqline(resid(fit.compVar3_BMI))
m = fit.compVar3_BMI
qqnorm(ranef(m)$fam_id[[1]],main='(1|fam_id) effect');qqline(ranef(m)$fam_id[[1]])

# add adjBMI.cubic 
data.fit.compVar3_BMI = data.fit.compVar3_BMI %>% 
  mutate(compVar3.adjBMI = resid(fit.compVar3_BMI) )

## visualize
#*[cubic bmi-fit is OK]
data.fit.compVar3_BMI%>% left_join(covdata %>% select(uniqueID,BMI,age)) %>% 
  mutate(IID = str_split(uniqueID,"_",simplify = T)[,2]) %>%
  mutate(IID = as.integer(IID))%>%
  left_join(cov %>% select(IID,contains("geno"))) %>%
  ggplot(aes(x=BMI,y=compVar3.adjBMI,color=Sex)) + #geom_point()+ 
  geom_smooth(method="gam") #+ geom_smooth(method = 'lm',col='red',fill='red',alpha=0.2)

#*[linear bmi-fit is not enough]
data.fit.compVar3_BMI%>% left_join(covdata %>% select(uniqueID,BMI,age)) %>% 
  ggplot(aes(x=BMI,y=compVar3.adjBMI.linear)) + #geom_point()+ 
  geom_smooth(method="gam") + geom_smooth(method = 'lm',col='red',fill='red',alpha=0.2)

#*[however, residuals are associated with fitted values !!! non-constant residuals]
data.frame(fitted.fit=fitted(fit.compVar3_BMI),
           resid.fit=resid(fit.compVar3_BMI)) %>% 
  ggplot(aes(x=fitted.fit,y=resid.fit)) + geom_smooth(method="gam")

## remove one outlier 
#pheno$compVar3[which.min(pheno$compVar3)] <- NA

phenodata = pheno_2026mar26 %>% 
  left_join(data.fit.compVar3_BMI %>% 
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>% 
             
              mutate(FID2 = str_split(uniqueID,"_",simplify = T)[,1])%>% 
              mutate(FID2 = str_remove(FID2,"FAM"))%>%
              mutate(FID2 = as.integer(FID2)) %>%
              select(IID,FID2,compVar3.adjBMI,uniqueID))
phenodata %>% head#938 ->935

# merge datasets ----
geno_pcs = fread("~/Library/CloudStorage/OneDrive-Personal/SYS_database/info_used_for_gwas.txt/sys_ado1_covariates_to_be_reordered.txt")

merge_data = phenodata %>%
  left_join(covdata %>% select(-age) %>% 
              mutate(IID=str_split(uniqueID,"_",simplify = T)[,2]) %>%
              mutate(IID=as.integer(IID))) %>%
  left_join(geno_pcs %>% select(IID,contains("genoPC")),join_by(IID)) %>%
  mutate(fam_id = str_split(uniqueID,"_",simplify = T)[,1],
         sub_id = str_split(uniqueID,"_",simplify = T)[,2])
merge_data %>% describe %>% select(1:4)#1010

## pheno
## fam
fam = merge_data %>% select(FID=FID2,IID)
# write_tsv(fam,file.path(inputdir,"fam.input"))

## pheno2 and cov2 are the ones after removing outliers
pheno = merge_data %>% 
  select(FID=IID,IID=IID,compVar3.adjBMI,compVar2.rnt ) %>%
  mutate(compVar3.adjBMI.rnt = inormal(compVar3.adjBMI))
head(pheno)
pheno %>% describe
# write_tsv(pheno,file.path(inputdir,"pheno2.input"))

## cov
cov = merge_data %>% select(FID=IID,IID=IID,age.c,age.c2,BMI,Sex, starts_with("geno")) %>%
  mutate(bmi2=(BMI/2)^2,bim3=(BMI/2)^3)
cov %>% head
cov %>% describe
cov = cov %>% mutate(Sex = as.numeric(factor(Sex, levels = c("F","M")))) #f=1, m=2
# write_tsv(cov,file.path(inputdir,"cov2.input"))

## phenonames
# cat(names(pheno%>% select(-IID,-FID)),sep='\n',
#     file = file.path(inputdir,'phenoNames.input'))
ggcorrplot(cor(pheno%>% select(-IID,-FID),use='p'),lab=T)

## move the input files to SciNet: '/scratch/jshinb/run_gwas_regenie_20260324_zp'##

# assoc with the top SNP ----
load('results/adj.metabo.Rdata')
compvar = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos/pheno.input")
covdata = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/GWAS_input_comp_vars_WM.PC1_metPC1_metPC3_525metabos/cov.input")
snp_dosages = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_dosages.csv")
snp_dosages = snp_dosages %>% left_join(
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/compVar_WM.PC1_Met.PC3_20260421_dosages.csv"))
snp_dosages_info = rbind(
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/gwas_results/PC1_mean_imputed/NormWM_METABOLON_ADO1_2025-02-12_snps_info.csv") %>% select(-1),
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/compVar_WM.PC1_Met.PC3_20260421_snps_info.csv")
) %>% arrange(row.id)

compvar = compvar %>% 
  left_join(
    snp_dosages %>%
      rename(IID = sub_id) %>% 
      select(IID,rs174566,rs174546) %>%
      rename(rs174566_G=rs174566,rs174546_T=rs174546)
  ) %>% 
  left_join(
    adj.metabo %>% select(uniqueID,`TAG50:3FA18:2`,`TAG50:4FA20:4`,
                          `TAG54:7FA16:1`,`LPC(18:2)`) %>% 
      left_join(covdata_all %>% select(uniqueID,fam_id,BMI)) %>%
      mutate(IID = as.integer(str_split(uniqueID,"_",simplify = T)[,2]))) 

# compvar = compvar %>% 
#   left_join(
#     PC2.Met.scores %>% mutate(IID = as.integer(str_split(uniqueID,"_",simplify = T)[,2]))
#     )
compvar %>% select(rs174546_T,
                   rs174566_G,
                   compVar1,compVar2,PC1.WM=WM.PC1.rnt,
                   PC1.Met = Met.PC1, 
                   # PC2.Met = PC2.Met, 
                   PC3.Met = Met.PC3, 
                   `TAG50:3FA18:2`,
                   `TAG50:4FA20:4`,
                   `TAG54:7FA16:1`,
                   `LPC(18:2)`,
                   Volume.WM.adjAgeSexICV= Volume.WM.imp.adjAgeSexICV,#Volume.WM.adjAgeSexICV.rnt,
                   T1wSI.WM.adjAgeSex=T1wSI.WM.imp.adjAgeSex,#T1wSI.WM.adjAgeSex.rnt,
                   MTR.WM.adjAgeSex=MTR.WM.imp.adjAgeSex,BMI#MTR.WM.adjAgeSex.rnt,
) %>%
  mutate(compVar2 = -compVar2) %>% #flipped
  cor(use='p') %>% ggcorrplot(lab = T)

compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age)) %>% 
  filter(!is.na(rs174546_T)) %>%
  # mutate(rs174546_T = factor(rs174546_T)) %>%
  
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(`TAG54:7FA16:1`))) +
  # geom_boxplot() + 
  geom_point() + geom_smooth(method="lm") + stat_cor()

compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age,fam_id)) %>% 
  filter(!is.na(rs174546_T)) %>%
  mutate(rs174546_T = factor(rs174546_T)) %>%
  
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(`TAG54:7FA16:1`))) +
  geom_boxplot()

compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age,fam_id)) %>% 
  filter(!is.na(rs174546_T)) %>%
  mutate(rs174546_T = factor(rs174546_T)) %>%
  
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(`LPC(18:2)`))) +
  geom_boxplot()

dev.off()

names(d_merge)
analdat = compvar%>% left_join(covdata) %>%
  left_join((fread(file.path(wd,"data","covdata_ados.txt"))) %>%
              mutate(IID=as.integer(str_split(uniqueID,"_",simplify = T)[,2])) %>%
              select(IID,age,fam_id))
names(analdat)
fit = lmer(inormal(`TAG54:7FA16:1`)~as.factor(rs174546_T)+genoPC1+genoPC2+genoPC3+genoPC4+(1|fam_id),
           data=analdat)

fit.recessive = lmer(inormal(`LPC(18:2)`)~as.factor(rs174546_T==2)+genoPC1+genoPC2+genoPC3+genoPC4+(1|fam_id),
                     data=analdat)
fit.additive = lmer(inormal(`LPC(18:2)`)~rs174546_T+genoPC1+genoPC2+genoPC3+genoPC4+(1|fam_id),
                    data=analdat)
summary(fit.additive)
summary(fit.recessive)

analdat %>%
  filter(!is.na(rs174546_T)) %>%
  mutate(rs174546_T = factor(rs174546_T)) %>%
  # ggplot(aes(x=rs174546_T,y=inormal(`TAG50:4FA20:4`))) +
  ggplot(aes(x=rs174546_T,y=inormal(MTR.WM.adjAgeSex))) +
  geom_boxplot()

fit.additive = lm(inormal(T1wSI.WM.imp.adjAgeSex.rnt)~rs174546_T,
                  data=analdat)
summary(fit.additive)
fit.recessive = lm(inormal(MTR.WM.imp.adjAgeSex)~as.factor(rs174546_T==2),
                   data=analdat)
summary(fit.recessive)


fit.recessive = lm(inormal(WM.PC1)~rs174546_T,
                   data=analdat)

fit.recessive = lm(inormal(WM.PC1)~rs174546_T,
                   data=analdat)
summary(fit.recessive)
