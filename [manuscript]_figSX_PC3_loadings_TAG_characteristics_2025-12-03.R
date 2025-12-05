rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)

TAG_fa = fread("brain_tagc_table.txt")
head(TAG_fa)
TAG_fa = TAG_fa %>% mutate(metabolite = str_remove(metabolite,'_nmol/mL'))
head(TAG_fa)
TAG_fa = TAG_fa %>% mutate(ARA=ifelse(faC_db == '20:4',"ARA-present",'ARA-absent'),
                           DHA=ifelse(faC_db == '22:6','DHA-present','DHA-absent'),
                           dbl.cb.bonds = NA,
                           total.n.cb = NA)
dbl.cb.bonds.levels = c('0','1-2','3-6','>6')
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db==0] <- '0'
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db>=1&TAG_fa$tag_db<=2] <- '1-2'
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db>=3 & TAG_fa$tag_db<=6] <- '3-6'
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db>6] <- '>6'
TAG_fa = TAG_fa %>% 
  mutate(dbl.cb.bonds = factor(dbl.cb.bonds,levels=dbl.cb.bonds.levels))
table(TAG_fa %>% select(tag_db,dbl.cb.bonds))

table(TAG_fa$tagC)
TAG_fa$total.n.cb[TAG_fa$tagC <48] = "<48"
TAG_fa$total.n.cb[TAG_fa$tagC>=48 & TAG_fa$tag_db<=53] = "48-53"
TAG_fa$total.n.cb[TAG_fa$tagC>=54 & TAG_fa$tag_db<=58] = "54-58"
TAG_fa$total.n.cb[TAG_fa$tagC>60] = ">58"
TAG_fa = TAG_fa %>% 
  mutate(total.n.cb = factor(total.n.cb,levels=c("<48","48-53","54-58",">58")))
TAG_fa %>% select(tagC,total.n.cb) %>% describe %>% select(min,max)

table(TAG_fa$DHA)

PC_loadings_616 = readxl::read_xlsx("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/outputs/MetabolitePC_loadings_616_2025-05-26_mod_nov14.xlsx")
PC_loadings_616 = data.table(PC_loadings_616)
head(PC_loadings_616)
TAG_fa = TAG_fa %>% select(metabolite,ARA,DHA,dbl.cb.bonds,total.n.cb,tag_db,tagC) %>% 
  rename(metabolite_id = metabolite) %>%
  left_join(PC_loadings_616,join_by(metabolite_id))
head(TAG_fa)
TAG_fa = TAG_fa %>% filter(!is.na(PC3))
TAG_fa %>% describe()
TAG_fa[['ARA.or.DHA']] = ""
TAG_fa$ARA.or.DHA[TAG_fa$ARA=="ARA-present"] <- "ARA"
TAG_fa$ARA.or.DHA[TAG_fa$DHA=="DHA-present"] <- "DHA"
table(TAG_fa$ARA.or.DHA)

# generate plots ----
my_colors <- c("ARA" = "#0072B2", "DHA" = "#D55E00") # Blue and Orange
# my_colors <- c("ARA" = "#702963", "DHA" = "#296370") # Blue and Orange
p_n_cb = ggplot()  + 
  geom_point(data=TAG_fa,aes(x=tagC,y=PC3),alpha=0.3) + 
  geom_point(data=TAG_fa %>% filter(ARA.or.DHA!=""),
             aes(x=tagC,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  geom_smooth(data=TAG_fa, aes(x=tagC,y=PC3),method="gam",color = "gray30") +
  xlab('carbons') + theme_bw() +
  scale_x_continuous(breaks = seq(40, 60, by = 4))
p_n_cb

p_n_db = ggplot()  + 
  geom_point(data=TAG_fa,aes(x=tag_db,y=PC3),alpha=0.3) + 
  geom_point(data=TAG_fa %>% filter(ARA.or.DHA!=""),
             aes(x=tag_db,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) +
  geom_hline(yintercept = 0, linewidth=0.5) + 
  geom_smooth(data=TAG_fa, aes(x=tag_db,y=PC3),method="gam",color = "gray30") +
  xlab('double bonds') + theme_bw() + 
  scale_x_continuous(breaks = seq(0, 12, by = 2)) 
p_n_db

# p_n_db = TAG_fa %>% ggplot(aes(x=tag_db,y=PC3)) + 
#   geom_point(alpha=0.3) + geom_smooth(method="gam") +
#   xlab('double bonds')
sample_sizes_ARA <- 
  TAG_fa %>% select(ARA,PC3) %>%
  group_by(ARA) %>%
  summarise(n = n())

boxplot.ymax = 0.75
p_ARA = TAG_fa %>% 
  ggplot(aes(x=ARA,y=PC3,fill=ARA)) +
  scale_fill_manual(values=c("grey30","#0072B2")) + 
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  geom_jitter(width=0.2,alpha=0.3,size=1) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  geom_text(data = sample_sizes_ARA, 
            # aes(x = ARA, y = 1.05, label = paste0("n=", n)),
            aes(x = ARA, y = boxplot.ymax, label = paste0("(",n," TAGs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = 5) +  # Adjust text size
  xlab(NULL)+ theme_bw() + theme(legend.position="none") 
p_ARA

sample_sizes_DHA <- 
  TAG_fa %>% select(DHA,PC3) %>%
  group_by(DHA) %>%
  summarise(n = n())

p_DHA = TAG_fa %>% 
  ggplot(aes(x=DHA,y=PC3,fill=DHA)) +
  scale_fill_manual(values=c("grey30","#D55E00")) + 
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  geom_jitter(width=0.2,alpha=0.3,size=1) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  geom_text(data = sample_sizes_DHA, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = DHA, y = boxplot.ymax, label = paste0("(",n," TAGs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = 5) +  # Adjust text size
  xlab(NULL)+ theme_bw() + theme(legend.position="none") 
p_DHA

p_ARA + p_DHA +  p_n_db + p_n_cb + 
  plot_layout(axes = "collect",guides="collect") &
  coord_cartesian(ylim = boxplot.ymax*c(-0.85,1))

# create PNG file in 'output' folder ----
ggsave('outputs/figSX_PC3_loadings_TAG.png',
       height=13.03*1.75*0.75,width=17.07*2**0.65,units='cm',
       dpi=300)
#smooth plots ----
## CEs ----
CE_fa = PC_loadings_616 %>% filter(class=="CE") %>% 
  select(metabolite_id) %>% 
  mutate(ceC = str_remove(metabolite_id,"CE[(]")) %>%
  mutate(ceC = str_split(ceC,":",simplify = T)[,1]) %>%
  mutate(ceC = as.numeric(ceC)) %>%
  mutate(ce_db = str_remove(metabolite_id,"CE[(]")) %>%
  mutate(ce_db = str_remove(metabolite_id,"[)]")) %>%
  mutate(ce_db = as.numeric(str_split(ce_db,":",simplify = T)[,2])) %>%
  mutate(ARA = ifelse(str_detect(metabolite_id,"20:4"),"ARA-present","ARA-absent")) %>%
  mutate(DHA = ifelse(str_detect(metabolite_id,"22:6"),"DHA-present","DHA-absent"))
CE_fa[['ARA.or.DHA']] = ""
CE_fa$ARA.or.DHA[CE_fa$ARA=="ARA-present"] <- "ARA"
CE_fa$ARA.or.DHA[CE_fa$DHA=="DHA-present"] <- "DHA"
CE_fa = CE_fa %>% left_join(PC_loadings_616,join_by(metabolite_id))

p_ce_cb = 
  ggplot()  + 
  geom_point(data=CE_fa,aes(x=ceC,y=PC3),alpha=0.3) + 
  geom_point(data=CE_fa %>% filter(ARA.or.DHA!=""),
             aes(x=ceC,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  # geom_smooth(data=CE_fa, aes(x=ceC,y=PC3),
  #             method="gam",
  #             formula = y ~ s(x, bs = "cs",k = 3),
  #             color = "gray30") +
  xlab('carbons') + theme_bw() #+  scale_x_continuous(breaks = seq(40, 60, by = 4))
p_ce_cb

p_ce_db = 
  ggplot()  + 
  geom_point(data=CE_fa,aes(x=ce_db,y=PC3),alpha=0.3) + 
  geom_point(data=CE_fa %>% filter(ARA.or.DHA!=""),
             aes(x=ce_db,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  # geom_smooth(data=CE_fa, aes(x=ce_db,y=PC3),
  #             method="gam",
  #             formula = y ~ s(x, bs = "cs",k = 3),
  #             color = "gray30") +
  xlab('double bonds') + theme_bw() #+  scale_x_continuous(breaks = seq(40, 60, by = 4))
p_ce_db

## PLs ----
PL_fa = 
  PC_loadings_616 %>% filter(`lipid class`=="Phospholipids") %>% 
  select(metabolite_id) %>% 
  mutate(pl_fa = str_split(metabolite_id,"[(]",simplify = T)[,2]) %>%
  mutate(pl_fa = str_remove(pl_fa,"[)]")) %>%
  mutate(pl_fa = str_remove(pl_fa,"P")) %>% 
  mutate(pl_fa1 = str_split(pl_fa,"/",simplify = T)[,1]) %>%
  mutate(pl_fa2 = str_split(pl_fa,"/",simplify = T)[,2]) %>%
  mutate(ARA = ifelse((pl_fa1=="20:4"|pl_fa2=="20:4"),"ARA-present","ARA-absent")) %>%
  mutate(DHA = ifelse((pl_fa1=="22:6"|pl_fa2=="22:6"),"DHA-present","DHA-absent")) %>%
  mutate(plC = as.numeric(str_split(pl_fa1,":",simplify = T)[,1]) + as.numeric(str_split(pl_fa2,":",simplify = T)[,1])) %>%
  mutate(pl_db = as.numeric(str_split(pl_fa1,":",simplify = T)[,2]) + as.numeric(str_split(pl_fa2,":",simplify = T)[,2]))
PL_fa[['ARA.or.DHA']] = ""
PL_fa$ARA.or.DHA[PL_fa$ARA=="ARA-present"] <- "ARA"
PL_fa$ARA.or.DHA[PL_fa$DHA=="DHA-present"] <- "DHA"
PL_fa = PL_fa %>% left_join(PC_loadings_616,join_by(metabolite_id))
p_pl_cb = 
  ggplot()  + 
  geom_point(data=PL_fa,aes(x=plC,y=PC3),alpha=0.3) + 
  geom_point(data=PL_fa %>% filter(ARA.or.DHA!=""),
             aes(x=plC,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  # geom_smooth(data=PL_fa, aes(x=plC,y=PC3),
  #             method="gam",
  #             formula = y ~ s(x, bs = "cs",k = 3),
  #             color = "gray30") +
  xlab('carbons') + theme_bw() #+  scale_x_continuous(breaks = seq(40, 60, by = 4))
p_pl_cb

p_pl_db = 
  ggplot() + 
  geom_point(data=PL_fa,aes(x=pl_db,y=PC3),alpha=0.3) + 
  geom_point(data=PL_fa %>% filter(ARA.or.DHA!=""),
             aes(x=pl_db,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5) + 
  # geom_smooth(data=PL_fa, aes(x=pl_db,y=PC3),
  #             method="gam",
  #             formula = y ~ s(x, bs = "cs",k = 3),
  #             color = "gray30") +
  xlab('double bonds') + theme_bw() #+  scale_x_continuous(breaks = seq(40, 60, by = 4))
p_pl_db

(p_ce_cb + p_ce_db+ p_pl_cb + p_pl_db) + plot_layout(guides = 'collect',axes = 'collect',ncol=2)

## generate PNG plot ----
ggsave('outputs/figSX_PC3_CE_PL_loadings.png',
       height=13.03*1.75*0.75,width=17.07*2**0.65,units='cm',
       dpi=300)

# boxplots ----
boxplot.ymax=0.675
box.n.size = 4
## CEs: boxplot - carbons ----
sample_sizes_ce_cb <- 
  CE_fa %>% mutate(x=as.character(ceC)) %>%
  select(x,PC3) %>%
  group_by(x) %>%
  summarise(n = n())
box_ce_cb = 
  CE_fa %>% mutate(x=as.character(ceC)) %>%
  ggplot(aes(x=x,y=PC3))+
  geom_boxplot(alpha=0.3,outlier.shape = NA) +
  geom_point(alpha=0.3) + 
  geom_point(data=CE_fa %>%  mutate(x=as.character(ceC)) %>% filter(ARA.or.DHA!=""),
             aes(x=x,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  theme_bw() + xlab("carbons") +
  geom_hline(yintercept = 0, linewidth=0.5, color="grey30") +
  geom_text(data = sample_sizes_ce_cb, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = x, y =boxplot.ymax , label = paste0("(",n," CEs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = box.n.size) +  # Adjust text size
  xlab('carbons')+ theme_bw() + theme(legend.position="none") 
## CEs: boxplot - double bonds ----
sample_sizes_ce_db <- 
  CE_fa %>% mutate(x=as.character(ce_db)) %>%
  select(x,PC3) %>%
  group_by(x) %>%
  summarise(n = n())
box_ce_db = CE_fa %>% mutate(x=as.character(ce_db)) %>%
  ggplot(aes(x=x,y=PC3))+
  geom_boxplot(alpha=0.3,outlier.shape = NA) +
  geom_point(alpha=0.3) + 
  geom_point(data=CE_fa %>%  mutate(x=as.character(ce_db)) %>% filter(ARA.or.DHA!=""),
             aes(x=x,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5, color="grey30") +
  geom_text(data = sample_sizes_ce_db, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = x, y =boxplot.ymax , label = paste0("(",n," CEs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = box.n.size) +  # Adjust text size
  xlab('double bonds')+ theme_bw() + theme(legend.position="none") 
## PLs: boxplot - carbons ----
sample_sizes_pl_cb <- 
  PL_fa %>% mutate(x=as.character(plC)) %>%
  select(x,PC3) %>%
  group_by(x) %>%
  summarise(n = n())
box_pl_cb = PL_fa %>% mutate(x=as.character(plC)) %>%
  ggplot(aes(x=x,y=PC3))+
  geom_boxplot(alpha=0.3,outlier.shape = NA) +
  geom_point(alpha=0.3) + 
  geom_point(data=PL_fa %>%  mutate(x=as.character(plC)) %>% filter(ARA.or.DHA!=""),
             aes(x=x,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  theme_bw() + xlab("carbons") +
  geom_hline(yintercept = 0, linewidth=0.5, color="grey30") +
  geom_text(data = sample_sizes_pl_cb, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = x, y =boxplot.ymax , label = paste0("(",n," PLs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = box.n.size) +  # Adjust text size
  xlab('carbons')+ theme_bw() + theme(legend.position="none") 
## PLs: boxplot - double bonds ----
sample_sizes_pl_db <- 
  PL_fa %>% mutate(x=as.character(pl_db)) %>%
  select(x,PC3) %>%
  group_by(x) %>%
  summarise(n = n())
box_pl_db = 
  PL_fa %>% mutate(x=as.character(pl_db)) %>%
  ggplot(aes(x=x,y=PC3))+
  geom_boxplot(alpha=0.3,outlier.shape = NA) +
  geom_point(alpha=0.3) + 
  geom_point(data=PL_fa %>%  mutate(x=as.character(pl_db)) %>% filter(ARA.or.DHA!=""),
             aes(x=x,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5, color="grey30") +
  geom_text(data = sample_sizes_pl_db, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = x, y =boxplot.ymax , label = paste0("(",n," PLs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = box.n.size) +  # Adjust text size
  xlab('double bonds')+ theme_bw() + theme(legend.position="none") 
(box_ce_cb + box_ce_db+ box_pl_cb + box_pl_db) + plot_layout(guides = 'collect',axes = 'collect',ncol=2)
ggsave('outputs/figSX_PC3_CE_PL_loadings_boxplot.png',
       height=13.03*1.75*0.75,width=17.07*2**0.65,units='cm',
       dpi=300)
sample_sizes_pl_db <- 
  PL_fa %>% mutate(x=as.character(ARA.or.DHA)) %>%
  select(x,PC3) %>%
  group_by(x) %>%
  summarise(n = n())
PL_fa %>% 
  # mutate(x=as.character(pl_db)) %>%
  mutate(x=as.character(ARA.or.DHA)) %>%
  ggplot(aes(x=x,y=PC3))+ 
  geom_boxplot(alpha=0.3,outlier.shape = NA) +
  geom_point(alpha=0.3) + 
  geom_point(data=PL_fa %>%  
               mutate(x=ARA.or.DHA) %>%
               filter(ARA.or.DHA!=""),
             aes(x=x,y=PC3,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_hline(yintercept = 0, linewidth=0.5, color="grey30") +
  geom_text(data = sample_sizes_pl_db, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = x, y =boxplot.ymax , label = paste0("(",n," PLs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = box.n.size) +  # Adjust text size
  xlab(NULL)+ theme_bw() + theme(legend.position="none") 
