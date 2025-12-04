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
TAG_fa = TAG_fa %>% filter(!is.na(PC1))
TAG_fa %>% describe()
TAG_fa[['ARA.or.DHA']] = ""
TAG_fa$ARA.or.DHA[TAG_fa$ARA=="ARA-present"] <- "ARA"
TAG_fa$ARA.or.DHA[TAG_fa$DHA=="DHA-present"] <- "DHA"
table(TAG_fa$ARA.or.DHA)

my_colors <- c("ARA" = "#0072B2", "DHA" = "#D55E00") # Blue and Orange
# my_colors <- c("ARA" = "#702963", "DHA" = "#296370") # Blue and Orange
p_n_cb = ggplot()  + 
  geom_point(data=TAG_fa,aes(x=tagC,y=PC1),alpha=0.3) + 
  geom_point(data=TAG_fa %>% filter(ARA.or.DHA!=""),
             aes(x=tagC,y=PC1,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_smooth(data=TAG_fa, aes(x=tagC,y=PC1),method="gam",color = "gray30") +
  xlab('carbons') + theme_bw() +
  scale_x_continuous(breaks = seq(40, 60, by = 4))
p_n_cb

p_n_db = ggplot()  + 
  geom_point(data=TAG_fa,aes(x=tag_db,y=PC1),alpha=0.3) + 
  geom_point(data=TAG_fa %>% filter(ARA.or.DHA!=""),
             aes(x=tag_db,y=PC1,color=ARA.or.DHA),
             alpha=0.7,pch=21,stroke=1.5,size=2) + 
  scale_color_manual(values=my_colors) + 
  geom_smooth(data=TAG_fa, aes(x=tag_db,y=PC1),method="gam",color = "gray30") +
  xlab('double bonds') + theme_bw() + 
  scale_x_continuous(breaks = seq(0, 12, by = 2))
p_n_db

# p_n_db = TAG_fa %>% ggplot(aes(x=tag_db,y=PC1)) + 
#   geom_point(alpha=0.3) + geom_smooth(method="gam") +
#   xlab('double bonds')
sample_sizes_ARA <- 
  TAG_fa %>% select(ARA,PC1) %>%
  group_by(ARA) %>%
  summarise(n = n())

p_ARA = TAG_fa %>% 
  ggplot(aes(x=ARA,y=PC1,fill=ARA)) +
  scale_fill_manual(values=c("grey30","#0072B2")) + 
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  geom_jitter(width=0.2,alpha=0.3,size=1) + 
  geom_text(data = sample_sizes_ARA, 
            # aes(x = ARA, y = 1.05, label = paste0("n=", n)),
            aes(x = ARA, y = 1.05, label = paste0("(",n," TAGs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = 5) +  # Adjust text size
  xlab(NULL)+ theme_bw() + theme(legend.position="none") 
p_ARA

sample_sizes_DHA <- 
  TAG_fa %>% select(DHA,PC1) %>%
  group_by(DHA) %>%
  summarise(n = n())

p_DHA = TAG_fa %>% 
  ggplot(aes(x=DHA,y=PC1,fill=DHA)) +
  scale_fill_manual(values=c("grey30","#D55E00")) + 
  geom_boxplot(outlier.shape=NA,alpha=0.5) + 
  geom_jitter(width=0.2,alpha=0.3,size=1) + 
  geom_text(data = sample_sizes_DHA, 
            # aes(x = DHA, y = 1.05, label = paste0("n=", n)),
            aes(x = DHA, y = 1.05, label = paste0("(",n," TAGs)")),
            vjust = 1.5,  # Adjust vertical position of text
            size = 5) +  # Adjust text size
  xlab(NULL)+ theme_bw() + theme(legend.position="none") 
p_DHA

p_ARA + p_DHA +  p_n_db + p_n_cb + 
  plot_layout(axes = "collect",guides="collect")

ggsave('outputs/figS3_PC1_loadings_TAG.png',
       height=13.03*1.75*0.75,width=17.07*2**0.65,units='cm',
       dpi=300)

# part not used ----
not.used = function(){
  #********************************************#
  # 1. Load libraries
  #********************************************#
  library(tidyverse)
  library(patchwork)
  
  #********************************************#
  # 2. Load your data
  #********************************************#
  # Replace with your file
  # df <- read_csv("your_file.csv")
  
  #********************************************#
  # 3. RENAME COLUMNS (EDIT HERE)
  #********************************************#
  # Replace the right-hand side with your column names
  df <- TAG_fa %>%
    rename(
      PC1_loading     = PC1,
      ARA_present     = ARA,       # binary indicator (0/1 or TRUE/FALSE)
      DHA_present     = DHA,       # binary indicator
      double_bonds    = tag_db,
      total_carbons   = tagC
    )
  sum(table(df$dbl.cb.bonds))
  
  #********************************************#
  # 4. Create PC1 bins
  #********************************************#
  df <- df %>%
    mutate(
      PC1_bin = cut(
        PC1_loading,
        breaks = seq(0, 1, by = 0.1),
        include.lowest = TRUE,
        right = TRUE,
        labels = sprintf("%.1f–%.1f", seq(0,0.9,0.1), seq(0.1,1,0.1))
      )
    )
  
  #********************************************#
  # 5. Create category groups
  #********************************************#
  
  # Double bond groups
  df <- df %>%
    mutate(
      db_group = case_when(
        double_bonds == 0 ~ "0",
        double_bonds >=1 & double_bonds <=2 ~ "1-2",
        double_bonds >=3 & double_bonds <=6 ~ "3–6",
        double_bonds > 6 ~ ">6"
      )
    ) %>% mutate(db_group = factor(db_group, levels = c("0", "1-2", "3–6", ">6")))
  
  # Carbon groups
  df <- df %>%
    mutate(
      C_group = case_when(
        total_carbons < 50 ~ "<50",
        total_carbons >= 50 & total_carbons <= 60 ~ "50–60",
        total_carbons > 60 ~ ">60"
      ),
      C_group = factor(C_group, levels = c("<50","50–60",">60"))
    )
  
  #********************************************#
  # 6. Compute proportions by bin and category
  #********************************************#
  # ARA
  df_ara <- df %>%
    group_by(PC1_bin) %>%
    summarise(proportion = mean(ARA_present)) %>%
    mutate(category = "ARA", feature = "ARA")
  
  # DHA
  df_dha <- df %>%
    group_by(PC1_bin) %>%
    summarise(proportion = mean(DHA_present)) %>%
    mutate(category = "DHA", feature = "DHA")
  
  # Double bonds
  df_db <- df %>%
    group_by(PC1_bin, db_group) %>%
    summarise(proportion = n() / sum(n()), .groups = "drop") %>%
    mutate(feature = "Double Bonds")
  
  # Carbon groups
  df_carb <- df %>%
    group_by(PC1_bin, C_group) %>%
    summarise(proportion = n() / sum(n()), .groups = "drop") %>%
    mutate(feature = "Carbon Groups")
  
  #********************************************#
  # 7. Combine all datasets into long format
  #********************************************#
  df_all <- bind_rows(
    df_ara %>% rename(category = category),
    df_dha %>% rename(category = category),
    df_db  %>% rename(category = db_group),
    df_carb %>% rename(category = C_group)
  )
  
  #********************************************#
  # 8. Faceted bar plot (4 panels)
  #********************************************#
  pc1.loading.levels = c("0.0–0.1", "0.1–0.2", "0.2–0.3", "0.3–0.4", "0.4–0.5", 
                         "0.5–0.6", "0.6–0.7", "0.7–0.8", "0.8–0.9", "0.9–1.0")
  df_ara$PC1_bin = factor(df_ara$PC1_bin, levels = pc1.loading.levels)
  p_bar1 <- df_ara %>%ggplot(aes(x = PC1_bin, y = proportion)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "PC1 loading bin", y = "Proportion of TAGs with ARA")
  p_bar1
  tab1 = as.data.frame.matrix(table(df$PC1_bin,df$ARA_present)) %>% mutate(n_bin = `FALSE` + `TRUE`)
  tab1 = tab1 %>% mutate(x=rownames(tab1),y=(`TRUE`/n_bin)*1.05, label = paste(`TRUE`,n_bin,sep="/"))
  p_bar1 = p_bar1 + geom_text(data=tab1,aes(x=x,y=y,label=label))
  
  p_bar2 <- ggplot(df_dha, aes(x = PC1_bin, y = proportion)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "PC1 loading bin", y = "Proportion of TAGs with DHA")
  p_bar2
  tab2 = as.data.frame.matrix(table(df$PC1_bin,df$DHA_present)) %>% mutate(n_bin = `FALSE` + `TRUE`)
  tab2 = tab2 %>% mutate(x=rownames(tab2),y=(`TRUE`/n_bin)*1.05, label = paste(`TRUE`,n_bin,sep="/"))
  p_bar2 = p_bar2 + geom_text(data=tab2,aes(x=x,y=y,label=label))
  
  #
  tab3 = as.data.frame.matrix(table(df%>% select(PC1_bin,dbl.cb.bonds)))
  tab3 = tab3 %>% mutate(n_bin = apply(tab3,1,sum))
  tab3 = tab3 %>% mutate(PC1_bin = rownames(tab3))
  tab3_prop = tab3 %>% select(PC1_bin) %>% bind_cols(tab3[,1:4]/tab3$n_bin)
  tab3_prop_long  = tab3_prop %>% 
    pivot_longer(
      cols = c('0','1','2-6','>6'),
      names_to = 'dbl.cb.bond',
      values_to = "proportion"
    ) %>% mutate(
      dbl.cb.bond = factor(dbl.cb.bond, levels = c("0", "1", "2-6", ">6")))
  
  
  p_bar3 = ggplot(tab3_prop_long, aes(fill=dbl.cb.bond, y=proportion, x=PC1_bin)) + 
    geom_bar(position="fill", stat="identity")
  
  #
  tab4 = as.data.frame.matrix(table(df%>% select(PC1_bin,total.n.cb)))
  tab4 = tab4 %>% mutate(n_bin = apply(tab4,1,sum))
  tab4 = tab4 %>% mutate(PC1_bin = rownames(tab4))
  tab4_prop = tab4 %>% select(PC1_bin) %>% bind_cols(tab4[,1:3]/tab4$n_bin)
  tab4_prop_long  = tab4_prop %>% 
    pivot_longer(
      cols = c('<50','50-60','>60'),
      names_to = 'total.n.cb',
      values_to = "proportion"
    ) %>% 
    mutate(total.n.cb = factor(total.n.cb, levels = c('<50','50-60','>60')))
  
  p_bar4 = ggplot(tab4_prop_long, aes(fill=total.n.cb, y=proportion, x=PC1_bin)) + 
    geom_bar(position="fill", stat="identity")
  #
  (p_bar1 + p_bar2 + p_bar3 + p_bar4 ) + plot_layout(ncol=2,nrow=2) + theme_bw() &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #********************************************#
  # 9. Heatmap
  #********************************************#
  p_heat <- ggplot(df_all, aes(x = PC1_bin, y = category, fill = proportion)) +
    geom_tile(color = "white") +
    facet_wrap(~feature, scales = "free_y") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "PC1 loading bin", y = "Category", fill = "Proportion")
  
  #********************************************#
  # 10. Line plot option
  #********************************************#
  p_line <- ggplot(df_all, aes(x = PC1_bin, y = proportion,
                               color = category, group = category)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~feature, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "PC1 loading bin", y = "Proportion")
  
  #********************************************#
  # 11. Print all figures
  #********************************************#
  p_bar
  p_heat
  p_line
}
