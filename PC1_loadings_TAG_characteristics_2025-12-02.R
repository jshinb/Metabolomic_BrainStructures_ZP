
TAG_fa = fread("brain_tagc_table.txt")
head(TAG_fa)
TAG_fa = TAG_fa %>% mutate(metabolite = str_remove(metabolite,'_nmol/mL'))
head(TAG_fa)
TAG_fa = TAG_fa %>% mutate(contains.ARA=ifelse(faC_db == '20:4',TRUE,FALSE),
                           contains.DHA=ifelse(faC_db == '22:6',TRUE,FALSE),
                           dbl.cb.bonds = NA,
                           total.n.cb = NA)
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db==0] <- '0'
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db==1] <- '1'
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db>=2 & TAG_fa$tag_db<=6] <- '2-6'
TAG_fa$dbl.cb.bonds[TAG_fa$tag_db>6] <- '>6'
TAG_fa = TAG_fa %>% mutate(dbl.cb.bonds = factor(dbl.cb.bonds,levels=c('0','1','2-6','>6')))
table(TAG_fa %>% select(tag_db,dbl.cb.bonds))

table(TAG_fa$tagC)
TAG_fa$total.n.cb[TAG_fa$tagC<50] = "<50"
TAG_fa$total.n.cb[TAG_fa$tagC>=50 & TAG_fa$tag_db<=50] = "50-60"
TAG_fa$total.n.cb[TAG_fa$tagC>60] = ">60"
TAG_fa = TAG_fa %>% mutate(total.n.cb = factor(total.n.cb,levels=c("<50","50-60",">60")))
TAG_fa %>% select(tagC,total.n.cb) %>% describe %>% select(min,max)

table(TAG_fa$contains.DHA)

PC_loadings_616 = readxl::read_xlsx("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/outputs/MetabolitePC_loadings_616_2025-05-26_mod_nov14.xlsx")
PC_loadings_616 = data.table(PC_loadings_616)
head(PC_loadings_616)
TAG_fa = TAG_fa %>% select(metabolite,contains.ARA,contains.DHA,dbl.cb.bonds,total.n.cb,tag_db,tagC) %>% 
  rename(metabolite_id = metabolite) %>%
  left_join(PC_loadings_616,join_by(metabolite_id))
head(TAG_fa)

##############################################
# 1. Load libraries
##############################################
library(tidyverse)
library(patchwork)

##############################################
# 2. Load your data
##############################################
# Replace with your file
# df <- read_csv("your_file.csv")

##############################################
# 3. RENAME COLUMNS (EDIT HERE)
##############################################
# Replace the right-hand side with your column names
df <- TAG_fa %>%
  rename(
    PC1_loading     = PC1,
    ARA_present     = contains.ARA,       # binary indicator (0/1 or TRUE/FALSE)
    DHA_present     = contains.DHA,       # binary indicator
    double_bonds    = tag_db,
    total_carbons   = tagC
  )
df = df %>% filter(!is.na(class))

##############################################
# 4. Create PC1 bins
##############################################
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

##############################################
# 5. Create category groups
##############################################

# Double bond groups
df <- df %>%
  mutate(
    db_group = case_when(
      double_bonds == 0 ~ "0",
      double_bonds == 1 ~ "1",
      double_bonds >= 2 & double_bonds <= 6 ~ "2–6",
      double_bonds > 6 ~ ">6"
    ),
    db_group = factor(db_group, levels = c("0", "1", "2–6", ">6"))
  )

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

##############################################
# 6. Compute proportions by bin and category
##############################################
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

##############################################
# 7. Combine all datasets into long format
##############################################
df_all <- bind_rows(
  df_ara %>% rename(category = category),
  df_dha %>% rename(category = category),
  df_db  %>% rename(category = db_group),
  df_carb %>% rename(category = C_group)
)

##############################################
# 8. Faceted bar plot (4 panels)
##############################################
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

#########################################sum##############################################
# 9. Heatmap
##############################################
p_heat <- ggplot(df_all, aes(x = PC1_bin, y = category, fill = proportion)) +
  geom_tile(color = "white") +
  facet_wrap(~feature, scales = "free_y") +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "PC1 loading bin", y = "Category", fill = "Proportion")

##############################################
# 10. Line plot option
##############################################
p_line <- ggplot(df_all, aes(x = PC1_bin, y = proportion,
                             color = category, group = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~feature, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "PC1 loading bin", y = "Proportion")

##############################################
# 11. Print all figures
##############################################
p_bar
p_heat
p_line
