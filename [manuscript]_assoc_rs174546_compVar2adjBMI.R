#-----------------------------------------------------------------@
# Created: 2026-10-23
# Test genetic association between FADS1 locus vs. compVar2adjBMI (PCA of WM.PC1 and Met.PC3)
#-----------------------------------------------------------------@

library(data.table)
library(readr)
library(dplyr)
library(ggplot2)

wd='~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP'
setwd(wd)

d_compVar2adjBMI = fread("gwas_results/regenie_results_20261023/compVar2adjBMI/gwas_compVar2adjBMI_ADO1_chr11_compVar3.adjBMI.rnt.regenie.gz")
head(d_compVar2adjBMI)
d_compVar2adjBMI %>% arrange(-LOG10P)
d_compVar2adjBMI %>% filter(ID=='rs174546') %>% mutate(pval = 10^(-LOG10P))
d_compVar2adjBMI %>% 
  ggplot()
topSNP = 'rs174546'

## 
LDproxies = fread("gwas_results/regenie_results_20261023/compVar2adjBMI/proxy_rs174546_2026-10-23.txt")
LDproxies  %>% head

rs174546_LDproxies = LDproxies%>% filter(R2>0.6) %>% pull(RS_Number) #61
d_compVar2adjBMI_rs174546_LDproxies = 
  d_compVar2adjBMI %>% filter(ID %in% rs174546_LDproxies) #52

d_compVar2adjBMI_rs174546_LDproxies %>%
  ggplot(aes(x=BETA, y=LOG10P)) + 
  geom_point()

## merge with LD
d_compVar2adjBMI %>% 
  right_join(LDproxies %>% rename(ID = RS_Number)) %>% 
  select(ID, BETA, LOG10P, R2) %>%
  mutate(
    R2_group = cut(
      R2,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
    )
  ) %>%
  ggplot(aes(x = BETA, y = LOG10P, color = R2_group)) + 
  geom_point(size = 1, alpha=0.5) +
  scale_color_manual(
    values = c(
      "0–0.2"   = "royalblue",
      "0.2–0.4" = "cyan2",
      "0.4–0.6" = "green3",
      "0.6–0.8" = "orange",
      "0.8–1.0" = "red"
    ),
    name = expression(r^2)
  ) +
  geom_point(
    data = ~ subset(.x, ID == topSNP),
    color = "darkred",
    size = 2,
    shape = 21
  ) +
  ggrepel::geom_text_repel(
    data = ~ subset(.x, ID == topSNP),
    aes(label = ID),
    color = "black",
    nudge_y = 0.15,
    box.padding = 0.5,
    point.padding = 0.5
  )
