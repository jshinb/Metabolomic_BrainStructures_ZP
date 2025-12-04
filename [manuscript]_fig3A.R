rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)
plot.df.volcano_main = fread('results/plot.df.volcano_main_2025-12-02.tsv')
metabodata_ado = fread('data/metabodata_ados.txt')
head(metabodata_ado)
names(metabodata_ado) = str_remove(names(metabodata_ado),"_nmol/mL")
length(intersect(names(metabodata_ado),unique(plot.df.volcano_main$metabolite_id)))

metabodata_ado = subset(metabodata_ado,select=c('uniqueID',unique(plot.df.volcano_main$metabolite_id)))

load('results/quad.lipids_NormWM_MTR_WMvol_assoc.Rdata')
siglipids = c(quad.lipids$quad1$metabolite_id,quad.lipids$quad2$metabolite_id,quad.lipids$quad4$metabolite_id)#616

cor.lipids = cor(metabodata_ado[,-1],use='p')
# load package
library(pheatmap)
pheatmap(cor.lipids)

load('results/adj.metabo.Rdata')
tail(names(adj.metabo))
cor.lipids_signif = cor(subset(adj.metabo,select=siglipids),use='p')
p_heatmap = pheatmap(cor.lipids_signif,fontsize = 0.1)
p_heatmap
ggsave('results/cor_heatmap_616_metabo_adjAgeSex_ados_2025-12-02.png',
       plot = p_heatmap,
       # width=11,height=8.5,units='in',
       width=22.5,height=13.35,units='cm',
       dpi=300)

# histogram ----
diag(cor.lipids_signif) = NA
r.metabo.adj = data.frame(r=c(cor.lipids_signif))
r.metabo.adj = na.omit(r.metabo.adj)

p_histogram = r.metabo.adj %>% ggplot(aes(x=r)) + 
  geom_histogram(binwidth = 0.1,color="white")

## Calculate quantiles
quantiles <- quantile(iris$Sepal.Length, probs = c(0.25, 0.5, 0.75))

## Create a data frame to help with labeling the quantile lines
quantile_df <- data.frame(
  quantile = quantiles,
  label = c("25th percentile", "50th percentile", "75th percentile")
)


## Estimate density
dens <- density(r.metabo.adj$r)
dd <- with(dens, data.frame(x, y))

## Calculate quantiles with bounding values
quantiles2 <- quantile(r.metabo.adj$r, probs = c(0, 0.25, 0.5, 0.75, 1))
#         0%        25%        50%        75%       100% 
# -0.4131992  0.2678381  0.4612028  0.6467400  0.9988038 

## Assign each observations to a quantile range
dd$quantile_range <- with(dd, cut(x, breaks = quantiles2, 
                                  labels = c("0-25%", "25-50%", "50-75%", "75-100%"),
                                  include.lowest = TRUE))
table(dd$quantile_range)

## Plotting
p_hist = ggplot() +
  geom_histogram(data = r.metabo.adj, aes(x =r),#aes(x =r, y = ..density..), 
                 binwidth = 0.1, fill = "grey", alpha = 0.25, color = "black") +
  # geom_line(data = dd, aes(x = x, y = y), size = 1.5) +
  # geom_ribbon(data = dd, aes(x = x, ymax = y, ymin = 0, fill = quantile_range), alpha = 0.5) +
  # scale_fill_brewer(palette = "Dark2", name = "Data Quantile Range") +
  
  # labs(title = "Distribution of pairwise Pearson correlation coefficients for age- and sex-adjusted metabolite levels",
       # x = "r", y = "Probability Density") +
  theme_minimal() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) + xlab('Correlation coefficient') + ylab('Frequency') + 
  geom_vline(xintercept =  quantiles2,color="darkred", linetype=2)
p_hist
ggsave('results/cor_histogram_616_metabo_adjAgeSex_ados_2025-12-02.png',
       plot = p_hist,
       width=6.5*1.0,height=3.5*1.0,units='in',
       dpi=300)
quantiles2
