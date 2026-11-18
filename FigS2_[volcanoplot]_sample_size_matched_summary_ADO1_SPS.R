# ----------------------------@
# Supplementary Figure S2
# Created: 2026-10-27
# Modified: 2026-11-17
# ----------------------------@

# 0. Setup ----
wd <- "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP"
setwd(wd)

library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
select = dplyr::select
describe = psych::describe

source(file.path(
  wd,
  "99_scripts",
  "[function]_define_fit_function_2025-01-09.R"
))

# 1. Results data frame ----
metabodata_info <- fread("data/metabodata_info_ados.txt")

load("results/sample_size_matched_summary.Rdata", verbose = TRUE)

traits <- c(
  "Volume.WM.adjAgeSexICV",
  "T1wSI.WM.adjAgeSex",
  "MTR.WM.adjAgeSex"
)

# 2. Add class information ----
volcano_main_all <- matched_summary %>%
  rename(metabolite_id = metabolite) %>%
  left_join(
    metabodata_info %>%
      mutate(metabolite_id = str_remove(metabolite_id, "_nmol/mL")) %>%
      select(metabolite_id, class),
    by = "metabolite_id"
  ) %>%
  mutate(
    trait = factor(trait, levels = traits),
    generation = "ado1_resample"
  )

# 3. Volcano plots ----
point.size <- 1
generation_plot <- "ado1_resample"

## 3A. adolescents (500 resamples) ----
load("results/signif_level_used_for_Fig1_apr2026.Rdata", verbose = TRUE)
print(signif.level)

p <- vector("list", length(traits))

for (ti in seq_along(traits)) {
  trait_plot <- traits[ti]
  
  p[[ti]] <- volcano_main_all %>%
    mutate(
      fdr.signif = if_else(
        mean_p < signif.level,
        "p.fdr<0.05",
        "p.fdr>=0.05"
      ),
      fdr.signif = factor(
        fdr.signif,
        levels = c("p.fdr<0.05", "p.fdr>=0.05")
      )
    ) %>%
    filter(trait == trait_plot, generation == generation_plot) %>%
    ggplot(aes(
      x = median_beta,
      y = -log10(mean_p),
      shape = fdr.signif,
      fill = fdr.signif,
      color = fdr.signif
    )) +
    geom_point(alpha = 0.5, size = point.size, stroke = 1) +
    scale_shape_manual(values = c(19, 21)) +
    scale_color_manual(values = c("#073B4C", "grey")) +
    scale_fill_manual(values = c("#073B4C", "grey")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(
        fill = scales::alpha("#f8edeb", 0.2),
        color = NA,
        linewidth = 0.5,
        linetype = "solid"
      ),
      panel.grid.major = element_line(
        linewidth = 0.5,
        linetype = "solid",
        colour = "grey95"
      ),
      panel.grid.minor = element_line(
        linewidth = 0.25,
        linetype = "solid",
        colour = "grey95"
      ),
      axis.title = element_text(size = 12),
      panel.border = element_rect(
        colour = "grey85",
        fill = NA,
        linewidth = 1
      )
    ) +
    geom_hline(
      yintercept = -log10(signif.level),
      linetype = 2,
      color = "darkred"
    ) +
    geom_vline(xintercept = 0, color = "darkgrey") +
    ggtitle(trait_plot)
}

p_resample_ados_adults <- wrap_plots(p, ncol = 3) +
  plot_layout(axes = "collect", guides = "collect") &
  coord_cartesian(ylim = c(0, 13), xlim = c(-0.3, 0.3)) &
  ggtitle(NULL) &
  xlab(NULL) &
  ylab(NULL) &
  theme(legend.position = "none")

ggsave(
  paste0(
    "results/[manuscript]_sample_size_matched_volcanoplots_ados_",
    Sys.Date(),
    ".png"
  ),
  plot = p_resample_ados_adults &
    theme(axis.text = element_text(size = 50)),
  width = 17.91 * 0.8,
  height = 13.26 * 0.8 * 0.5,
  units = "in",
  dpi = 300
)

## 3B. adults (same as those from the primary analyses) ----
load("results/separate_Fig1_volcanoplots.Rdata", verbose = TRUE)

p_resample_ados_adults <-
  volcano_WM.vol_adults +
  volcano_WM.T1wSI_adults +
  volcano_WM.MTR_adults +
  plot_layout(ncol = 3, axes = "collect", guides = "collect") &
  coord_cartesian(ylim = c(0, 13), xlim = c(-0.3, 0.3)) &
  ggtitle(NULL) &
  xlab(NULL) &
  ylab(NULL) &
  theme(legend.position = "none")

ggsave(
  paste0(
    "results/[manuscript]_sample_size_matched_volcanoplots_adults_",
    Sys.Date(),
    ".png"
  ),
  plot = p_resample_ados_adults &
    theme(axis.text = element_text(size = 50)),
  width = 17.91 * 0.8,
  height = 13.26 * 0.8 * 0.5,
  units = "in",
  dpi = 300
)

# 4. Pie plots ----
group <- "ados_resample"

for (trait_plot in traits) {
  plot.df <- volcano_main_all %>% filter(trait == trait_plot)
  
  ind_sig <- plot.df$mean_p < signif.level
  ind_pos <- plot.df$median_beta > 0
  
  nmetabo <- c(
    sum(ind_sig),
    sum(!ind_sig),
    sum(ind_sig & !ind_pos),
    sum(ind_sig & ind_pos)
  )
  
  mylabels <- c(
    "p<0.0189",
    "p>=0.0189",
    "beta<0; p<0.0189",
    "beta>0; p<0.0189"
  )
  
  df <- data.frame(
    group = mylabels,
    value = nmetabo
  )
  
  p_pi <- df %>%
    filter(group %in% c("p<0.0189", "p>=0.0189")) %>%
    ggplot(aes(x = "", y = value, fill = group)) +
    geom_bar(
      stat = "identity",
      width = 1,
      color = "white"
    ) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("#9babbe", "#dbd9da")) +
    theme_void() +
    theme(legend.position = "none")
  
  pie_fname <- paste0(
    "[manuscript]_Fig1_pie_volcano_",
    trait_plot,
    "_",
    group,
    "_",
    Sys.Date(),
    ".png"
  )
  
  ggsave(
    filename = pie_fname,
    plot = p_pi,
    width = 6,
    height = 5.5,
    units = "in",
    dpi = 300
  )
}
