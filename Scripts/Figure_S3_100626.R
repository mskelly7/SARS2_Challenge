# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# Supplemental Figure S3 – microbiome diversity & bacterial burden
# Last updated: 10/06/2026

script_summary <- "This script generates Figure S3 summarizing longitudinal microbiome alpha diversity and bacterial burden in nasal and 
throat samples. Nasal 16S rRNA gene sequencing and throat metagenomic sequencing data are imported as phyloseq objects and used to calculate 
Shannon diversity and observed richness (number of ASVs or taxa) across longitudinal timepoints, stratified by infection status. Bacterial 
burden is summarized using log10-transformed qPCR values across timepoints for nasal and throat samples."
cat(script_summary)

remove(list=ls())
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(123)

project_dir <- "________________________"
data_dir <- file.path(project_dir,"Data")
fig_dir <- file.path(project_dir,"Figures")
suppfig1_dir <- file.path(fig_dir,"FigureS3")

dir.create(suppfig1_dir, recursive=TRUE, showWarnings=FALSE)

all_timepoint <- c("M1","1AM","1PM","2AM","2PM","3AM","3PM","5AM","7AM","14AM")
timepoint_labels <- c(
  "M1"="-1am",
  "1AM"="1am",
  "1PM"="1pm",
  "2AM"="2am",
  "2PM"="2pm",
  "3AM"="3am",
  "3PM"="3pm",
  "5AM"="5am",
  "7AM"="7am",
  "14AM"="14am"
)
infection_colors <- c("Uninfected"="#388fcc","Infected"="#f68c34")
dataset_labels <- c("MTS"="Nasal","Throat"="Throat")

load_phy_mts <- function() {
  readRDS(file.path(data_dir,"phy.mts.16s.pruned.rds"))
}

load_phy_throat <- function() {
  phy.temp <- readRDS(file.path(data_dir,"phy.throat.ngs.pruned.rds"))
  otu_mat <- t(as(otu_table(phy.temp), "matrix"))
  colnames(otu_mat) <- paste0("Taxon", colnames(otu_mat))
  tax_mat <- as(tax_table(phy.temp), "matrix")
  rownames(tax_mat) <- paste0("Taxon", rownames(tax_mat))
  phyloseq(
    otu_table(otu_mat, taxa_are_rows=FALSE),
    sample_data(data.frame(sample_data(phy.temp))),
    tax_table(tax_mat)
  )
}

build_richness_df <- function(phy, dataset_label) {
  tab_richness <- estimate_richness(phy)[, c("Observed","Shannon")]
  tab_meta <- as.data.frame(sample_data(phy)[, c("studyID","infected","timepoint")])
  tab_meta$infection_status <- ifelse(tab_meta$infected == "N", "Uninfected", "Infected")
  tab_meta <- tab_meta[, c("studyID","infection_status","timepoint")]
  tab_richness <- merge(tab_meta, tab_richness, by=0) %>% column_to_rownames(var="Row.names")
  colnames(tab_richness) <- c("studyID","infection_status","timepoint","Richness","Shannon")
  tab_richness$infection_status <- factor(tab_richness$infection_status, levels=c("Uninfected","Infected"))
  tab_richness$timepoint <- factor(tab_richness$timepoint, levels=all_timepoint)
  tab_richness$dataset <- dataset_label
  tab_richness
}

build_qpcr_df <- function(phy, dataset_label) {
  df_metadata <- data.frame(sample_data(phy))
  df <- df_metadata[, c("sampleID","studyID","sampletype","timepoint","infected","qpcr_log10")]
  df$infection_status <- ifelse(df$infected == "N", "Uninfected", "Infected")
  df$infection_status <- factor(df$infection_status, levels=c("Uninfected","Infected"))
  df$timepoint <- factor(df$timepoint, levels=all_timepoint)
  df <- df[!is.na(df$qpcr_log10), ]
  df$sampletype <- dataset_label
  df
}

base_boxplot_theme <- theme_bw() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=11),
    strip.text=element_text(size=11, face="bold"),
    strip.background=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=10),
    panel.grid.major.x=element_blank(),
    panel.grid.minor=element_blank()
  )

phy_mts <- load_phy_mts()
phy_throat <- load_phy_throat()

mts_tab_richness <- build_richness_df(phy_mts, "MTS")
throat_tab_richness <- build_richness_df(phy_throat, "Throat")
all_tab_richness <- rbind(mts_tab_richness, throat_tab_richness)
all_tab_richness$dataset <- factor(all_tab_richness$dataset, levels=c("MTS","Throat"))

df_shannon <- all_tab_richness[, c("studyID","infection_status","timepoint","Shannon","dataset")]
df_shannon_long <- pivot_longer(df_shannon, cols="Shannon")

p_shannon <- df_shannon_long %>%
  ggplot() +
  geom_boxplot(aes(x=timepoint, y=value, color=infection_status)) +
  scale_color_manual(values=infection_colors) +
  scale_x_discrete(labels=timepoint_labels) +
  labs(x=NULL, y="Shannon diversity") +
  facet_wrap(~dataset, scales="free_y", labeller=labeller(dataset=dataset_labels)) +
  scale_y_continuous(limits=c(0,6)) +
  base_boxplot_theme

df_richness <- all_tab_richness[, c("studyID","infection_status","timepoint","Richness","dataset")]
df_richness_long <- pivot_longer(df_richness, cols="Richness")

richness_limits <- data.frame(
  dataset=factor(
    c("MTS","MTS","Throat","Throat"),
    levels=c("MTS","Throat")
  ),
  timepoint=factor(
    c("M1","M1","M1","M1"),
    levels=all_timepoint
  ),
  infection_status=factor(
    c("Uninfected","Uninfected","Uninfected","Uninfected"),
    levels=c("Uninfected","Infected")
  ),
  value=c(0,500,0,10000),
  name="Richness"
)

p_richness <- df_richness_long %>%
  ggplot() +
  geom_blank(
    data=richness_limits,
    aes(x=timepoint,y=value)
  ) +
  geom_boxplot(aes(x=timepoint,y=value,color=infection_status)) +
  scale_color_manual(values=infection_colors) +
  scale_x_discrete(labels=timepoint_labels) +
  labs(x=NULL,y="Observed richness") +
  facet_wrap(
    ~dataset,
    scales="free_y",
    labeller=labeller(dataset=dataset_labels)
  ) +
  scale_y_continuous(
    expand=expansion(mult=c(0.02,0.05))
  ) +
  base_boxplot_theme

df_mts_qpcr <- build_qpcr_df(phy_mts, "MTS")
df_throat_qpcr <- build_qpcr_df(phy_throat, "Throat")
df_ml <- rbind(df_mts_qpcr, df_throat_qpcr)
df_ml$sampletype <- factor(df_ml$sampletype, levels=c("MTS","Throat"))

p_qpcr <- df_ml %>%
  ggplot() +
  geom_boxplot(aes(x=timepoint, y=qpcr_log10, color=infection_status)) +
  scale_color_manual(values=infection_colors) +
  scale_x_discrete(labels=timepoint_labels) +
  labs(x=NULL, y=expression(Log[10]~qPCR)) +
  facet_wrap(~sampletype, scales="free_y", labeller=labeller(sampletype=dataset_labels)) +
  scale_y_continuous(limits=c(0,6)) +
  base_boxplot_theme

ggsave(file.path(suppfig1_dir,"FigureS3a.png"), p_shannon, width=7, height=4, dpi=600, bg="white")
ggsave(file.path(suppfig1_dir,"FigureS3b.png"), p_richness, width=7, height=4, dpi=600, bg="white")
ggsave(file.path(suppfig1_dir,"FigureS3c.png"), p_qpcr, width=7, height=4, dpi=600, bg="white")

figure_s3 <- cowplot::plot_grid(
  p_shannon,
  p_richness,
  p_qpcr,
  NULL,
  ncol=2,
  labels=c("A","B","C",""),
  label_size=12,
  align="hv"
)

ggsave(
  filename=file.path(fig_dir,"FigureS3.png"),
  plot=figure_s3,
  width=15,
  height=8,
  dpi=600,
  bg="white"
)

density_by_site <- df_ml %>%
  dplyr::group_by(sampletype) %>%
  dplyr::summarise(
    n=dplyr::n(),
    mean_log10_qpcr=mean(qpcr_log10,na.rm=TRUE),
    sd_log10_qpcr=sd(qpcr_log10,na.rm=TRUE),
    median_log10_qpcr=median(qpcr_log10,na.rm=TRUE),
    iqr_log10_qpcr=IQR(qpcr_log10,na.rm=TRUE),
    .groups="drop"
  )
cat("\n=== Mean bacterial density by site ===\n")
print(density_by_site,n=Inf)

density_fold_change <- density_by_site %>%
  dplyr::select(sampletype,mean_log10_qpcr) %>%
  tidyr::pivot_wider(
    names_from=sampletype,
    values_from=mean_log10_qpcr
  ) %>%
  dplyr::mutate(
    log10_difference_Throat_minus_MTS=Throat - MTS,
    fold_difference_Throat_vs_MTS=10^log10_difference_Throat_minus_MTS
  )
cat("\n=== Fold difference in mean bacterial density: Throat vs Nasal ===\n")
print(density_fold_change,n=Inf)

# Across all samples, mean bacterial density was approximately 51.4-fold higher in throat compared with nasal samples (Δlog10 density = 1.71).