# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# Figure S8 - nasal multi-omic analyses by viral load
# Figure S9 - throat multi-omic analyses by viral load
# Last updated: 07/11/2026

script_summary <- "This script generates Figures S8 and S9 summarizing longitudinal multi-omic signatures associated with viral load
among infected participants using previously fit DIABLO models and precomputed immune-module scores. Final nasal and throat viral-load
DIABLO models at the selected design strength (0.25) are imported across longitudinal timepoints, and component 1 microbiome and
immune-module loadings are extracted. Loadings are converted to signed values such that features associated with high viral load
are displayed as positive and features associated with low viral load are displayed as negative. Figure S8 presents nasal analyses,
including model performance, signed component 1 feature loadings, and microbe-host correlation networks. Figure S9 presents throat
model performance and signed component 1 feature loadings. The plotting architecture matches the finalized Figure 5/S7 script."
cat(script_summary)

set.seed(1234)
remove(list=ls())
version

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(mixOmics)
library(cowplot)
library(ggplotify)
library(grid)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)

project_dir <- "________________________"
data_dir <- file.path(project_dir,"Data")
diablo_results_dir <- file.path(project_dir,"DIABLO","results_chaussabel")

design_strength <- 0.25
design_label <- paste0("design_",gsub("\\.","p",as.character(design_strength)))

fig_dir <- file.path(project_dir,"Figures")
figs8_dir <- file.path(fig_dir,"FigureS8")
figs8_design_dir <- file.path(figs8_dir,design_label)
figs9_dir <- file.path(fig_dir,"FigureS9")
figs9_design_dir <- file.path(figs9_dir,design_label)

dir.create(figs8_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figs8_design_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figs9_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figs9_design_dir,recursive=TRUE,showWarnings=FALSE)

diablo_vl_dir <- file.path(diablo_results_dir,"viral_load")
diablo_vl_mts_dir <- file.path(diablo_vl_dir,"mts")
diablo_vl_throat_dir <- file.path(diablo_vl_dir,"throat")

gene_sets <- readRDS(file.path(data_dir,"chaussabel_modules_17.rds"))
module_order <- names(gene_sets)

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each=2),c("AM","PM")))[-c(9,11,13)]

loading_col_fun <- colorRamp2(c(-1,0,1),c("#388fcc","white","#F9A65A"))
ann_colors <- list(view=c(`Species`="#C7B9E2",`Immune modules`="#BFE3C0"))

perf_time_labels <- setNames(all_timepoint,all_timepoint)
perf_time_labels["M1"] <- "-1am"
perf_time_labels <- setNames(tolower(perf_time_labels),names(perf_time_labels))
perf_time_labels["M1"] <- "-1am"

perf_colors <- list(
  view=c(
    "Microbiome"="#7E6BBE",
    "Immune modules"="#5FAE6E"
  )
)

site_map <- list(
  nasal=list(
    site_label="MTS",
    base_dir=diablo_vl_mts_dir,
    prefix="MTS",
    ps_path=file.path(data_dir,"phy.mts.16s.pruned.species.rds"),
    ssgsea_path=file.path(data_dir,"ssgsea.rds")
  ),
  throat=list(
    site_label="Throat",
    base_dir=diablo_vl_throat_dir,
    prefix="Throat",
    ps_path=file.path(data_dir,"phy.throat.ngs.pruned.species.rds"),
    ssgsea_path=file.path(data_dir,"ssgsea.rds")
  )
)

# Define helper functions

normalize_strain_tokens <- function(x){
  x <- as.character(x)
  x <- gsub("sp\\.\\s*","sp. ",x,perl=TRUE)
  x <- gsub("sp\\.\\s*HMT","sp. HMT",x,perl=TRUE)
  x <- gsub("\\bKHUD\\s+","KHUD",x,perl=TRUE)
  x <- gsub("HMT-\\s*([0-9]+)","HMT\\1",x,perl=TRUE)
  x <- gsub("^([A-Z][a-zA-Z-]+)\\s+\\1;*s?$","\\1 sp.",x,perl=TRUE)
  x <- gsub("^([A-Z][a-zA-Z-]+);*$","\\1 sp.",x,perl=TRUE)
  x
}

escape_plotmath_text <- function(x){
  gsub("'","\\\\'",x,fixed=FALSE)
}

abbrev_genus <- function(genus){
  if(is.na(genus) || !nzchar(genus)) return(genus)
  paste0(substr(genus,1,1),".")
}

species_label_plotmath <- function(label,abbreviate_genus=TRUE){
  label <- normalize_strain_tokens(label)
  parts <- strsplit(trimws(label),"\\s+")[[1]]
  if(length(parts)==0) return("''")
  genus <- parts[1]
  rest <- parts[-1]
  always_abbrev <- c("Fusobacterium")
  if(length(rest)==0 || grepl("^;+?$",paste(rest,collapse=""))){
    if(tolower(genus) %in% tolower(always_abbrev)) genus <- abbrev_genus(genus)
    return(paste0("italic('",escape_plotmath_text(genus),"')~'sp.'"))
  }
  if(rest[1]=="sp."){
    if(tolower(genus) %in% tolower(always_abbrev)) genus <- abbrev_genus(genus)
    tail_txt <- paste(rest[-1],collapse=" ")
    tail_txt <- escape_plotmath_text(tail_txt)
    full_label <- paste(genus,"sp.",tail_txt)
    if(full_label=="Anaerococcus sp. HMT290") tail_txt <- paste0(tail_txt,"  ")
    if(nzchar(tail_txt)){
      return(paste0("italic('",escape_plotmath_text(genus),"')~'sp.'~'",tail_txt,"'"))
    }else{
      return(paste0("italic('",escape_plotmath_text(genus),"')~'sp.'"))
    }
  }
  if(tolower(genus) %in% tolower(always_abbrev)){
    genus <- abbrev_genus(genus)
  }else if(abbreviate_genus){
    genus <- abbrev_genus(genus)
  }
  species_epithet <- rest[1]
  tail_txt <- paste(rest[-1],collapse=" ")
  if(grepl("^[a-z]+$",species_epithet)){
    if(nzchar(tail_txt)){
      return(paste0("italic('",escape_plotmath_text(paste(genus,species_epithet)),"')~'",escape_plotmath_text(tail_txt),"'"))
    }else{
      return(paste0("italic('",escape_plotmath_text(paste(genus,species_epithet)),"')"))
    }
  }
  paste0("italic('",escape_plotmath_text(genus),"')~'",escape_plotmath_text(paste(rest,collapse=" ")),"'")
}

rnaseq_label_plotmath <- function(label){
  lbl <- as.character(label)
  if(grepl("^M[0-9]+\\.[0-9]+_",lbl)){
    lbl <- clean_module_label(lbl)
  }else{
    lbl <- normalize_strain_tokens(lbl)
  }
  paste0("'",escape_plotmath_text(lbl),"'")
}

feature_label_plotmath <- function(name,block,abbreviate_genus=TRUE){
  if(block=="Microbiome"){
    species_label_plotmath(name,abbreviate_genus=abbreviate_genus)
  }else{
    rnaseq_label_plotmath(name)
  }
}

display_time_label <- function(x){
  ifelse(x=="M1","-1am",tolower(x))
}

read_diablo_model <- function(base_dir,prefix,time,design_lab=NULL){
  if(is.null(design_lab)) design_lab <- design_label
  candidates <- c(
    file.path(base_dir,paste0(prefix,"_",time,"_",design_lab,"_final_model.rds")),
    file.path(base_dir,paste0(prefix,"_",time,"_final_model.rds"))
  )
  fp <- candidates[file.exists(candidates)][1]
  if(is.na(fp)) stop(paste0("DIABLO model file not found. Tried: ",paste(candidates,collapse="; ")))
  readRDS(fp)
}

get_comp1_loadings_tables <- function(model){
  p1 <- plotLoadings(model,comp=1,contrib="max",method="median",block="microbiome",plot=FALSE)
  df_microbiome <- p1$microbiome
  df_microbiome$feature <- rownames(df_microbiome)
  df_microbiome$view <- "microbiome"
  df_microbiome <- df_microbiome[,c("feature","importance","view","GroupContrib")]
  df_microbiome <- df_microbiome[
    !grepl("^Bacteria([; ]|$)",df_microbiome$feature) &
      !grepl("unclassified|uncultured",df_microbiome$feature,ignore.case=TRUE),,
    drop=FALSE
  ]
  colnames(df_microbiome) <- c("feature","loading","view","vl_cat")
  p2 <- plotLoadings(model,comp=1,contrib="max",method="median",block="rnaseq",plot=FALSE)
  df_rnaseq <- p2$rnaseq
  df_rnaseq$feature <- rownames(df_rnaseq)
  df_rnaseq$view <- "rnaseq"
  df_rnaseq <- df_rnaseq[,c("feature","importance","view","GroupContrib")]
  colnames(df_rnaseq) <- c("feature","loading","view","vl_cat")
  list(microbiome=df_microbiome,rnaseq=df_rnaseq)
}

build_loading_long_df <- function(base_dir,prefix,times=all_timepoint){
  list_loading <- vector("list",length(times))
  names(list_loading) <- times
  for(time in times){
    final_model <- read_diablo_model(base_dir,prefix,time)
    loading_tables <- get_comp1_loadings_tables(final_model)
    df_all <- rbind(loading_tables$microbiome,loading_tables$rnaseq)
    rownames(df_all) <- NULL
    df_all$loading <- ifelse(df_all$vl_cat=="low",-abs(df_all$loading),abs(df_all$loading))
    list_loading[[time]] <- df_all
  }
  bind_rows(list_loading,.id="Time")
}

species_labels_consecutive_genus <- function(labels){
  labels <- as.character(labels)
  out <- character(length(labels))
  prev_genus <- NA_character_
  for(i in seq_along(labels)){
    x <- normalize_strain_tokens(labels[i])
    parts <- strsplit(trimws(x),"\\s+")[[1]]
    genus <- if(length(parts)>=1) parts[1] else NA_character_
    use_abbrev <- !is.na(genus) && !is.na(prev_genus) && genus==prev_genus
    out[i] <- species_label_plotmath(x,abbreviate_genus=use_abbrev)
    prev_genus <- genus
  }
  out
}

get_heatmap_plain_labels <- function(loadings_long_df){
  df_annot_microbiome <- loadings_long_df %>%
    filter(view=="microbiome") %>%
    dplyr::select(feature,view) %>%
    distinct() %>%
    arrange(feature) %>%
    column_to_rownames("feature")
  df_annot_rnaseq <- loadings_long_df %>%
    filter(view=="rnaseq") %>%
    dplyr::select(feature,view) %>%
    distinct() %>%
    separate(feature,into=c("key","value"),sep="_",extra="merge",remove=FALSE) %>%
    arrange(value,key) %>%
    dplyr::select(-key,-value) %>%
    column_to_rownames("feature")
  annotation_row <- rbind(df_annot_microbiome,df_annot_rnaseq)
  annotation_row$view <- factor(
    annotation_row$view,
    levels=c("microbiome","rnaseq"),
    labels=c("Species","Immune modules")
  )
  loading_matrix <- loadings_long_df %>%
    dplyr::select(feature,Time,loading) %>%
    pivot_wider(names_from=Time,values_from=loading) %>%
    column_to_rownames("feature")
  loading_matrix <- loading_matrix[rownames(annotation_row),,drop=FALSE]
  raw_lab <- rownames(loading_matrix)
  lab_row <- ifelse(
    annotation_row$view=="Species",
    normalize_strain_tokens(gsub("_"," ",raw_lab)),
    clean_module_label(raw_lab)
  )
  lab_row <- gsub("^Corynebacterium\\s+Corynebacterium;?s$","Corynebacterium sp.",lab_row)
  lab_row <- gsub("^Streptococcus\\s+Streptococcus;?s$","Streptococcus sp.",lab_row)
  normalize_strain_tokens(lab_row)
}

get_accuracy_df <- function(site_label,base_dir,prefix){
  perf_results <- data.frame(Time=character(),Pred_dist=character(),BER=numeric(),stringsAsFactors=FALSE)
  for(time in all_timepoint){
    candidates <- c(
      file.path(base_dir,paste0(prefix,"_",time,"_",design_label,"_perf.rds")),
      file.path(base_dir,paste0(prefix,"_",time,"_perf.rds"))
    )
    path <- candidates[file.exists(candidates)][1]
    if(is.na(path)) stop(paste0("Performance file not found. Tried: ",paste(candidates,collapse="; ")))
    perf_diablo <- readRDS(path)
    perf_results <- rbind(
      perf_results,
      data.frame(Time=time,Pred_dist="Max",BER=perf_diablo$WeightedVote.error.rate$max.dist["Overall.BER","comp1"]),
      data.frame(Time=time,Pred_dist="Centroids",BER=perf_diablo$WeightedVote.error.rate$centroids.dist["Overall.BER","comp1"]),
      data.frame(Time=time,Pred_dist="Mahalanobis",BER=perf_diablo$WeightedVote.error.rate$mahalanobis.dist["Overall.BER","comp1"])
    )
  }
  perf_results$Time <- factor(perf_results$Time,levels=all_timepoint)
  perf_results$Pred_dist <- factor(perf_results$Pred_dist,levels=c("Mahalanobis","Centroids","Max"))
  df_maha <- perf_results %>% filter(Pred_dist=="Mahalanobis")
  df_maha$Accuracy <- 1 - df_maha$BER
  df_maha$Site <- site_label
  df_maha
}

get_auc_df <- function(site_label,base_dir,prefix){
  auc_results <- data.frame(Time=numeric(),Time_label=character(),View=character(),AUC=numeric(),stringsAsFactors=FALSE)
  i <- 0
  for(time in all_timepoint){
    model <- read_diablo_model(base_dir,prefix,time)
    i <- i + 1
    auroc_res <- auroc(model,roc.comp=1)
    auc_results <- rbind(
      auc_results,
      data.frame(Time=i,Time_label=time,View="Microbiome",AUC=auroc_res$microbiome[["comp1"]][,"AUC"]),
      data.frame(Time=i,Time_label=time,View="RNA-Seq",AUC=auroc_res$rnaseq[["comp1"]][,"AUC"])
    )
  }
  auc_results$Time_label <- factor(auc_results$Time_label,levels=all_timepoint)
  auc_results$View <- factor(
    auc_results$View,
    levels=c("Microbiome","RNA-Seq"),
    labels=c("Microbiome","Immune modules")
  )
  auc_results$Site <- site_label
  auc_results
}

make_perf_theme <- function(
    x_text_size=9.0,
    y_text_size=9.0,
    y_title_size=11,
    legend_text_size=10,
    plot_margin=margin(8,10,8,10)
){
  theme_bw() +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      axis.text.x=element_text(angle=45,hjust=1,size=x_text_size,color="black"),
      axis.text.y=element_text(size=y_text_size,color="black"),
      axis.title.y=element_text(size=y_title_size,face="bold",margin=margin(r=3)),
      axis.title.x=element_blank(),
      legend.title=element_blank(),
      legend.text=element_text(size=legend_text_size),
      legend.position="bottom",
      legend.key.width=unit(0.7,"cm"),
      legend.key.height=unit(0.25,"cm"),
      plot.margin=plot_margin
    )
}

make_auc_plot <- function(auc_df,show_legend=TRUE,perf_theme=make_perf_theme()){
  ggplot(auc_df,aes(x=Time_label,y=AUC,group=View,color=View)) +
    geom_line(linewidth=0.65) +
    geom_point(size=1.8) +
    scale_color_manual(values=perf_colors$view) +
    scale_x_discrete(labels=perf_time_labels) +
    coord_cartesian(ylim=c(0.5,1.0)) +
    labs(x=NULL,y="AUC") +
    perf_theme +
    theme(
      legend.position=if(show_legend) "top" else "none",
      legend.direction="horizontal",
      legend.box.margin=margin(0,0,-6,0),
      legend.spacing.x=unit(0.05,"cm"),
      legend.key.width=unit(0.45,"cm"),
      legend.key.height=unit(0.20,"cm")
    )
}

make_accuracy_plot <- function(acc_df,perf_theme=make_perf_theme()){
  ggplot(acc_df,aes(x=Time,y=Accuracy,group=1)) +
    geom_line(linewidth=0.65,color="#619CFF") +
    geom_point(size=2.0,color="#619CFF") +
    scale_x_discrete(labels=perf_time_labels) +
    coord_cartesian(ylim=c(0,1)) +
    labs(x=NULL,y="Accuracy") +
    perf_theme +
    theme(legend.position="none")
}

build_accuracy_panel <- function(
    auc_df,
    acc_df,
    show_auc_legend=TRUE,
    auc_theme=make_perf_theme(),
    acc_theme=make_perf_theme(),
    rel_heights=c(1,0.90)
){
  plot_grid(
    make_auc_plot(auc_df,show_legend=show_auc_legend,perf_theme=auc_theme),
    make_accuracy_plot(acc_df,perf_theme=acc_theme),
    ncol=1,
    align="v",
    rel_heights=rel_heights
  )
}

build_heatmap_plot <- function(loadings_long_df,
                               heatmap_timepoints=NULL,
                               row_label_width=NULL,
                               fontsize_row=9.0,
                               fontsize_col=9.0,
                               row_title_fontsize=11.0,
                               draw_padding=grid::unit(c(0,0,10,0),"pt")){
  if(!is.null(heatmap_timepoints)){
    loadings_long_df <- loadings_long_df %>%
      dplyr::filter(Time %in% heatmap_timepoints)
  }
  df_annot_microbiome <- loadings_long_df %>%
    filter(view=="microbiome") %>%
    dplyr::select(feature,view) %>%
    distinct() %>%
    arrange(feature) %>%
    column_to_rownames("feature")
  df_annot_rnaseq <- loadings_long_df %>%
    filter(view=="rnaseq") %>%
    dplyr::select(feature,view) %>%
    distinct() %>%
    mutate(feature=factor(feature,levels=module_order)) %>%
    arrange(feature) %>%
    column_to_rownames("feature")
  annotation_row <- rbind(df_annot_microbiome,df_annot_rnaseq)
  annotation_row$view <- factor(
    annotation_row$view,
    levels=c("microbiome","rnaseq"),
    labels=c("Species","Immune modules")
  )
  loading_matrix <- loadings_long_df %>%
    dplyr::select(feature,Time,loading) %>%
    pivot_wider(names_from=Time,values_from=loading) %>%
    column_to_rownames("feature")
  colnames(loading_matrix) <- display_time_label(colnames(loading_matrix))
  loading_matrix <- loading_matrix[rownames(annotation_row),,drop=FALSE]
  raw_lab <- rownames(loading_matrix)
  lab_row <- ifelse(
    annotation_row$view=="Species",
    normalize_strain_tokens(gsub("_"," ",raw_lab)),
    clean_module_label(raw_lab)
  )
  lab_row <- gsub("^Corynebacterium\\s+Corynebacterium;?s$","Corynebacterium sp.",lab_row)
  lab_row <- gsub("^Streptococcus\\s+Streptococcus;?s$","Streptococcus sp.",lab_row)
  lab_row <- normalize_strain_tokens(lab_row)
  labels_row_pm <- character(length(lab_row))
  species_idx <- annotation_row$view=="Species"
  labels_row_pm[species_idx] <- species_labels_consecutive_genus(lab_row[species_idx])
  labels_row_pm[!species_idx] <- vapply(
    lab_row[!species_idx],
    rnaseq_label_plotmath,
    character(1)
  )
  htmp <- ComplexHeatmap::pheatmap(
    as.matrix(loading_matrix),
    labels_row=parse(text=labels_row_pm),
    main="",
    fontsize_row=fontsize_row,
    fontsize_col=fontsize_col,
    angle_col="45",
    cluster_cols=FALSE,
    cluster_rows=FALSE,
    annotation_row=annotation_row,
    annotation_legend=FALSE,
    annotation_names_row=FALSE,
    show_colnames=TRUE,
    color=loading_col_fun,
    annotation_colors=ann_colors,
    row_split=annotation_row$view,
    row_names_max_width=row_label_width,
    heatmap_legend_param=list(title=NULL),
    scale="none",
    border_color="grey60",
    na_col="gray80",
    row_title_gp=grid::gpar(
      fontsize=row_title_fontsize,
      fontface="bold"
    )
  )
  grob <- grid::grid.grabExpr({
    grid::grid.newpage()
    ComplexHeatmap::draw(
      htmp,
      show_heatmap_legend=FALSE,
      show_annotation_legend=FALSE,
      padding=draw_padding
    )
  })
  list(
    plot=ggplotify::as.ggplot(grob),
    matrix=loading_matrix
  )
}

make_loading_legend_plot <- function(
    legend_fontsize=9.5,
    legend_label_fontsize=7.5,
    legend_width=0.22,
    legend_height=0.18,
    legend_x=0.50,
    bar_y=0.40,
    label_y=0.16,
    title_y=0.72
){
  grob <- grid::grid.grabExpr({
    grid::grid.newpage()
    bar_cols <- loading_col_fun(seq(-1,1,length.out=200))
    grid::grid.raster(
      matrix(bar_cols,nrow=1),
      x=grid::unit(legend_x,"npc"),
      y=grid::unit(bar_y,"npc"),
      width=grid::unit(legend_width,"npc"),
      height=grid::unit(legend_height,"npc"),
      interpolate=TRUE
    )
    tick_vals <- c(-1,-0.5,0,0.5,1)
    tick_x <- legend_x + legend_width * ((tick_vals + 1) / 2 - 0.5)
    grid::grid.text(
      label=c("-1.0","-0.5","0","0.5","1.0"),
      x=grid::unit(tick_x,"npc"),
      y=grid::unit(label_y,"npc"),
      gp=grid::gpar(fontsize=legend_label_fontsize)
    )
    grid::grid.text(
      "Low viral load",
      x=grid::unit(legend_x - legend_width / 2 - 0.015,"npc"),
      y=grid::unit(title_y,"npc"),
      just="right",
      gp=grid::gpar(fontsize=legend_fontsize,fontface="bold")
    )
    grid::grid.text(
      "High viral load",
      x=grid::unit(legend_x + legend_width / 2 + 0.015,"npc"),
      y=grid::unit(title_y,"npc"),
      just="left",
      gp=grid::gpar(fontsize=legend_fontsize,fontface="bold")
    )
  })
  ggplotify::as.ggplot(grob)
}

get_comp1_selected_features <- function(model,block1="microbiome",block2="rnaseq"){
  loading_tables <- get_comp1_loadings_tables(model)
  df1 <- loading_tables[[block1]]
  df2 <- loading_tables[[block2]]
  if(is.null(df1) || is.null(df2)) stop("Requested blocks not found in comp 1 loading tables.")
  if(nrow(df1)==0 || nrow(df2)==0) stop("No selected variables for one of the blocks.")
  list(block1=trimws(df1$feature),block2=trimws(df2$feature))
}

make_log10_species_matrix_for_timepoint <- function(timepoint,selected_species,ps_path,pseudocount=1e-6,infected_only=FALSE){
  selected_species <- trimws(selected_species)
  ps_obj <- readRDS(ps_path)
  meta_df <- as.data.frame(sample_data(ps_obj))
  if(!("timepoint" %in% colnames(meta_df))) stop("Column 'timepoint' not found in sample_data")
  if(infected_only && !("infected" %in% colnames(meta_df))) stop("Column 'infected' not found in sample_data")
  samples_tp <- rownames(meta_df[meta_df$timepoint==timepoint,,drop=FALSE])
  if(infected_only){
    infected_status <- tolower(trimws(as.character(meta_df$infected)))
    infected_samples <- rownames(meta_df[infected_status %in% c("y","yes","true","1","infected","positive","pos"),,drop=FALSE])
    samples_tp <- intersect(samples_tp,infected_samples)
  }
  if(length(samples_tp)==0) stop(paste0("No samples remain for timepoint ",timepoint," after filtering"))
  ps_tp <- prune_samples(samples_tp,ps_obj)
  ps_rel <- transform_sample_counts(ps_tp,function(x) x/sum(x))
  otu_df <- as.data.frame(otu_table(ps_rel))
  if(!taxa_are_rows(ps_rel)) otu_df <- t(otu_df)
  otu_mat <- as.matrix(otu_df)
  tax_df <- as.data.frame(tax_table(ps_rel))
  tax_df <- tax_df[rownames(otu_mat),,drop=FALSE]
  tax_df$feature <- trimws(as.character(tax_df$Species))
  tax_df$feature[is.na(tax_df$feature) | tax_df$feature==""] <- rownames(tax_df)[is.na(tax_df$feature) | tax_df$feature==""]
  feature_norm <- normalize_strain_tokens(tax_df$feature)
  selected_norm <- normalize_strain_tokens(selected_species)
  out <- vapply(seq_along(selected_species),function(ii){
    sp <- selected_species[ii]
    sp_norm <- selected_norm[ii]
    idx <- which(tax_df$feature==sp)
    if(length(idx)==0) idx <- which(feature_norm==sp_norm)
    if(length(idx)==0) idx <- which(rownames(tax_df)==sp)
    if(length(idx)==0) idx <- which(grepl(sp_norm,feature_norm,ignore.case=TRUE,fixed=TRUE))
    if(length(idx)==0) return(rep(NA_real_,ncol(otu_mat)))
    rel_abund <- colSums(otu_mat[idx,,drop=FALSE])
    log10(rel_abund + pseudocount)
  },numeric(ncol(otu_mat)))
  out <- as.data.frame(out)
  rownames(out) <- colnames(otu_mat)
  colnames(out) <- selected_species
  out
}

make_ssgsea_module_matrix_for_timepoint <- function(
    timepoint,
    selected_modules,
    site_label,
    ssgsea_path=file.path(data_dir,"ssgsea.rds")
){
  ssgsea <- readRDS(ssgsea_path)
  ssgsea_tp <- ssgsea %>%
    dplyr::filter(
      site == site_label,
      timepoint == !!timepoint
    )
  missing_modules <- setdiff(
    selected_modules,
    unique(ssgsea_tp$module)
  )
  if(length(missing_modules) > 0){
    stop(
      paste(
        "Selected modules not found:",
        paste(missing_modules,collapse=", ")
      )
    )
  }
  ssgsea_tp %>%
    dplyr::select(sample,module,ssgsea_scaled) %>%
    tidyr::pivot_wider(
      names_from=module,
      values_from=ssgsea_scaled
    ) %>%
    tibble::column_to_rownames("sample") %>%
    as.data.frame() %>%
    dplyr::select(all_of(selected_modules))
}

build_microbe_host_correlation_edge_df_from_comp1 <- function(
    model,
    timepoint,
    ps_path,
    site_label,
    ssgsea_path=file.path(data_dir,"ssgsea.rds"),
    top_n=25,
    abs_cutoff=NULL,
    cor_method="pearson",
    pseudocount=1e-6,
    block1="microbiome",
    block2="rnaseq",
    infected_only=TRUE
){
  selected <- get_comp1_selected_features(model,block1=block1,block2=block2)
  species_mat <- make_log10_species_matrix_for_timepoint(
    timepoint=timepoint,
    selected_species=selected$block1,
    ps_path=ps_path,
    pseudocount=pseudocount,
    infected_only=infected_only
  )
  module_mat <- make_ssgsea_module_matrix_for_timepoint(
    timepoint=timepoint,
    selected_modules=selected$block2,
    site_label=site_label,
    ssgsea_path=ssgsea_path
  )
  common_samples <- intersect(rownames(species_mat),rownames(module_mat))
  if(length(common_samples)<3) stop("Fewer than 3 overlapping samples were found")
  species_mat <- species_mat[common_samples,,drop=FALSE]
  module_mat <- module_mat[common_samples,,drop=FALSE]
  edge_list <- list()
  k <- 1
  for(sp in colnames(species_mat)){
    for(mod in colnames(module_mat)){
      x <- species_mat[[sp]]
      y <- module_mat[[mod]]
      ok <- is.finite(x) & is.finite(y)
      if(sum(ok)>=3){
        rho <- suppressWarnings(cor(x[ok],y[ok],method=cor_method,use="complete.obs"))
        if(is.finite(rho)){
          edge_list[[k]] <- data.frame(
            from=sp,
            to=mod,
            weight=rho,
            w_abs=abs(rho),
            from_block=block1,
            to_block=block2,
            stringsAsFactors=FALSE
          )
          k <- k + 1
        }
      }
    }
  }
  df <- bind_rows(edge_list)
  if(nrow(df)==0) return(df)
  if(!is.null(abs_cutoff)) df <- df[df$w_abs>=abs_cutoff,,drop=FALSE]
  df <- df[order(df$w_abs,decreasing=TRUE),,drop=FALSE]
  if(!is.null(top_n)) df <- df[seq_len(min(nrow(df),top_n)),,drop=FALSE]
  rownames(df) <- NULL
  df
}

plot_signed_network_panel <- function(edge_df,title=NULL,pos_col="#E54B4B",neg_col="#3B63FF",
                                      node_micro="#C7B9E2",node_rna="#BFE3C0",
                                      node_size=3.0,label_size=3.0,title_size=12,
                                      plot_y=0,plot_height=0.93,title_y=0.96){
  if(nrow(edge_df)==0) return(ggplot2::ggplot()+ggplot2::theme_void())
  edge_df$sign <- ifelse(edge_df$weight>=0,"Positive","Negative")
  edge_df$w_abs <- abs(edge_df$weight)
  edge_df$layout_w <- 1/edge_df$w_abs
  edge_df$layout_w[!is.finite(edge_df$layout_w)] <- 1e6
  node_df1 <- unique(edge_df[,c("from","from_block"),drop=FALSE])
  colnames(node_df1) <- c("name","block")
  node_df2 <- unique(edge_df[,c("to","to_block"),drop=FALSE])
  colnames(node_df2) <- c("name","block")
  node_df <- unique(rbind(node_df1,node_df2))
  node_df$block <- ifelse(node_df$block=="microbiome","Microbiome",ifelse(node_df$block=="rnaseq","RNA-seq",node_df$block))
  node_df$label_parse <- mapply(feature_label_plotmath,node_df$name,node_df$block,MoreArgs=list(abbreviate_genus=TRUE),USE.NAMES=FALSE)
  g <- igraph::graph_from_data_frame(edge_df,directed=FALSE,vertices=node_df)
  g_tbl <- as_tbl_graph(g)
  if(nrow(node_df)>=6 && nrow(edge_df)>=6){
    lay <- suppressWarnings(create_layout(g_tbl,layout="stress",weights=layout_w))
  }else{
    lay <- create_layout(g_tbl,layout="fr",weights=w_abs)
  }
  lay_df <- as.data.frame(lay)
  p <- ggraph(lay) +
    geom_edge_link(aes(edge_width=w_abs,edge_colour=sign),alpha=0.8,lineend="round",show.legend=FALSE) +
    geom_node_point(aes(fill=block),shape=21,size=node_size,stroke=0.7,color="black",show.legend=FALSE) +
    ggrepel::geom_text_repel(
      data=lay_df,
      aes(x=x,y=y,label=label_parse),
      parse=TRUE,
      size=label_size,
      box.padding=0.3,
      point.padding=0.2,
      force=1.5,
      force_pull=0.5,
      max.overlaps=Inf,
      max.iter=40000,
      max.time=10,
      segment.color="grey20",
      segment.size=0.18,
      min.segment.length=0.1,
      bg.color="white",
      bg.r=0.12,
      inherit.aes=FALSE,
      show.legend=FALSE
    ) +
    scale_edge_width(range=c(0.2,1.2),limits=c(0,1),guide="none") +
    scale_edge_colour_manual(values=c("Positive"=pos_col,"Negative"=neg_col),guide="none") +
    scale_fill_manual(values=c("Microbiome"=node_micro,"RNA-seq"=node_rna),guide="none") +
    scale_x_continuous(expand=expansion(mult=0.04)) +
    scale_y_continuous(expand=expansion(mult=0.04)) +
    coord_cartesian(clip="off") +
    theme_void() +
    theme(plot.margin=ggplot2::margin(0,0,0,0))
  if(is.null(title) || !nzchar(title)) return(p)
  cowplot::ggdraw() +
    cowplot::draw_plot(p,x=0,y=plot_y,width=1,height=plot_height) +
    cowplot::draw_label(title,x=0.5,y=title_y,hjust=0.5,vjust=1,fontface="bold",size=title_size)
}

build_diablo_signed_network_panel_plot <- function(base_dir,prefix,time,title=NULL,abs_cutoff=NULL,top_n=25,block1="microbiome",block2="rnaseq",
                                                   ps_path,site_label,ssgsea_path=file.path(data_dir,"ssgsea.rds"),
                                                   cor_method="pearson",pseudocount=1e-6,pos_col="#E54B4B",neg_col="#3B63FF",
                                                   node_size=3.0,label_size=3.0,title_size=12,
                                                   plot_y=0,plot_height=0.93,title_y=0.96){
  model <- read_diablo_model(base_dir,prefix,time)
  edge_df <- build_microbe_host_correlation_edge_df_from_comp1(
    model=model,
    timepoint=time,
    ps_path=ps_path,
    site_label=site_label,
    ssgsea_path=ssgsea_path,
    top_n=top_n,
    abs_cutoff=abs_cutoff,
    cor_method=cor_method,
    pseudocount=pseudocount,
    block1=block1,
    block2=block2
  )
  plot_signed_network_panel(
    edge_df=edge_df,
    title=title,
    pos_col=pos_col,
    neg_col=neg_col,
    node_size=node_size,
    label_size=label_size,
    title_size=title_size,
    plot_y=plot_y,
    plot_height=plot_height,
    title_y=title_y
  )
}

build_network_panel_vertical <- function(base_dir,prefix,times,titles,abs_cutoff=NULL,top_n=25,
                                         ps_path,site_label,ssgsea_path=file.path(data_dir,"ssgsea.rds"),cor_method="pearson",
                                         node_size=3.0,network_label_size=3.0,title_size=12,
                                         plot_y=0,plot_height=0.93,title_y=0.96,
                                         plot_margin=margin(0,0,0,0,"mm"),
                                         panel_labels=NULL,panel_label_size=10){
  plot_list <- lapply(times,function(tt){
    build_diablo_signed_network_panel_plot(
      base_dir=base_dir,
      prefix=prefix,
      time=tt,
      title=unname(titles[tt]),
      abs_cutoff=abs_cutoff,
      top_n=top_n,
      ps_path=ps_path,
      site_label=site_label,
      ssgsea_path=ssgsea_path,
      cor_method=cor_method,
      node_size=node_size,
      label_size=network_label_size,
      title_size=title_size,
      plot_y=plot_y,
      plot_height=plot_height,
      title_y=title_y
    ) +
      theme(plot.margin=plot_margin)
  })
  cowplot::plot_grid(
    plotlist=plot_list,
    ncol=1,
    align="v",
    axis="lr",
    rel_heights=rep(1,length(plot_list)),
    labels=panel_labels,
    label_size=panel_label_size,
    label_fontface="bold"
  )
}

clean_module_label <- function(x){
  gsub("_"," ",x)
}

# Generate input data

mts_loading <- build_loading_long_df(diablo_vl_mts_dir,"MTS")
throat_loading <- build_loading_long_df(diablo_vl_throat_dir,"Throat")

mts_loading_filt <- mts_loading
throat_loading_filt <- throat_loading

all_heatmap_labels <- c(
  get_heatmap_plain_labels(mts_loading_filt),
  get_heatmap_plain_labels(throat_loading_filt)
)

shared_row_label_width <- ComplexHeatmap::max_text_width(
  all_heatmap_labels,
  gp=grid::gpar(fontsize=9.5)
) + grid::unit(2,"mm")

mts_acc <- get_accuracy_df("Nasal",diablo_vl_mts_dir,"MTS")
throat_acc <- get_accuracy_df("Throat",diablo_vl_throat_dir,"Throat")
mts_auc <- get_auc_df("Nasal",diablo_vl_mts_dir,"MTS")
throat_auc <- get_auc_df("Throat",diablo_vl_throat_dir,"Throat")

# Construct model-performance panels

fig_s8a <- build_accuracy_panel(
  auc_df=mts_auc,
  acc_df=mts_acc,
  show_auc_legend=TRUE,
  auc_theme=make_perf_theme(
    x_text_size=8.0,
    y_text_size=8.0,
    y_title_size=9.5,
    legend_text_size=9,
    plot_margin=margin(6,15,4,8)
  ),
  acc_theme=make_perf_theme(
    x_text_size=8,
    y_text_size=8,
    y_title_size=9.5,
    plot_margin=margin(4,15,6,8)
  ),
  rel_heights=c(1,0.90)
)

fig_s9a <- cowplot::plot_grid(
  build_accuracy_panel(
    auc_df=throat_auc,
    acc_df=throat_acc,
    show_auc_legend=TRUE,
    auc_theme=make_perf_theme(
      x_text_size=8.0,
      y_text_size=8.0,
      y_title_size=9.5,
      legend_text_size=8.7,
      plot_margin=margin(6,15,4,8)
    ),
    acc_theme=make_perf_theme(
      x_text_size=8,
      y_text_size=8,
      y_title_size=9.5,
      plot_margin=margin(4,15,6,8)
    ),
    rel_heights=c(1,0.90)
  ),
  NULL,
  ncol=1,
  rel_heights=c(2,1)
)

# Construct heatmap panels

mts_heat <- build_heatmap_plot(
  mts_loading_filt,
  heatmap_timepoints=NULL,
  row_label_width=shared_row_label_width,
  fontsize_row=8.0,
  fontsize_col=8.0,
  row_title_fontsize=9.5,
  draw_padding=grid::unit(c(0,0,10,0),"pt")
)

fig_s8b_legend <- make_loading_legend_plot(
  legend_fontsize=8.8,
  legend_label_fontsize=7.3,
  legend_width=0.25,
  legend_height=0.35,
  legend_x=0.50,
  bar_y=0.55,
  label_y=0.13,
  title_y=0.55
)

fig_s8b <- cowplot::plot_grid(
  fig_s8b_legend,
  mts_heat$plot,
  ncol=1,
  rel_heights=c(0.05,1.00),
  align="v"
)

throat_heat <- build_heatmap_plot(
  throat_loading_filt,
  heatmap_timepoints=NULL,
  row_label_width=shared_row_label_width,
  fontsize_row=8.0,
  fontsize_col=8.0,
  row_title_fontsize=9.5,
  draw_padding=grid::unit(c(0,0,10,0),"pt")
)

fig_s9b_legend <- make_loading_legend_plot(
  legend_fontsize=9.0,
  legend_label_fontsize=7.5,
  legend_width=0.25,
  legend_height=0.35,
  legend_x=0.50,
  bar_y=0.55,
  label_y=0.13,
  title_y=0.55
)

fig_s9b <- cowplot::plot_grid(
  fig_s9b_legend,
  throat_heat$plot,
  ncol=1,
  rel_heights=c(0.06,1.00),
  align="v"
)

# Construct Figure S8 nasal network panels

network_timepoint_vl <- c("1AM","7AM")
network_titles_vl <- c(`1AM`="1am",`7AM`="7am")

fig_s8c <- build_network_panel_vertical(
  base_dir=diablo_vl_mts_dir,
  prefix="MTS",
  times=network_timepoint_vl,
  titles=network_titles_vl,
  abs_cutoff=NULL,
  top_n=25,
  ps_path=file.path(data_dir,"phy.mts.16s.pruned.species.rds"),
  site_label="MTS",
  ssgsea_path=file.path(data_dir,"ssgsea.rds"),
  cor_method="pearson",
  node_size=2.3,
  network_label_size=2.6,
  title_size=10.0,
  plot_height=0.92,
  title_y=0.98,
  plot_margin=margin(0,0,0,0,"mm"),
  panel_labels=c("C","D"),
  panel_label_size=11
)

# Assemble Figure S8

figure_s8 <- cowplot::plot_grid(
  fig_s8a,
  fig_s8b,
  fig_s8c,
  ncol=3,
  labels=c("A","B",""),
  label_size=11,
  label_fontface="bold",
  rel_widths=c(0.80,1.05,0.80),
  align="h"
)

ggsave(filename=file.path(fig_dir,"FigureS8.png"),plot=figure_s8,width=11.5,height=7.0,dpi=600,bg="white")
ggsave(filename=file.path(figs8_dir,"FigureS8a.png"),plot=fig_s8a,width=3.5,height=7.0,dpi=600,bg="white")
ggsave(filename=file.path(figs8_dir,"FigureS8b.png"),plot=fig_s8b,width=4.5,height=7.0,dpi=600,bg="white")
ggsave(filename=file.path(figs8_dir,"FigureS8c.png"),plot=fig_s8c,width=3.5,height=7.0,dpi=600,bg="white")

# Assemble Figure S9

figure_s9 <- cowplot::plot_grid(
  fig_s9a,
  fig_s9b,
  ncol=2,
  labels=c("A","B"),
  label_size=11,
  label_fontface="bold",
  rel_widths=c(0.75,1.02),
  align="h"
)

ggsave(
  filename=file.path(fig_dir,"FigureS9.png"),
  plot=figure_s9,
  width=8.0,
  height=7.0,
  dpi=600,
  bg="white"
)

ggsave(filename=file.path(figs9_dir,"FigureS9a.png"),plot=fig_s9a,width=3.5,height=7.0,dpi=600,bg="white")
ggsave(filename=file.path(figs9_dir,"FigureS9b.png"),plot=fig_s9b,width=5.5,height=7.0,dpi=600,bg="white")

# Save individual networks for all timepoints

for(site_name in names(site_map)){
  base_dir <- site_map[[site_name]]$base_dir
  prefix <- site_map[[site_name]]$prefix
  ps_path <- site_map[[site_name]]$ps_path
  site_label <- site_map[[site_name]]$site_label
  ssgsea_path <- site_map[[site_name]]$ssgsea_path
  for(tp in all_timepoint){
    tp_label <- display_time_label(tp)
    p <- build_diablo_signed_network_panel_plot(
      base_dir=base_dir,
      prefix=prefix,
      time=tp,
      title=tp_label,
      abs_cutoff=NULL,
      top_n=25,
      ps_path=ps_path,
      site_label=site_label,
      ssgsea_path=ssgsea_path,
      cor_method="pearson"
    )
    out_dir <- ifelse(site_name=="nasal",figs8_design_dir,figs9_design_dir)
    ggsave(
      filename=file.path(out_dir,paste0(prefix,"_",tp,"_network.png")),
      plot=p,
      width=4,
      height=3.5,
      dpi=600,
      bg="white"
    )
  }
}