# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# Figure 4 - nasal immune module expression by symptom burden
# Figure S6 - throat immune module expression by symptom burden
# Last updated: 10/07/2026

script_summary <- "This script generates Figures 4 and S6, which summarize baseline and longitudinal immune module expression patterns by symptom burden 
among infected participants in nasal and throat samples, respectively. Previously computed Chaussabel ssGSEA module scores generated from TMM-normalized 
logCPM RNA-seq data, scaled to a 0-1 range within each anatomical site, are imported and used for all downstream analyses. Analyses are restricted to infected 
participants with low or high symptom burden. Mean module expression profiles are computed separately for low- and high-symptom participants across longitudinal 
timepoints. Radar plots display baseline (day -1am) immune module expression differences between low- and high-symptom groups, while heatmaps display 
within-group percent changes in scaled module expression relative to baseline over time. Longitudinal trajectory plots are generated for selected immune 
modules to visualize temporal patterns in individual participants and group means. Baseline module-specific differences shown in radar plots are estimated
from the longitudinal linear mixed-effects model as the high-versus-low symptom burden contrast at day -1am, with false discovery rate correction
across modules within the baseline timepoint. Overall baseline immune profile differences between symptom burden groups are additionally assessed using 
Euclidean-distance PERMANOVA performed on scaled ssGSEA values. Longitudinal changes in module expression are evaluated using linear mixed-effects models of 
scaled ssGSEA values adjusted for age and sex, incorporating timepoint, symptom burden category, and subject-level random effects. Heatmap star annotations 
are based on FDR correction performed within each symptom-burden group and timepoint across immune modules. The script also performs timepoint-specific 
comparisons between symptom burden groups and interaction analyses to identify differential longitudinal immune trajectories between low- and 
high-symptom participants."
cat(script_summary)

remove(list=ls())
version
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(fmsb)
library(cowplot)
library(ggplotify)
library(grid)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(vegan)

set.seed(123)

project_dir <- "________________________"
data_dir <- file.path(project_dir,"Data")
fig_dir <- file.path(project_dir,"Figures")
diablo_data_dir <- file.path(project_dir,"DIABLO","data")
fig4_dir <- file.path(fig_dir,"Figure4")
figs6_dir <- file.path(fig_dir,"FigureS6")

dir.create(fig4_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(figs6_dir,showWarnings=FALSE,recursive=TRUE)

gene_sets <- readRDS(file.path(data_dir,"chaussabel_modules_17.rds"))
module_order <- names(gene_sets)
print(module_order)

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each=2),c("AM","PM")))[-c(9,11,13)]
cols_symptoms <- c("#C75A5A","#5E86C2")
cols_symptoms_named <- setNames(cols_symptoms,c("low","high"))
timepoint_labs <- c("M1"="-1am","1AM"="1am","1PM"="1pm","2AM"="2am","2PM"="2pm","3AM"="3am","3PM"="3pm","5AM"="5am","7AM"="7am","14AM"="14am")

ssgsea <- readRDS(file.path(data_dir,"ssgsea.rds"))

prep_ssgsea_site <- function(ssgsea_df,site_filter){
  df <- ssgsea_df %>%
    dplyr::filter(site==site_filter)
  if(!"sx_group" %in% names(df)){
    if("sx_score_cat" %in% names(df)){
      df$sx_group <- df$sx_score_cat
    }else{
      stop("Neither sx_group nor sx_score_cat found in ssgsea object.")
    }
  }
  if(!"infected" %in% names(df)){
    if("infection_status" %in% names(df)){
      df$infected <- df$infection_status
    }else{
      stop("Neither infected nor infection_status found in ssgsea object.")
    }
  }
  df %>%
    dplyr::mutate(
      timepoint=factor(timepoint,levels=all_timepoint),
      infected=factor(as.character(infected),levels=c("uninfected","infected")),
      sx_group=factor(as.character(sx_group),levels=c("low","high")),
      sex=factor(sex)
    )
}

ssgsea_long_mts <- prep_ssgsea_site(ssgsea,"MTS")
ssgsea_long_throat <- prep_ssgsea_site(ssgsea,"Throat")

# Define helper functions

compute_mean_ssgsea <- function(ssgsea_long_df,filter_expr){
  ssgsea_long_df %>%
    dplyr::filter({{filter_expr}}) %>%
    dplyr::group_by(module,timepoint) %>%
    dplyr::summarise(
      avg_ssgsea_score=mean(ssgsea_scaled,na.rm=TRUE),
      .groups="drop"
    )
}

build_heatmap_matrix <- function(combined_df){
  if(is.null(combined_df) || nrow(combined_df)==0){
    stop("combined_df is empty or NULL in build_heatmap_matrix")
  }
  mat <- combined_df %>%
    dplyr::select(module,timepoint,avg_ssgsea_score) %>%
    tidyr::pivot_wider(names_from=timepoint,values_from=avg_ssgsea_score) %>%
    tibble::column_to_rownames("module") %>%
    as.matrix()
  common_modules <- intersect(module_order,rownames(mat))
  mat <- mat[common_modules,,drop=FALSE]
  if(nrow(mat)==0){
    stop("No modules in module_order matched rownames(mat). Check names(gene_sets) and module_order.")
  }
  if(!"M1" %in% colnames(mat)){
    stop("M1 not found in columns for baseline")
  }
  baseline_M1 <- mat[,"M1"]
  mat <- sweep(mat,1,baseline_M1,"-")
  mat <- sweep(mat,1,baseline_M1,"/")
  mat <- mat * 100
  mat <- mat[,colnames(mat)!="M1",drop=FALSE]
  colnames(mat) <- tolower(colnames(mat))
  lab_row <- rownames(mat)
  lab_row <- gsub("_"," ",lab_row)
  lab_row <- gsub("cell[s]?$","cells",lab_row,ignore.case=TRUE)
  list(mat=mat,lab_row=lab_row)
}

build_radar_data <- function(combined_low,combined_high,low_lab="Low",high_lab="High"){
  combined_low$group <- low_lab
  combined_high$group <- high_lab
  combined_df_group <- dplyr::bind_rows(combined_low,combined_high)
  combined_df_group$timepoint <- factor(combined_df_group$timepoint,levels=all_timepoint)
  combined_df_group$group <- factor(combined_df_group$group,levels=c(low_lab,high_lab))
  m1_df <- combined_df_group %>%
    dplyr::filter(timepoint=="M1") %>%
    dplyr::select(group,module,avg_ssgsea_score) %>%
    tidyr::pivot_wider(names_from=module,values_from=avg_ssgsea_score)
  low_row <- which(m1_df$group==low_lab)
  high_row <- which(m1_df$group==high_lab)
  if(length(low_row)!=1 || length(high_row)!=1){
    stop("Could not uniquely identify Low/High rows for radar_data")
  }
  m1_mat <- as.data.frame(m1_df[,setdiff(colnames(m1_df),"group"),drop=FALSE])
  common_modules <- intersect(module_order,colnames(m1_mat))
  if(length(common_modules)==0){
    stop("No modules in module_order matched columns in radar data. Check names(gene_sets), module_order, and module_label_map.")
  }
  m1_mat <- m1_mat[,common_modules,drop=FALSE]
  plot_data <- rbind(
    m1_mat[low_row,,drop=FALSE],
    m1_mat[high_row,,drop=FALSE]
  )
  rownames(plot_data) <- c(low_lab,high_lab)
  colnames(plot_data) <- colnames(plot_data)
  keep_cols <- colSums(!is.na(plot_data)) > 0
  plot_data <- plot_data[,keep_cols,drop=FALSE]
  max_val <- 1.0
  min_val <- 0.0
  radar_mat <- rbind(
    max=rep(max_val,ncol(plot_data)),
    min=rep(min_val,ncol(plot_data)),
    plot_data
  )
  list(
    radar_data=as.data.frame(radar_mat),
    max_val=max_val,
    min_val=min_val
  )
}

plot_radar_generic <- function(
    radar_obj,
    group_labels,
    colors_border=cols_symptoms,
    module_radius=1.13,
    module_cex=0.60,
    legend_cex=0.60,
    numeric_cex=0.55,
    legend_pos=c(0.70,1.35),
    module_label_colors=NULL,
    module_label_fontface=NULL,
    module_label_box=NULL
){
  op <- par(bg="white",mar=c(0,2.5,1.5,0),xaxs="i",yaxs="i")
  on.exit(par(op))
  plot.new()
  square_side <- min(par("plt")[2]-par("plt")[1],par("plt")[4]-par("plt")[3])
  par(plt=c(0.5-square_side/2,0.5+square_side/2,0.5-square_side/2,0.5+square_side/2))
  rd <- as.data.frame(radar_obj$radar_data)
  if(nrow(rd)<4) stop("radar_data must have rows: max,min,group1,group2")
  n_var <- ncol(rd)
  full_mods <- colnames(rd)
  short_mods <- sapply(strsplit(full_mods,"_"),`[`,1)
  axis_values <- c(0.0,0.2,0.4,0.6,0.8,1.0)
  seg <- length(axis_values)-1
  axis_radii <- c(0.12,0.29,0.45,0.70,0.87,1.04)
  radarchart(
    rd,
    axistype=1,
    vlabels=rep("",n_var),
    vlcex=0,
    caxislabels=rep("",seg+1),
    seg=seg,
    pcol=colors_border,
    plwd=2,
    plty=1,
    cglcol="grey30",
    cglty=2,
    cglwd=1,
    axislabcol="grey30",
    title=""
  )
  for(i in seq_along(axis_values)){
    v <- axis_values[i]
    label <- if(abs(v)<1e-9) "0" else sprintf("%.1f",v)
    x <- axis_radii[i] * cos(pi/2+0.18)
    y <- axis_radii[i] * sin(pi/2+0.18)
    text(x,y,labels=label,cex=numeric_cex,col="grey30",adj=c(0.5,0.5))
  }
  angles <- seq(0,2*pi,length.out=n_var+1)[1:n_var] + pi/2
  for(i in seq_len(n_var)){
    mod_full <- full_mods[i]
    col_i <- if(!is.null(module_label_colors) && mod_full %in% names(module_label_colors)){
      module_label_colors[mod_full]
    }else{
      "black"
    }
    font_i <- if(!is.null(module_label_fontface) && mod_full %in% names(module_label_fontface)){
      module_label_fontface[mod_full]
    }else{
      1L
    }
    x <- module_radius * cos(angles[i])
    y <- module_radius * sin(angles[i])
    box_i <- !is.null(module_label_box) &&
      mod_full %in% names(module_label_box) &&
      isTRUE(module_label_box[[mod_full]])
    if(box_i){
      rect(
        xleft=x - 0.110,
        xright=x + 0.110,
        ybottom=y - 0.040,
        ytop=y + 0.040,
        border="black",
        lwd=0.8
      )
    }
    text(
      x,
      y,
      labels=short_mods[i],
      cex=module_cex,
      col=col_i,
      font=font_i
    )
  }
  legend(
    x=legend_pos[1],
    y=legend_pos[2],
    legend=group_labels,
    col=colors_border,
    lwd=2,
    bty="n",
    cex=legend_cex,
    y.intersp=2.0
  )
}

plot_heatmap_generic <- function(
    mat,
    lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=TRUE,
    show_legend=TRUE,
    row_fontsize=7,
    col_fontsize=7,
    main_title=NULL,
    stars_mat=NULL,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
){
  color_fun <- colorRamp2(c(heat_min,heat_mid,heat_max),c("blue","white","red"))
  ht <- Heatmap(
    mat,
    name=legend_title,
    col=color_fun,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    rect_gp=gpar(col="grey60",lwd=0.6),
    na_col="gray",
    show_row_names=show_row_names,
    row_labels=lab_row,
    row_names_gp=gpar(fontsize=row_fontsize),
    show_column_names=TRUE,
    column_names_gp=gpar(fontsize=col_fontsize),
    column_names_rot=45,
    column_title=main_title,
    column_title_gp=gpar(fontsize=7.8,fontface="bold",vjust=0.2),
    heatmap_legend_param=list(
      legend_height=unit(6,"mm"),
      legend_width=unit(0.5,"mm"),
      title=legend_title,
      at=heat_breaks,
      labels=paste0(heat_breaks,"%"),
      labels_gp=gpar(fontsize=6.5),
      title_gp=gpar(fontsize=7,fontface="bold")
    ),
    cell_fun=function(j,i,x,y,w,h,fill){
      if(!is.null(stars_mat)){
        rn <- rownames(mat)[i]
        cn <- colnames(mat)[j]
        lab <- stars_mat[rn,cn]
        if(!is.na(lab) && nzchar(lab)){
          grid.text(lab,x,y,gp=gpar(col="grey20",fontsize=8,fontface="bold"))
        }
      }
    }
  )
  draw(
    ht,
    heatmap_legend_side="right",
    show_heatmap_legend=show_legend,
    padding=unit(c(1,5,3,3),"pt")
  )
}

fit_module_lmm <- function(df){
  lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID),data=df)
}

add_module_fdr <- function(df,group_vars,p_col="p",q_col="q"){
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate("{q_col}" := p.adjust(.data[[p_col]],method="BH")) %>%
    dplyr::ungroup()
}

run_module_mixed_models <- function(ssgsea_long){
  df_inf <- ssgsea_long %>%
    dplyr::filter(infected=="infected",sx_group %in% c("low","high")) %>%
    dplyr::mutate(sx_group=factor(as.character(sx_group),levels=c("low","high")))
  module_nest <- df_inf %>%
    dplyr::group_by(module) %>%
    tidyr::nest()
  module_nest <- module_nest %>%
    dplyr::mutate(
      model=purrr::map(data,~fit_module_lmm(.x)),
      emm_time_by_sx=purrr::map(model,~emmeans(.x,~timepoint | sx_group)),
      emm_sx_by_time=purrr::map(model,~emmeans(.x,~sx_group | timepoint))
    )
  within_contrasts <- module_nest %>%
    dplyr::mutate(
      within_vs_M1=purrr::map(
        emm_time_by_sx,
        ~{
          ct <- contrast(.x,method="trt.vs.ctrl",ref="M1")
          out <- summary(ct,infer=TRUE) %>% as.data.frame()
          out$contrast_type <- "within_vs_baseline"
          out
        }
      )
    ) %>%
    dplyr::select(module,within_vs_M1) %>%
    tidyr::unnest(within_vs_M1)
  between_contrasts <- module_nest %>%
    dplyr::mutate(
      between_sx=purrr::map(
        emm_sx_by_time,
        ~{
          ct <- contrast(.x,method="revpairwise")
          out <- summary(ct,infer=TRUE) %>% as.data.frame()
          out$contrast_type <- "between_groups"
          out
        }
      )
    ) %>%
    dplyr::select(module,between_sx) %>%
    tidyr::unnest(between_sx)
  all_contrasts <- dplyr::bind_rows(within_contrasts,between_contrasts) %>%
    dplyr::rename(
      estimate=estimate,
      se=SE,
      df=df,
      t=t.ratio,
      p=p.value
    ) %>%
    dplyr::mutate(
      term_clean=dplyr::case_when(
        contrast_type=="within_vs_baseline" ~ contrast,
        contrast_type=="between_groups" ~ paste(timepoint,contrast,sep=":")
      )
    )
  all_contrasts
}

make_star_matrix <- function(sig_df,heat_obj,group_flag,alpha=c(0.20,0.05,0.01)){
  heat_mat <- heat_obj$mat
  modules <- rownames(heat_mat)
  tps <- colnames(heat_mat)
  res <- sig_df %>%
    dplyr::filter(sx_group==group_flag) %>%
    dplyr::mutate(
      module=trimws(as.character(module)),
      timepoint_delta=trimws(tolower(timepoint_delta))
    ) %>%
    dplyr::filter(
      module %in% modules,
      timepoint_delta %in% tps
    ) %>%
    dplyr::mutate(
      star=dplyr::case_when(
        q < alpha[3] ~ "***",
        q < alpha[2] ~ "**",
        q < alpha[1] ~ "*",
        TRUE ~ ""
      )
    )
  stars <- matrix(
    "",
    nrow=nrow(heat_mat),
    ncol=ncol(heat_mat),
    dimnames=dimnames(heat_mat)
  )
  if(nrow(res)>0){
    for(k in seq_len(nrow(res))){
      m <- res$module[k]
      tp <- res$timepoint_delta[k]
      if(!is.na(m) && !is.na(tp)) stars[m,tp] <- res$star[k]
    }
  }
  stars
}

plot_radar_mts_sx <- function(){
  plot_radar_generic(
    radar_obj=radar_mts_obj,
    group_labels=c("Low symptom burden","High symptom burden"),
    colors_border=cols_symptoms,
    legend_pos=c(0.35,1.33),
    module_label_colors=mts_label_style$col,
    module_label_fontface=mts_label_style$font,
    module_label_box=mts_label_style$box
  )
}

plot_heat_mts_low <- function(show_row_names=FALSE,show_legend=FALSE){
  plot_heatmap_generic(
    mat=heat_mts_low_obj$mat,
    lab_row=heat_mts_low_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="Low symptom burden",
    stars_mat=stars_mts_low,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_mts_high <- function(show_row_names=TRUE,show_legend=TRUE){
  plot_heatmap_generic(
    mat=heat_mts_high_obj$mat,
    lab_row=heat_mts_high_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="High symptom burden",
    stars_mat=stars_mts_high,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_radar_throat_sx <- function(){
  plot_radar_generic(
    radar_obj=radar_throat_obj,
    group_labels=c("Low symptom burden","High symptom burden"),
    colors_border=cols_symptoms,
    legend_pos=c(0.35,1.33),
    module_label_colors=throat_label_style$col,
    module_label_fontface=throat_label_style$font,
    module_label_box=throat_label_style$box
  )
}

plot_heat_throat_low <- function(show_row_names=FALSE,show_legend=FALSE){
  plot_heatmap_generic(
    mat=heat_throat_low_obj$mat,
    lab_row=heat_throat_low_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="Low symptom burden",
    stars_mat=stars_throat_low,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_throat_high <- function(show_row_names=TRUE,show_legend=TRUE){
  plot_heatmap_generic(
    mat=heat_throat_high_obj$mat,
    lab_row=heat_throat_high_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="High symptom burden",
    stars_mat=stars_throat_high,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

extract_baseline_lmm_results <- function(all_module_results,site_label){
  all_module_results %>%
    dplyr::filter(
      contrast_type=="between_groups",
      timepoint=="M1"
    ) %>%
    dplyr::mutate(
      site=site_label,
      contrast="High - Low",
      timepoint=as.character(timepoint)
    ) %>%
    add_module_fdr(group_vars=c("site","timepoint","contrast")) %>%
    dplyr::select(site,module,timepoint,contrast,estimate,se,df,t,p,q) %>%
    dplyr::arrange(q)
}

build_module_label_style <- function(
    baseline_results,
    bold_q = 0.20,
    box_q = 0.10
){
  modules <- as.character(module_order)
  label_box  <- setNames(rep(FALSE, length(modules)), modules)
  label_cols <- setNames(rep("black", length(modules)), modules)
  label_font <- setNames(rep(1L, length(modules)), modules)
  tmp <- baseline_results %>%
    dplyr::mutate(module = as.character(module))
  for(k in seq_len(nrow(tmp))){
    m <- tmp$module[k]
    q <- tmp$q[k]
    if(is.na(q) || !(m %in% modules)) next
    if(q < box_q){
      label_cols[m] <- "black"
      label_font[m] <- 2L
      label_box[m]  <- TRUE
    }else if(q < bold_q){
      label_cols[m] <- "black"
      label_font[m] <- 2L
    }
  }
  list(col = label_cols, font = label_font, box = label_box)
}

plot_module_trajectory <- function(traj_df,mod,show_legend=FALSE,show_y_title=FALSE){
  df_plot <- traj_df %>% dplyr::filter(module_short==mod)
  if(nrow(df_plot)==0) stop("Module not found in trajectory data: ",mod)
  df_mean <- df_plot %>%
    dplyr::group_by(timepoint,sx_group) %>%
    dplyr::summarise(mean_ssgsea=mean(ssgsea_scaled,na.rm=TRUE),.groups="drop")
  ggplot(df_plot,aes(x=timepoint,y=ssgsea_scaled,color=sx_group)) +
    geom_point(alpha=0.15,size=1.1,position=position_jitter(width=0.35,height=0)) +
    geom_line(data=df_mean,aes(y=mean_ssgsea,group=sx_group),linewidth=0.55) +
    geom_point(data=df_mean,aes(y=mean_ssgsea),size=1.3) +
    scale_color_manual(
      values=cols_symptoms_named,
      breaks=c("low","high"),
      labels=c("Low","High")
    ) +
    scale_x_discrete(limits=names(timepoint_labs),labels=unname(timepoint_labs),drop=FALSE) +
    scale_y_continuous(
      breaks=scales::pretty_breaks(n=4),
      expand=expansion(mult=c(0.08,0.08))
    ) +
    labs(
      title=mod,
      x=NULL,
      y=if(show_y_title) "ssGSEA" else NULL,
      color="Symptom burden"
    ) +
    theme_bw(base_size=8) +
    theme(
      text=element_text(color="black"),
      plot.title=element_text(face="bold",hjust=0.5,size=8,margin=margin(b=1)),
      legend.position=if(show_legend) "right" else "none",
      legend.title=element_text(size=7,face="bold"),
      legend.text=element_text(size=7),
      legend.key.height=unit(0.35,"cm"),
      legend.key.width=unit(0.65,"cm"),
      axis.text.x=element_text(angle=45,hjust=1,size=6.5),
      axis.text.y=element_text(size=6.5),
      axis.title.y=if(show_y_title) element_text(size=7.5,face="bold",margin=margin(r=4)) else element_blank(),
      axis.ticks.y=element_line(linewidth=0.3),
      axis.ticks.x=element_line(linewidth=0.3),
      axis.ticks.length=unit(2.5,"pt"),
      panel.grid=element_blank(),
      plot.margin=margin(8,0,1,4)
    )
}

# Build figure objects


combined_mts_low <- compute_mean_ssgsea(
  ssgsea_long_mts,
  infected=="infected" & sx_group=="low"
)
combined_mts_high <- compute_mean_ssgsea(
  ssgsea_long_mts,
  infected=="infected" & sx_group=="high"
)
combined_throat_low <- compute_mean_ssgsea(
  ssgsea_long_throat,
  infected=="infected" & sx_group=="low"
)
combined_throat_high <- compute_mean_ssgsea(
  ssgsea_long_throat,
  infected=="infected" & sx_group=="high"
)

radar_mts_obj <- build_radar_data(
  combined_mts_low,
  combined_mts_high,
  low_lab="Low",
  high_lab="High"
)
radar_throat_obj <- build_radar_data(
  combined_throat_low,
  combined_throat_high,
  low_lab="Low",
  high_lab="High"
)

heat_mts_low_obj <- build_heatmap_matrix(combined_mts_low)
heat_mts_high_obj <- build_heatmap_matrix(combined_mts_high)
heat_throat_low_obj <- build_heatmap_matrix(combined_throat_low)
heat_throat_high_obj <- build_heatmap_matrix(combined_throat_high)

# Fit linear mixed models
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID))

all_module_results_mts <- run_module_mixed_models(ssgsea_long_mts)
all_module_results_throat <- run_module_mixed_models(ssgsea_long_throat)

# Baseline symptom-burden differences from the longitudinal mixed model
# Contrast: High - Low at M1
# FDR: BH across modules within the baseline timepoint

baseline_mts_results <- extract_baseline_lmm_results(
  all_module_results_mts,
  site_label = "Nasal"
)

baseline_throat_results <- extract_baseline_lmm_results(
  all_module_results_throat,
  site_label = "Throat"
)

mts_label_style <- build_module_label_style(baseline_mts_results)
throat_label_style <- build_module_label_style(baseline_throat_results)

# Fit linear mixed models for change in ssGSEA values from baseline among infected participants
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID))

heatmap_sig_mts <- all_module_results_mts %>%
  dplyr::filter(contrast_type=="within_vs_baseline") %>%
  dplyr::mutate(
    site="Nasal",
    timepoint_delta=sub(" - M1","",contrast)
  ) %>%
  add_module_fdr(group_vars=c("site","sx_group","timepoint_delta"))

heatmap_sig_throat <- all_module_results_throat %>%
  dplyr::filter(contrast_type=="within_vs_baseline") %>%
  dplyr::mutate(
    site="Throat",
    timepoint_delta=sub(" - M1","",contrast)
  ) %>%
  add_module_fdr(group_vars=c("site","sx_group","timepoint_delta"))

stars_mts_low <- make_star_matrix(heatmap_sig_mts,heat_mts_low_obj,"low")
stars_mts_high <- make_star_matrix(heatmap_sig_mts,heat_mts_high_obj,"high")
stars_throat_low <- make_star_matrix(heatmap_sig_throat,heat_throat_low_obj,"low")
stars_throat_high <- make_star_matrix(heatmap_sig_throat,heat_throat_high_obj,"high")

# Longitudinal module trajectories (absolute scaled ssGSEA; infected participants only)

modules_panel_4c <- c("M1.2","M4.6","M5.7","M3.6")
modules_panel_s6c <- c("M1.2","M3.4","M4.13","M5.1")

traj_df_mts <- ssgsea_long_mts %>%
  dplyr::filter(infected=="infected",sx_group %in% c("low","high")) %>%
  dplyr::mutate(
    module_short=stringr::str_extract(as.character(module),"^M[0-9]+\\.[0-9]+"),
    sx_group=factor(as.character(sx_group),levels=c("low","high")),
    timepoint=factor(as.character(timepoint),levels=names(timepoint_labs))
  )

traj_df_throat <- ssgsea_long_throat %>%
  dplyr::filter(infected=="infected",sx_group %in% c("low","high")) %>%
  dplyr::mutate(
    module_short=stringr::str_extract(as.character(module),"^M[0-9]+\\.[0-9]+"),
    sx_group=factor(as.character(sx_group),levels=c("low","high")),
    timepoint=factor(as.character(timepoint),levels=names(timepoint_labs))
  )

traj_df_mts %>%
  dplyr::filter(module_short %in% modules_panel_4c) %>%
  dplyr::group_by(module_short) %>%
  dplyr::summarise(
    min_ssgsea=min(ssgsea_scaled,na.rm=TRUE),
    max_ssgsea=max(ssgsea_scaled,na.rm=TRUE),
    range_ssgsea=max_ssgsea-min_ssgsea,
    n=dplyr::n(),
    .groups="drop"
  ) %>%
  dplyr::arrange(match(module_short,modules_panel_4c)) %>%
  print()

traj_df_throat %>%
  dplyr::filter(module_short %in% modules_panel_s6c) %>%
  dplyr::group_by(module_short) %>%
  dplyr::summarise(
    min_ssgsea=min(ssgsea_scaled,na.rm=TRUE),
    max_ssgsea=max(ssgsea_scaled,na.rm=TRUE),
    range_ssgsea=max_ssgsea-min_ssgsea,
    n=dplyr::n(),
    .groups="drop"
  ) %>%
  dplyr::arrange(match(module_short,modules_panel_s6c)) %>%
  print()

legend_module_4c <- modules_panel_4c[length(modules_panel_4c)]
legend_4c <- cowplot::get_legend(
  plot_module_trajectory(
    traj_df=traj_df_mts,
    mod=legend_module_4c,
    show_legend=TRUE,
    show_y_title=FALSE
  )
)

panel_4c_plots <- purrr::map2(
  modules_panel_4c,
  seq_along(modules_panel_4c),
  ~plot_module_trajectory(
    traj_df=traj_df_mts,
    mod=.x,
    show_y_title=.y==1
  )
)

fig_4c <- plot_grid(
  plotlist=c(panel_4c_plots,list(legend_4c)),
  ncol=length(modules_panel_4c) + 1,
  rel_widths=c(1.25,rep(1.16,3),0.60),
  align="h"
)

fig_4a_radar_raw <- suppressWarnings(as.ggplot(~ plot_radar_mts_sx()))

fig_4a <- ggdraw() +
  draw_plot(
    fig_4a_radar_raw,
    x=0,
    y=-0.02,
    width=1,
    height=1
  )

fig_4b_heat_low <- suppressWarnings(as.ggplot(
  ~ plot_heat_mts_low(show_row_names=FALSE,show_legend=FALSE)
)) + theme(plot.margin=margin(2,1,2,2))

fig_4b_heat_high <- suppressWarnings(as.ggplot(
  ~ plot_heat_mts_high(show_row_names=TRUE,show_legend=TRUE)
)) + theme(plot.margin=margin(2,1,2,0))

fig_4b_heats <- plot_grid(
  fig_4b_heat_low,
  fig_4b_heat_high,
  ncol=2,
  rel_widths=c(0.63,1.28)
)

fig_4b <- ggdraw() +
  draw_plot(fig_4b_heats,x=0.02,y=-0.015,width=0.96,height=1.01)

figure4_top <- plot_grid(
  fig_4a,
  fig_4b,
  ncol=2,
  rel_widths=c(0.9,1.6),
  labels=c("A","B"),
  label_size=8.5
)

figure4 <- plot_grid(
  figure4_top,
  fig_4c,
  ncol=1,
  rel_heights=c(0.95,0.55),
  labels=c("","C"),
  label_size=8.5
)

ggsave(file.path(fig4_dir,"Figure4a.png"),fig_4a,width=3.5,height=3.35,dpi=600,bg="white")
ggsave(file.path(fig4_dir,"Figure4b.png"),fig_4b,width=6,height=3.35,dpi=600,bg="white")
ggsave(file.path(fig4_dir,"Figure4c.png"),fig_4c,width=9.5,height=2.0,dpi=600,bg="white")
ggsave(file.path(fig_dir,"Figure4.png"),figure4,width=9.5,height=5.35,dpi=600,bg="white")

# Assemble Figure S6 panels - throat only

legend_module_s6c <- modules_panel_s6c[length(modules_panel_s6c)]
legend_s6c <- cowplot::get_legend(
  plot_module_trajectory(
    traj_df=traj_df_throat,
    mod=legend_module_s6c,
    show_legend=TRUE,
    show_y_title=FALSE
  )
)

panel_s6c_plots <- purrr::map2(
  modules_panel_s6c,
  seq_along(modules_panel_s6c),
  ~plot_module_trajectory(
    traj_df=traj_df_throat,
    mod=.x,
    show_y_title=.y==1
  )
)

fig_s6c <- plot_grid(
  plotlist=c(panel_s6c_plots,list(legend_s6c)),
  ncol=length(modules_panel_s6c) + 1,
  rel_widths=c(1.25,rep(1.16,3),0.60),
  align="h"
)

fig_s6a_radar_raw <- suppressWarnings(as.ggplot(~ plot_radar_throat_sx()))

fig_s6a <- ggdraw() +
  draw_plot(
    fig_s6a_radar_raw,
    x=0,
    y=-0.02,
    width=1,
    height=1
  )

fig_s6b_heat_low <- suppressWarnings(as.ggplot(
  ~ plot_heat_throat_low(show_row_names=FALSE,show_legend=FALSE)
)) + theme(plot.margin=margin(2,1,2,2))

fig_s6b_heat_high <- suppressWarnings(as.ggplot(
  ~ plot_heat_throat_high(show_row_names=TRUE,show_legend=TRUE)
)) + theme(plot.margin=margin(2,1,2,0))

fig_s6b_heats <- plot_grid(
  fig_s6b_heat_low,
  fig_s6b_heat_high,
  ncol=2,
  rel_widths=c(0.63,1.28)
)

fig_s6b <- ggdraw() +
  draw_plot(
    fig_s6b_heats,
    x=0.02,
    y=-0.015,
    width=0.96,
    height=1.01
  )

figure_s6_top <- plot_grid(
  fig_s6a,
  fig_s6b,
  ncol=2,
  rel_widths=c(0.9,1.6),
  labels=c("A","B"),
  label_size=8.5
)

figure_s6 <- plot_grid(
  figure_s6_top,
  fig_s6c,
  ncol=1,
  rel_heights=c(0.95,0.55),
  labels=c("","C"),
  label_size=8.5
)

ggsave(file.path(figs6_dir,"FigureS6a.png"),fig_s6a,width=3.5,height=3.35,dpi=600,bg="white")
ggsave(file.path(figs6_dir,"FigureS6b.png"),fig_s6b,width=6,height=3.35,dpi=600,bg="white")
ggsave(file.path(figs6_dir,"FigureS6c.png"),fig_s6c,width=9.5,height=2.0,dpi=600,bg="white")
ggsave(file.path(fig_dir,"FigureS6.png"),figure_s6,width=9.5,height=5.35,dpi=600,bg="white")

# Sanity checks for radar plot values (ssGSEA at baseline)

get_radar_values <- function(radar_obj){
  rd <- radar_obj$radar_data
  rd[-c(1,2),,drop=FALSE] %>%
    tibble::rownames_to_column("group") %>%
    tidyr::pivot_longer(
      -group,
      names_to="module",
      values_to="baseline_ssgsea"
    ) %>%
    tidyr::pivot_wider(
      names_from=group,
      values_from=baseline_ssgsea
    ) %>%
    dplyr::arrange(module)
}

cat("\n=== Nasal baseline radar values - symptom burden among infected participants ===\n")
print(get_radar_values(radar_mts_obj),n=Inf)

cat("\n=== Throat baseline radar values - symptom burden among infected participants ===\n")
print(get_radar_values(radar_throat_obj),n=Inf)

# Baseline symptom-burden differences from the longitudinal mixed model
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID))
# Contrast: High - Low at M1; FDR: BH across modules

cat("\n=== Baseline mixed-model contrast - nasal, High vs Low symptom burden among infected participants ===\n")
print(baseline_mts_results,n=Inf)

cat("\n=== Baseline mixed-model contrast - throat, High vs Low symptom burden among infected participants ===\n")
print(baseline_throat_results,n=Inf)

cat("\n=== Baseline mixed-model contrast - nasal, q < 0.20 ===\n")
print(baseline_mts_results %>% dplyr::filter(q < 0.20),n=Inf)

cat("\n=== Baseline mixed-model contrast - throat, q < 0.20 ===\n")
print(baseline_throat_results %>% dplyr::filter(q < 0.20),n=Inf)

# Overall immune profile difference at baseline by symptom burden among infected participants – nasal

set.seed(123)

permanova_mts_df <- ssgsea_long_mts %>%
  dplyr::filter(timepoint=="M1",infected=="infected",sx_group %in% c("low","high")) %>%
  dplyr::select(studyID,sx_group,age,sex,module,ssgsea_scaled) %>%
  tidyr::pivot_wider(
    names_from=module,
    values_from=ssgsea_scaled
  )

immune_mat <- permanova_mts_df %>%
  dplyr::select(-studyID,-sx_group,-age,-sex) %>%
  as.matrix()

permanova_mts_euclidean <- vegan::adonis2(
  immune_mat ~ sx_group + age + sex,
  data=permanova_mts_df,
  method="euclidean",
  permutations=10000,
  by="margin"
)

print(permanova_mts_euclidean)

# Overall immune profile difference at baseline by symptom burden among infected participants – throat

set.seed(123)

permanova_throat_df <- ssgsea_long_throat %>%
  dplyr::filter(timepoint=="M1",infected=="infected",sx_group %in% c("low","high")) %>%
  dplyr::select(studyID,sx_group,age,sex,module,ssgsea_scaled) %>%
  tidyr::pivot_wider(
    names_from=module,
    values_from=ssgsea_scaled
  )

immune_mat_throat <- permanova_throat_df %>%
  dplyr::select(-studyID,-sx_group,-age,-sex) %>%
  as.matrix()

permanova_throat_euclidean <- vegan::adonis2(
  immune_mat_throat ~ sx_group + age + sex,
  data=permanova_throat_df,
  method="euclidean",
  permutations=10000,
  by="margin"
)

print(permanova_throat_euclidean)

# Sanity checks for heatmaps (% ΔssGSEA)

check_heatmap_values <- function(heat_obj, group_label, site_label) {
  mat <- heat_obj$mat
  cat("\n====", site_label, "-", group_label, "====\n")
  rng <- range(mat, na.rm = TRUE)
  min_idx <- which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  max_idx <- which(mat == max(mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  cat("Range (% change from baseline): [", sprintf("%.1f%%", rng[1]), ", ", sprintf("%.1f%%", rng[2]), "]\n", sep = "")
  cat("Minimum: ", sprintf("%.1f%%", min(mat, na.rm = TRUE)), " (", rownames(mat)[min_idx[1]], ", ", colnames(mat)[min_idx[2]], ")\n", sep = "")
  cat("Maximum: ", sprintf("%.1f%%", max(mat, na.rm = TRUE)), " (", rownames(mat)[max_idx[1]], ", ", colnames(mat)[max_idx[2]], ")\n", sep = "")
  cat("Matrix (% change from baseline):\n")
  print(apply(round(mat, 1), c(1, 2), function(x) paste0(x, "%")), quote = FALSE)
  invisible(mat)
}

mts_low_heat_vals <- check_heatmap_values(heat_mts_low_obj, "Low symptom burden", "Nasal")
mts_high_heat_vals <- check_heatmap_values(heat_mts_high_obj, "High symptom burden", "Nasal")
throat_low_heat_vals <- check_heatmap_values(heat_throat_low_obj, "Low symptom burden", "Throat")
throat_high_heat_vals <- check_heatmap_values(heat_throat_high_obj, "High symptom burden", "Throat")

# Linear mixed effects models comparing later timepoints to baseline by symptom burden among infected participants
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID)); FDR: BH

cat("\n=== Nasal within-group changes, q < 0.20 ===\n")
print(
  heatmap_sig_mts %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(sx_group,timepoint_delta,module),
  n=Inf
)

cat("\n=== Throat within-group changes, q < 0.20 ===\n")
print(
  heatmap_sig_throat %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(sx_group,timepoint_delta,module),
  n=Inf
)

# Timepoint-specific symptom-burden differences in scaled ssGSEA
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID))
# Contrast: High - Low at each timepoint

symptom_timepoint_mts <- all_module_results_mts %>%
  dplyr::filter(contrast_type=="between_groups") %>%
  dplyr::mutate(
    site="Nasal",
    timepoint=as.character(timepoint),
    contrast="High - Low"
  ) %>%
  dplyr::arrange(module,timepoint)

symptom_timepoint_throat <- all_module_results_throat %>%
  dplyr::filter(contrast_type=="between_groups") %>%
  dplyr::mutate(
    site="Throat",
    timepoint=as.character(timepoint),
    contrast="High - Low"
  ) %>%
  dplyr::arrange(module,timepoint)

symptom_timepoint_all <- dplyr::bind_rows(
  symptom_timepoint_mts,
  symptom_timepoint_throat
) %>%
  add_module_fdr(group_vars=c("site","timepoint")) %>%
  dplyr::select(site,module,timepoint,contrast,estimate,se,df,t,p,q) %>%
  dplyr::arrange(site,timepoint,module)

cat("\n=== Symptom-burden differences in absolute module expression at each timepoint – q < 0.20 ===\n")
print(
  symptom_timepoint_all %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(site,timepoint,module),
  n=Inf
)

# Timepoint-specific interaction contrasts - differential changes in ssGSEA values from baseline by symptom burden
# Contrast: (High at timepoint - High at M1) - (Low at timepoint - Low at M1)

compute_delta_delta_sx <- function(ssgsea_long_df,site_label){
  df_inf <- ssgsea_long_df %>%
    dplyr::filter(infected=="infected",sx_group %in% c("low","high")) %>%
    dplyr::mutate(
      timepoint=factor(as.character(timepoint),levels=all_timepoint),
      sx_group=factor(as.character(sx_group),levels=c("low","high")),
      sex=factor(sex)
    )
  module_nest <- df_inf %>%
    dplyr::group_by(module) %>%
    tidyr::nest()
  mod_nest <- module_nest %>%
    dplyr::mutate(
      model=purrr::map(data,~lmer(ssgsea_scaled ~ age + sex + timepoint * sx_group + (1 | studyID),data=.x)),
      emm_joint=purrr::map(model,~emmeans(.x,~timepoint * sx_group)),
      dd=purrr::map(
        emm_joint,
        ~{
          emm <- .x
          grid <- as.data.frame(emm@grid)
          tps <- setdiff(all_timepoint,"M1")
          methods_list <- vector("list",length(tps))
          names(methods_list) <- tps
          for(i in seq_along(tps)){
            tp <- tps[i]
            coefs <- rep(0,nrow(grid))
            coefs[grid$timepoint==tp & grid$sx_group=="high"] <- 1
            coefs[grid$timepoint=="M1" & grid$sx_group=="high"] <- -1
            coefs[grid$timepoint==tp & grid$sx_group=="low"] <- -1
            coefs[grid$timepoint=="M1" & grid$sx_group=="low"] <- 1
            methods_list[[i]] <- coefs
          }
          res <- contrast(emm,method=methods_list)
          summary(res,infer=TRUE) %>% as.data.frame()
        }
      )
    )
  mod_nest %>%
    dplyr::select(module,dd) %>%
    tidyr::unnest(dd) %>%
    dplyr::rename(
      timepoint=contrast,
      estimate=estimate,
      se=SE,
      df=df,
      t=t.ratio,
      p=p.value
    ) %>%
    dplyr::mutate(site=site_label) %>%
    dplyr::select(site,module,timepoint,estimate,se,df,t,p) %>%
    add_module_fdr(group_vars=c("site","timepoint"))
}

interaction_delta_mts <- compute_delta_delta_sx(ssgsea_long_mts,"Nasal")
interaction_delta_throat <- compute_delta_delta_sx(ssgsea_long_throat,"Throat")

cat("\n=== Nasal ΔΔ, q < 0.20 ===\n")
print(
  interaction_delta_mts %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(timepoint,module),
  n=Inf
)

cat("\n=== Throat ΔΔ, q < 0.20 ===\n")
print(
  interaction_delta_throat %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(timepoint,module),
  n=Inf
)