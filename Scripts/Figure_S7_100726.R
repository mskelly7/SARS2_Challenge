# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# Figure S7 – nasal and throat immune module expression by viral load among infected participants
# Last updated: 10/07/2026

script_summary <- "This script generates Figure S7, which summarizes baseline and longitudinal immune module expression patterns by viral load
among infected participants in nasal and throat samples. Previously computed Chaussabel ssGSEA module scores generated from TMM-normalized
logCPM RNA-seq data, scaled to a 0-1 range within each anatomical site, are imported and used for all downstream analyses. Analyses are restricted
to infected participants with low or high viral load. Mean module expression profiles are computed separately for low- and high-viral-load
participants across longitudinal timepoints. Radar plots display baseline (day -1am) immune module expression differences between low- and
high-viral-load groups, while heatmaps display within-group percent changes in scaled module expression relative to baseline over time. Baseline
module-specific differences shown in radar plots are estimated from the longitudinal linear mixed-effects model as the high-versus-low viral load
contrast at day -1am, with false discovery rate correction across modules within the baseline timepoint. Overall baseline immune profile
differences between viral-load groups are additionally assessed using Euclidean-distance PERMANOVA performed on scaled ssGSEA values. Longitudinal
changes in module expression are evaluated using linear mixed-effects models of scaled ssGSEA values adjusted for age and sex, incorporating
timepoint, viral-load category, and subject-level random effects. Heatmap star annotations are based on FDR correction performed within each
viral-load group and timepoint across immune modules. The script also performs timepoint-specific comparisons between viral-load groups and
interaction analyses to identify differential longitudinal immune trajectories between low- and high-viral-load participants."
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
figs7_dir <- file.path(fig_dir,"FigureS7")

dir.create(figs7_dir,showWarnings=FALSE,recursive=TRUE)

gene_sets <- readRDS(file.path(data_dir,"chaussabel_modules_17.rds"))
module_order <- names(gene_sets)
print(module_order)

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each=2),c("AM","PM")))[-c(9,11,13)]
cols_viral <- c("#236343","#653875")
cols_viral_named <- setNames(cols_viral,c("low","high"))

ssgsea <- readRDS(file.path(data_dir,"ssgsea.rds"))

prep_ssgsea_site_vl <- function(ssgsea_df,site_filter){
  df <- ssgsea_df %>%
    dplyr::filter(site==site_filter)
  if(!"vl_cat" %in% names(df)){
    if("vl_group" %in% names(df)){
      df$vl_cat <- df$vl_group
    }else if("viral_load_group" %in% names(df)){
      df$vl_cat <- df$viral_load_group
    }else{
      stop("No viral-load grouping variable found. Expected vl_cat, vl_group, or viral_load_group.")
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
      infection_status=infected,
      vl_cat=factor(as.character(vl_cat),levels=c("low","high")),
      vl_group=factor(as.character(vl_cat),levels=c("low","high")),
      sex=factor(sex)
    )
}

ssgsea_long_mts_vl <- prep_ssgsea_site_vl(ssgsea,"MTS")
ssgsea_long_throat_vl <- prep_ssgsea_site_vl(ssgsea,"Throat")

compute_mean_ssgsea_vl <- function(ssgsea_long_df,vl_group_sel){
  ssgsea_long_df %>%
    dplyr::filter(infected=="infected",vl_group==vl_group_sel) %>%
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
  keep_cols <- colSums(!is.na(plot_data)) > 0
  plot_data <- plot_data[,keep_cols,drop=FALSE]
  radar_mat <- rbind(
    max=rep(1,ncol(plot_data)),
    min=rep(0,ncol(plot_data)),
    plot_data
  )
  list(
    radar_data=as.data.frame(radar_mat),
    max_val=1,
    min_val=0
  )
}

plot_radar_generic <- function(
    radar_obj,
    group_labels,
    colors_border=cols_viral,
    module_radius=1.13,
    module_cex=0.60,
    legend_cex=0.60,
    numeric_cex=0.55,
    legend_pos=c(0.50,1.28),
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

fit_module_lmm_vl <- function(df){
  lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID),data=df)
}

add_module_fdr <- function(df,group_vars,p_col="p",q_col="q"){
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate("{q_col}" := p.adjust(.data[[p_col]],method="BH")) %>%
    dplyr::ungroup()
}

run_module_mixed_models_vl <- function(ssgsea_long_df){
  df_inf <- ssgsea_long_df %>%
    dplyr::filter(infected=="infected",vl_group %in% c("low","high")) %>%
    dplyr::mutate(
      vl_group=factor(as.character(vl_group),levels=c("low","high")),
      timepoint=factor(as.character(timepoint),levels=all_timepoint),
      sex=factor(sex)
    )
  module_nest <- df_inf %>%
    dplyr::group_by(module) %>%
    tidyr::nest()
  module_nest <- module_nest %>%
    dplyr::mutate(
      model=purrr::map(data,~fit_module_lmm_vl(.x)),
      emm_time_by_vl=purrr::map(model,~emmeans(.x,~timepoint | vl_group)),
      emm_vl_by_time=purrr::map(model,~emmeans(.x,~vl_group | timepoint))
    )
  within_contrasts <- module_nest %>%
    dplyr::mutate(
      within_vs_M1=purrr::map(
        emm_time_by_vl,
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
      between_vl=purrr::map(
        emm_vl_by_time,
        ~{
          ct <- contrast(.x,method="revpairwise")
          out <- summary(ct,infer=TRUE) %>% as.data.frame()
          out$contrast_type <- "between_groups"
          out
        }
      )
    ) %>%
    dplyr::select(module,between_vl) %>%
    tidyr::unnest(between_vl)
  dplyr::bind_rows(within_contrasts,between_contrasts) %>%
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
}

make_star_matrix_vl <- function(sig_df,heat_obj,group_flag,alpha=c(0.20,0.05,0.01)){
  heat_mat <- heat_obj$mat
  modules <- rownames(heat_mat)
  tps <- colnames(heat_mat)
  res <- sig_df %>%
    dplyr::filter(vl_group==group_flag) %>%
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

compute_delta_delta_vl <- function(ssgsea_long_df,site_label){
  df_inf <- ssgsea_long_df %>%
    dplyr::filter(infected=="infected",vl_group %in% c("low","high")) %>%
    dplyr::mutate(
      timepoint=factor(as.character(timepoint),levels=all_timepoint),
      vl_group=factor(as.character(vl_group),levels=c("low","high")),
      sex=factor(sex)
    )
  module_nest <- df_inf %>%
    dplyr::group_by(module) %>%
    tidyr::nest()
  mod_nest <- module_nest %>%
    dplyr::mutate(
      model=purrr::map(
        data,
        ~lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID),data=.x)
      ),
      dd=purrr::map(
        model,
        ~{
          emm <- emmeans(.x,~timepoint * vl_group)
          grid <- as.data.frame(emm@grid)
          tps <- setdiff(levels(grid$timepoint),"M1")
          methods_list <- vector("list",length(tps))
          names(methods_list) <- tps
          for(i in seq_along(tps)){
            tp <- tps[i]
            coefs <- rep(0,nrow(grid))
            coefs[grid$timepoint==tp & grid$vl_group=="high"] <- 1
            coefs[grid$timepoint=="M1" & grid$vl_group=="high"] <- -1
            coefs[grid$timepoint==tp & grid$vl_group=="low"] <- -1
            coefs[grid$timepoint=="M1" & grid$vl_group=="low"] <- 1
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

build_module_label_style <- function(baseline_results){
  modules <- as.character(module_order)
  label_box <- setNames(rep(FALSE,length(modules)),modules)
  label_cols <- setNames(rep("black",length(modules)),modules)
  label_font <- setNames(rep(1L,length(modules)),modules)
  tmp <- baseline_results %>%
    dplyr::mutate(module=as.character(module))
  for(k in seq_len(nrow(tmp))){
    m <- tmp$module[k]
    q <- tmp$q[k]
    if(is.na(q) || !(m %in% modules)) next
    if(q < 0.10){
      label_cols[m] <- "black"
      label_font[m] <- 2L
      label_box[m] <- TRUE
    }else if(q < 0.20){
      label_cols[m] <- "black"
      label_font[m] <- 2L
    }
  }
  list(col=label_cols,font=label_font,box=label_box)
}

# Build figure objects

combined_mts_low_vl <- compute_mean_ssgsea_vl(ssgsea_long_mts_vl,"low")
combined_mts_high_vl <- compute_mean_ssgsea_vl(ssgsea_long_mts_vl,"high")
combined_throat_low_vl <- compute_mean_ssgsea_vl(ssgsea_long_throat_vl,"low")
combined_throat_high_vl <- compute_mean_ssgsea_vl(ssgsea_long_throat_vl,"high")

radar_mts_vl_obj <- build_radar_data(
  combined_mts_low_vl,
  combined_mts_high_vl,
  low_lab="Low",
  high_lab="High"
)

radar_throat_vl_obj <- build_radar_data(
  combined_throat_low_vl,
  combined_throat_high_vl,
  low_lab="Low",
  high_lab="High"
)

heat_mts_low_vl_obj <- build_heatmap_matrix(combined_mts_low_vl)
heat_mts_high_vl_obj <- build_heatmap_matrix(combined_mts_high_vl)
heat_throat_low_vl_obj <- build_heatmap_matrix(combined_throat_low_vl)
heat_throat_high_vl_obj <- build_heatmap_matrix(combined_throat_high_vl)

# Fit linear mixed models
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID))

all_module_results_mts_vl <- run_module_mixed_models_vl(ssgsea_long_mts_vl)
all_module_results_throat_vl <- run_module_mixed_models_vl(ssgsea_long_throat_vl)

# Baseline viral-load differences from the longitudinal mixed model
# Contrast: High - Low at M1
# FDR: BH across modules within the baseline timepoint

baseline_mts_vl_results <- extract_baseline_lmm_results(
  all_module_results_mts_vl,
  site_label="Nasal"
)

baseline_throat_vl_results <- extract_baseline_lmm_results(
  all_module_results_throat_vl,
  site_label="Throat"
)

mts_label_style_vl <- build_module_label_style(baseline_mts_vl_results)
throat_label_style_vl <- build_module_label_style(baseline_throat_vl_results)

# Fit linear mixed models for change in ssGSEA values from baseline among infected participants
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID))

heatmap_sig_mts_vl <- all_module_results_mts_vl %>%
  dplyr::filter(contrast_type=="within_vs_baseline") %>%
  dplyr::mutate(
    site="Nasal",
    timepoint_delta=sub(" - M1","",contrast)
  ) %>%
  add_module_fdr(group_vars=c("site","vl_group","timepoint_delta"))

heatmap_sig_throat_vl <- all_module_results_throat_vl %>%
  dplyr::filter(contrast_type=="within_vs_baseline") %>%
  dplyr::mutate(
    site="Throat",
    timepoint_delta=sub(" - M1","",contrast)
  ) %>%
  add_module_fdr(group_vars=c("site","vl_group","timepoint_delta"))

stars_mts_low_vl <- make_star_matrix_vl(heatmap_sig_mts_vl,heat_mts_low_vl_obj,"low")
stars_mts_high_vl <- make_star_matrix_vl(heatmap_sig_mts_vl,heat_mts_high_vl_obj,"high")
stars_throat_low_vl <- make_star_matrix_vl(heatmap_sig_throat_vl,heat_throat_low_vl_obj,"low")
stars_throat_high_vl <- make_star_matrix_vl(heatmap_sig_throat_vl,heat_throat_high_vl_obj,"high")

plot_radar_mts_vl <- function(){
  plot_radar_generic(
    radar_obj=radar_mts_vl_obj,
    group_labels=c("Low viral load","High viral load"),
    colors_border=cols_viral,
    module_label_colors=mts_label_style_vl$col,
    module_label_fontface=mts_label_style_vl$font,
    module_label_box=mts_label_style_vl$box
  )
}

plot_radar_throat_vl <- function(){
  plot_radar_generic(
    radar_obj=radar_throat_vl_obj,
    group_labels=c("Low viral load","High viral load"),
    colors_border=cols_viral,
    module_label_colors=throat_label_style_vl$col,
    module_label_fontface=throat_label_style_vl$font,
    module_label_box=throat_label_style_vl$box
  )
}

plot_heat_mts_low_vl <- function(show_row_names=FALSE,show_legend=FALSE){
  plot_heatmap_generic(
    mat=heat_mts_low_vl_obj$mat,
    lab_row=heat_mts_low_vl_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="Low viral load",
    stars_mat=stars_mts_low_vl,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_mts_high_vl <- function(show_row_names=TRUE,show_legend=TRUE){
  plot_heatmap_generic(
    mat=heat_mts_high_vl_obj$mat,
    lab_row=heat_mts_high_vl_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="High viral load",
    stars_mat=stars_mts_high_vl,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_throat_low_vl <- function(show_row_names=FALSE,show_legend=FALSE){
  plot_heatmap_generic(
    mat=heat_throat_low_vl_obj$mat,
    lab_row=heat_throat_low_vl_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="Low viral load",
    stars_mat=stars_throat_low_vl,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_throat_high_vl <- function(show_row_names=TRUE,show_legend=TRUE){
  plot_heatmap_generic(
    mat=heat_throat_high_vl_obj$mat,
    lab_row=heat_throat_high_vl_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="High viral load",
    stars_mat=stars_throat_high_vl,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

fig_s7a_radar_raw <- suppressWarnings(as.ggplot(~plot_radar_mts_vl()))

fig_s7a <- ggdraw() +
  draw_plot(
    fig_s7a_radar_raw,
    x=0,
    y=-0.02,
    width=1,
    height=1
  )

fig_s7b_heat_low <- suppressWarnings(as.ggplot(
  ~plot_heat_mts_low_vl(show_row_names=FALSE,show_legend=FALSE)
)) + theme(plot.margin=margin(2,1,2,2))

fig_s7b_heat_high <- suppressWarnings(as.ggplot(
  ~plot_heat_mts_high_vl(show_row_names=TRUE,show_legend=TRUE)
)) + theme(plot.margin=margin(2,1,2,0))

fig_s7b_heats <- plot_grid(
  fig_s7b_heat_low,
  fig_s7b_heat_high,
  ncol=2,
  rel_widths=c(0.63,1.28)
)

fig_s7b <- ggdraw() +
  draw_plot(
    fig_s7b_heats,
    x=0.02,
    y=-0.015,
    width=0.96,
    height=1.01
  )

fig_s7c_radar_raw <- suppressWarnings(as.ggplot(~plot_radar_throat_vl()))

fig_s7c <- ggdraw() +
  draw_plot(
    fig_s7c_radar_raw,
    x=0,
    y=-0.02,
    width=1,
    height=1
  )

fig_s7d_heat_low <- suppressWarnings(as.ggplot(
  ~plot_heat_throat_low_vl(show_row_names=FALSE,show_legend=FALSE)
)) + theme(plot.margin=margin(2,1,2,2))

fig_s7d_heat_high <- suppressWarnings(as.ggplot(
  ~plot_heat_throat_high_vl(show_row_names=TRUE,show_legend=TRUE)
)) + theme(plot.margin=margin(2,1,2,0))

fig_s7d_heats <- plot_grid(
  fig_s7d_heat_low,
  fig_s7d_heat_high,
  ncol=2,
  rel_widths=c(0.63,1.28)
)

fig_s7d <- ggdraw() +
  draw_plot(
    fig_s7d_heats,
    x=0.02,
    y=-0.015,
    width=0.96,
    height=1.01
  )

figure_s7 <- plot_grid(
  fig_s7a,
  fig_s7b,
  fig_s7c,
  fig_s7d,
  ncol=2,
  rel_widths=c(0.9,1.6),
  rel_heights=c(1,1),
  labels=c("A","B","C","D"),
  label_size=9.5
)

ggsave(file.path(figs7_dir,"FigureS7a.png"),fig_s7a,width=3.5,height=3.35,dpi=600,bg="white")
ggsave(file.path(figs7_dir,"FigureS7b.png"),fig_s7b,width=6.0,height=3.35,dpi=600,bg="white")
ggsave(file.path(figs7_dir,"FigureS7c.png"),fig_s7c,width=3.5,height=3.35,dpi=600,bg="white")
ggsave(file.path(figs7_dir,"FigureS7d.png"),fig_s7d,width=6.0,height=3.35,dpi=600,bg="white")
ggsave(file.path(fig_dir,"FigureS7.png"),figure_s7,width=9.5,height=6.7,dpi=600,bg="white")

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

cat("\n=== Nasal baseline radar values - viral load among infected participants ===\n")
print(get_radar_values(radar_mts_vl_obj),n=Inf)

cat("\n=== Throat baseline radar values - viral load among infected participants ===\n")
print(get_radar_values(radar_throat_vl_obj),n=Inf)

# Baseline viral-load differences from the longitudinal mixed model
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID))
# Contrast: High - Low at M1; FDR: BH across modules

cat("\n=== Baseline mixed-model contrast - nasal, High vs Low viral load among infected participants ===\n")
print(baseline_mts_vl_results,n=Inf)

cat("\n=== Baseline mixed-model contrast - throat, High vs Low viral load among infected participants ===\n")
print(baseline_throat_vl_results,n=Inf)

cat("\n=== Baseline mixed-model contrast - nasal, q < 0.20 ===\n")
print(
  baseline_mts_vl_results %>%
    dplyr::filter(q < 0.20),
  n=Inf
)

cat("\n=== Baseline mixed-model contrast - throat, q < 0.20 ===\n")
print(
  baseline_throat_vl_results %>%
    dplyr::filter(q < 0.20),
  n=Inf
)

# Overall immune profile difference at baseline by viral load among infected participants – nasal

set.seed(123)

permanova_mts_vl_df <- ssgsea_long_mts_vl %>%
  dplyr::filter(timepoint=="M1",infected=="infected",vl_group %in% c("low","high")) %>%
  dplyr::select(studyID,vl_group,age,sex,module,ssgsea_scaled) %>%
  tidyr::pivot_wider(
    names_from=module,
    values_from=ssgsea_scaled
  )

immune_mat_mts_vl <- permanova_mts_vl_df %>%
  dplyr::select(-studyID,-vl_group,-age,-sex) %>%
  as.matrix()

permanova_mts_vl_euclidean <- vegan::adonis2(
  immune_mat_mts_vl ~ vl_group + age + sex,
  data=permanova_mts_vl_df,
  method="euclidean",
  permutations=10000,
  by="margin"
)

print(permanova_mts_vl_euclidean)

# Overall immune profile difference at baseline by viral load among infected participants – throat

set.seed(123)

permanova_throat_vl_df <- ssgsea_long_throat_vl %>%
  dplyr::filter(timepoint=="M1",infected=="infected",vl_group %in% c("low","high")) %>%
  dplyr::select(studyID,vl_group,age,sex,module,ssgsea_scaled) %>%
  tidyr::pivot_wider(
    names_from=module,
    values_from=ssgsea_scaled
  )

immune_mat_throat_vl <- permanova_throat_vl_df %>%
  dplyr::select(-studyID,-vl_group,-age,-sex) %>%
  as.matrix()

permanova_throat_vl_euclidean <- vegan::adonis2(
  immune_mat_throat_vl ~ vl_group + age + sex,
  data=permanova_throat_vl_df,
  method="euclidean",
  permutations=10000,
  by="margin"
)

print(permanova_throat_vl_euclidean)

# Sanity checks for heatmaps (% ΔssGSEA)

check_heatmap_values <- function(heat_obj,label){
  mat <- heat_obj$mat
  cat("\n====",label,"====\n")
  rng <- range(mat,na.rm=TRUE)
  min_idx <- which(mat == min(mat,na.rm=TRUE),arr.ind=TRUE)[1,]
  max_idx <- which(mat == max(mat,na.rm=TRUE),arr.ind=TRUE)[1,]
  cat("Range (% change from baseline): [",sprintf("%.1f%%",rng[1]),", ",sprintf("%.1f%%",rng[2]),"]\n",sep="")
  cat("Minimum: ",sprintf("%.1f%%",min(mat,na.rm=TRUE))," (",rownames(mat)[min_idx[1]],", ",colnames(mat)[min_idx[2]],")\n",sep="")
  cat("Maximum: ",sprintf("%.1f%%",max(mat,na.rm=TRUE))," (",rownames(mat)[max_idx[1]],", ",colnames(mat)[max_idx[2]],")\n",sep="")
  cat("Matrix (% change from baseline):\n")
  print(apply(round(mat,1),c(1,2),function(x) paste0(x,"%")),quote=FALSE)
  invisible(mat)
}

mts_low_vl_heat_vals <- check_heatmap_values(heat_mts_low_vl_obj,"Nasal – Low viral load")
mts_high_vl_heat_vals <- check_heatmap_values(heat_mts_high_vl_obj,"Nasal – High viral load")
throat_low_vl_heat_vals <- check_heatmap_values(heat_throat_low_vl_obj,"Throat – Low viral load")
throat_high_vl_heat_vals <- check_heatmap_values(heat_throat_high_vl_obj,"Throat – High viral load")

# Linear mixed effects models comparing later timepoints to baseline by viral load among infected participants
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID)); FDR: BH

cat("\n=== Nasal within-group changes, q < 0.20 ===\n")
print(
  heatmap_sig_mts_vl %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(vl_group,timepoint_delta,module),
  n=Inf
)

cat("\n=== Throat within-group changes, q < 0.20 ===\n")
print(
  heatmap_sig_throat_vl %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(vl_group,timepoint_delta,module),
  n=Inf
)

# Timepoint-specific viral-load differences in scaled ssGSEA
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * vl_group + (1 | studyID))
# Contrast: High - Low at each timepoint

viral_timepoint_mts <- all_module_results_mts_vl %>%
  dplyr::filter(contrast_type=="between_groups") %>%
  dplyr::mutate(
    site="Nasal",
    timepoint=as.character(timepoint),
    contrast="High - Low"
  ) %>%
  dplyr::arrange(module,timepoint)

viral_timepoint_throat <- all_module_results_throat_vl %>%
  dplyr::filter(contrast_type=="between_groups") %>%
  dplyr::mutate(
    site="Throat",
    timepoint=as.character(timepoint),
    contrast="High - Low"
  ) %>%
  dplyr::arrange(module,timepoint)

viral_timepoint_all <- dplyr::bind_rows(
  viral_timepoint_mts,
  viral_timepoint_throat
) %>%
  add_module_fdr(group_vars=c("site","timepoint")) %>%
  dplyr::select(site,module,timepoint,contrast,estimate,se,df,t,p,q) %>%
  dplyr::arrange(site,timepoint,module)

cat("\n=== Viral-load differences at each timepoint – q < 0.20 ===\n")
print(
  viral_timepoint_all %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(site,timepoint,module),
  n=Inf
)

# Timepoint-specific interaction contrasts - differential changes in ssGSEA values from baseline by viral load
# Contrast: (High at timepoint - High at M1) - (Low at timepoint - Low at M1)

interaction_delta_mts_vl <- compute_delta_delta_vl(ssgsea_long_mts_vl,"Nasal")
interaction_delta_throat_vl <- compute_delta_delta_vl(ssgsea_long_throat_vl,"Throat")

cat("\n=== Nasal VL ΔΔ, q < 0.20 ===\n")
print(
  interaction_delta_mts_vl %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(timepoint,module),
  n=Inf
)

cat("\n=== Throat VL ΔΔ, q < 0.20 ===\n")
print(
  interaction_delta_throat_vl %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(timepoint,module),
  n=Inf
)