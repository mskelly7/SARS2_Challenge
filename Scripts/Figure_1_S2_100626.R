# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH 
# Figure 1 - nasal immune module expression by infection status
# Figure S2 - throat immune module expression by infection status
# Last updated: 10/07/2026

script_summary <- "This script generates Figures 1 and S2, which summarize baseline and longitudinal immune module expression patterns by infection
status in nasal and throat samples, respectively. Previously computed Chaussabel ssGSEA module scores generated from TMM-normalized logCPM RNA-seq data, 
scaled to a 0-1 range within each anatomical site, are imported and used for all downstream analyses. Mean module expression profiles are computed
separately for infected and uninfected participants across longitudinal timepoints. Radar plots display baseline (day -1am) immune module expression
differences between participants who subsequently became infected and those who remained uninfected, while heatmaps display within-group percent
changes in scaled module expression relative to baseline over time. Longitudinal trajectory plots are generated for selected immune modules to
visualize temporal patterns in individual participants and group means. Baseline module-specific differences shown in the radar plots are estimated
from longitudinal linear mixed-effects models as the infected-versus-uninfected contrast at day -1am, with false discovery rate correction across
modules within the baseline timepoint. Overall baseline immune profile differences are additionally assessed using Euclidean-distance PERMANOVA
performed on scaled ssGSEA values. Longitudinal changes in module expression are evaluated using linear mixed-effects models of scaled ssGSEA values
adjusted for age and sex and incorporating timepoint, infection status, and participant-level random effects. The script also performs timepoint-specific 
interaction analyses to identify differential longitudinal immune trajectories between infected and uninfected participants."
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
fig1_dir <- file.path(fig_dir,"Figure1")
figs2_dir <- file.path(fig_dir,"FigureS2")

dir.create(fig1_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figs2_dir, showWarnings = FALSE, recursive = TRUE)

gene_sets <- readRDS(file.path(data_dir,"chaussabel_modules_17.rds"))
module_order <- names(gene_sets)
print(module_order)

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each = 2),c("AM","PM")))[-c(9,11,13)]
cols_infection <- c("#388fcc","#f68c34")   # Uninfected, Infected

ssgsea <- readRDS(file.path(data_dir,"ssgsea.rds"))

ssgsea_long_mts <- ssgsea %>%
  dplyr::filter(site=="MTS") %>%
  dplyr::mutate(
    timepoint=factor(timepoint,levels=all_timepoint),
    infected=factor(infected,levels=c("uninfected","infected")),
    sex=factor(sex)
  )

ssgsea_long_throat <- ssgsea %>%
  dplyr::filter(site=="Throat") %>%
  dplyr::mutate(
    timepoint=factor(timepoint,levels=all_timepoint),
    infected=factor(infected,levels=c("uninfected","infected")),
    sex=factor(sex)
  )

# Define helper functions

compute_mean_ssgsea <- function(ssgsea_long_df,filter_expr) {
  ssgsea_long_df %>%
    dplyr::filter({{filter_expr}}) %>%
    dplyr::group_by(module,timepoint) %>%
    dplyr::summarise(
      avg_ssgsea_score=mean(ssgsea_scaled,na.rm=TRUE),
      .groups="drop"
    )
}

build_heatmap_matrix <- function(combined_df) {
  if (is.null(combined_df) || nrow(combined_df) == 0) {
    stop("combined_df is empty or NULL in build_heatmap_matrix")
  }
  mat <- combined_df %>%
    dplyr::select(module,timepoint,avg_ssgsea_score) %>%
    tidyr::pivot_wider(names_from=timepoint,values_from=avg_ssgsea_score) %>%
    tibble::column_to_rownames("module") %>%
    as.matrix()
  common_modules <- intersect(module_order,rownames(mat))
  mat <- mat[common_modules,,drop=FALSE]
  if (nrow(mat) == 0) {
    stop("No modules in module_order matched rownames(mat). Check names(gene_sets) and module_order.")
  }
  if (!"M1" %in% colnames(mat)) {
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

build_radar_data <- function(combined_uninf,combined_inf,uninf_lab="Uninfected",inf_lab="Infected") {
  combined_uninf$group <- uninf_lab
  combined_inf$group <- inf_lab
  combined_df_group <- dplyr::bind_rows(combined_uninf,combined_inf)
  combined_df_group$timepoint <- factor(combined_df_group$timepoint,levels=all_timepoint)
  combined_df_group$group <- factor(combined_df_group$group,levels=c(uninf_lab,inf_lab))
  m1_df <- combined_df_group %>%
    dplyr::filter(timepoint=="M1") %>%
    dplyr::select(group,module,avg_ssgsea_score) %>%
    tidyr::pivot_wider(names_from=module,values_from=avg_ssgsea_score)
  u_row <- which(m1_df$group==uninf_lab)
  i_row <- which(m1_df$group==inf_lab)
  if (length(u_row)!=1 || length(i_row)!=1) {
    stop("Could not uniquely identify Uninfected/Infected rows for radar_data")
  }
  m1_mat <- as.data.frame(m1_df[,setdiff(colnames(m1_df),"group"),drop=FALSE])
  common_modules <- intersect(module_order,colnames(m1_mat))
  if (length(common_modules)==0) {
    stop("No modules in module_order matched columns in radar data. Check names(gene_sets), module_order, and module_label_map.")
  }
  m1_mat <- m1_mat[,common_modules,drop=FALSE]
  plot_data <- rbind(
    m1_mat[u_row,,drop=FALSE],
    m1_mat[i_row,,drop=FALSE]
  )
  rownames(plot_data) <- c(uninf_lab,inf_lab)
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
    colors_border         = cols_infection,
    module_radius         = 1.13,
    module_cex            = 0.60,
    legend_cex            = 0.60,
    numeric_cex           = 0.55,
    legend_pos            = c(0.60,1.28),
    module_label_colors   = NULL,
    module_label_fontface = NULL,
    module_label_box      = NULL
) {
  op <- par(bg="white",mar=c(0,2.5,1.5,0),xaxs="i",yaxs="i")
  plot.new()
  square_side <- min(par("plt")[2]-par("plt")[1],par("plt")[4]-par("plt")[3])
  par(plt = c(0.5-square_side/2,0.5+square_side/2,0.5-square_side/2,0.5+square_side/2))
  rd <- as.data.frame(radar_obj$radar_data)
  if (nrow(rd) < 4) stop("radar_data must have rows: max,min,group1,group2")
  n_var <- ncol(rd)
  full_mods  <- colnames(rd)
  short_mods <- sapply(strsplit(full_mods,"_"),`[`,1)
  axis_values <- c(0.0,0.2,0.4,0.6,0.8,1.0)
  seg <- length(axis_values)-1
  axis_radii <- c(0.12,0.29,0.45,0.70,0.87,1.04)
  radarchart(
    rd,
    axistype    = 1,
    vlabels     = rep("",n_var),
    vlcex       = 0,
    caxislabels = rep("",seg+1),
    seg         = seg,
    pcol        = colors_border,
    plwd        = 2,
    plty        = 1,
    cglcol      = "grey30",
    cglty       = 2,
    cglwd       = 1,
    axislabcol  = "grey30",
    title       = ""
  )
  for (i in seq_along(axis_values)) {
    v <- axis_values[i]
    label <- if (abs(v) < 1e-9) "0" else sprintf("%.1f",v)
    x <- axis_radii[i] * cos(pi/2+0.18)
    y <- axis_radii[i] * sin(pi/2+0.18)
    text(x,y,labels = label,cex = numeric_cex,col = "grey30",adj = c(0.5,0.5))
  }
  angles <- seq(0,2*pi,length.out = n_var+1)[1:n_var] + pi/2
  for (i in seq_len(n_var)) {
    mod_full <- full_mods[i]
    col_i <- if (!is.null(module_label_colors) && mod_full %in% names(module_label_colors)) {
      module_label_colors[mod_full]
    } else {
      "black"
    }
    font_i <- if (!is.null(module_label_fontface) && mod_full %in% names(module_label_fontface)) {
      module_label_fontface[mod_full]
    } else {
      1L
    }
    x <- module_radius * cos(angles[i])
    y <- module_radius * sin(angles[i])
    box_i <- !is.null(module_label_box) &&
      mod_full %in% names(module_label_box) &&
      isTRUE(module_label_box[[mod_full]])
    if (box_i) {
      rect(
        xleft   = x - 0.110,
        xright  = x + 0.110,
        ybottom = y - 0.040,
        ytop    = y + 0.040,
        border  = "black",
        lwd     = 0.8
      )
    }
    text(
      x,
      y,
      labels = short_mods[i],
      cex    = module_cex,
      col    = col_i,
      font   = font_i
    )
  }
  legend(
    x         = legend_pos[1],
    y         = legend_pos[2],
    legend    = group_labels,
    col       = colors_border,
    lwd       = 2,
    bty       = "n",
    cex       = legend_cex,
    y.intersp = 2.0
  )
}

plot_heatmap_generic <- function(
    mat,
    lab_row,
    legend_title   = "% ΔssGSEA",
    show_row_names = TRUE,
    show_legend    = TRUE,
    row_fontsize   = 7,
    col_fontsize   = 7,
    main_title     = NULL,
    stars_mat      = NULL,
    heat_min       = -100,
    heat_mid       = 0,
    heat_max       = 100,
    heat_breaks    = seq(-100,100,by=50)
) {
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
    cell_fun=function(j,i,x,y,w,h,fill) {
      if (!is.null(stars_mat)) {
        rn <- rownames(mat)[i]
        cn <- colnames(mat)[j]
        lab <- stars_mat[rn,cn]
        if (!is.na(lab) && nzchar(lab)) {
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

fit_module_lmm <- function(df) {
  lmer(ssgsea_scaled ~ age + sex + timepoint * infected + (1 | studyID),data = df)
}

add_module_fdr <- function(df, group_vars, p_col = "p", q_col = "q") {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate("{q_col}" := p.adjust(.data[[p_col]], method = "BH")) %>%
    dplyr::ungroup()
}

run_module_mixed_models <- function(ssgsea_long) {
  module_nest <- ssgsea_long %>%
    dplyr::group_by(module) %>%
    tidyr::nest()
  module_nest <- module_nest %>%
    dplyr::mutate(
      model = purrr::map(data, ~fit_module_lmm(.x)),
      emm_time_by_inf = purrr::map(model, ~emmeans(.x, ~timepoint | infected)),
      emm_inf_by_time = purrr::map(model, ~emmeans(.x, ~infected | timepoint))
    )
  within_contrasts <- module_nest %>%
    dplyr::mutate(
      within_vs_M1 = purrr::map(
        emm_time_by_inf,
        ~{
          ct <- contrast(.x, method = "trt.vs.ctrl", ref = "M1")
          out <- summary(ct, infer = TRUE) %>% as.data.frame()
          out$contrast_type <- "within_vs_baseline"
          out
        }
      )
    ) %>%
    dplyr::select(module, within_vs_M1) %>%
    tidyr::unnest(within_vs_M1)
  between_contrasts <- module_nest %>%
    dplyr::mutate(
      between_inf = purrr::map(
        emm_inf_by_time,
        ~{
          ct <- contrast(.x, method = "revpairwise")
          out <- summary(ct, infer = TRUE) %>% as.data.frame()
          out$contrast_type <- "between_groups"
          out
        }
      )
    ) %>%
    dplyr::select(module, between_inf) %>%
    tidyr::unnest(between_inf)
  dplyr::bind_rows(within_contrasts, between_contrasts) %>%
    dplyr::rename(
      estimate = estimate,
      se = SE,
      df = df,
      t = t.ratio,
      p = p.value
    ) %>%
    dplyr::mutate(
      term_clean = dplyr::case_when(
        contrast_type == "within_vs_baseline" ~ contrast,
        contrast_type == "between_groups" ~ paste(timepoint, contrast, sep = ":")
      )
    )
}

make_star_matrix <- function(sig_df, heat_obj, infected_flag, alpha = c(0.20, 0.05, 0.01)) {
  heat_mat <- heat_obj$mat
  modules <- rownames(heat_mat)
  tps <- colnames(heat_mat)
  res <- sig_df %>%
    dplyr::filter(infected == infected_flag) %>%
    dplyr::mutate(
      module = trimws(as.character(module)),
      timepoint_delta = trimws(tolower(timepoint_delta))
    ) %>%
    dplyr::filter(
      module %in% modules,
      timepoint_delta %in% tps
    ) %>%
    dplyr::mutate(
      star = dplyr::case_when(
        q < alpha[3] ~ "***",
        q < alpha[2] ~ "**",
        q < alpha[1] ~ "*",
        TRUE ~ ""
      )
    )
  stars <- matrix(
    "",
    nrow = nrow(heat_mat),
    ncol = ncol(heat_mat),
    dimnames = dimnames(heat_mat)
  )
  if (nrow(res) > 0) {
    for (k in seq_len(nrow(res))) {
      m <- res$module[k]
      tp <- res$timepoint_delta[k]
      if (!is.na(m) && !is.na(tp)) stars[m, tp] <- res$star[k]
    }
  }
  stars
}

compute_delta_delta_inf <- function(ssgsea_long_df, site_label) {
  df_inf <- ssgsea_long_df %>%
    dplyr::filter(infected %in% c("uninfected", "infected")) %>%
    dplyr::mutate(
      inf_group = factor(as.character(infected), levels = c("uninfected", "infected"))
    )
  mod_nest <- df_inf %>%
    dplyr::group_by(module) %>%
    tidyr::nest()
  mod_nest <- mod_nest %>%
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmer(ssgsea_scaled ~ age + sex + timepoint * inf_group + (1 | studyID), data = .x)
      ),
      dd = purrr::map(
        model,
        ~{
          emm <- emmeans(.x, ~timepoint * inf_group)
          grid <- as.data.frame(emm@grid)
          tps <- setdiff(levels(grid$timepoint), "M1")
          methods_list <- vector("list", length(tps))
          names(methods_list) <- tps
          for (i in seq_along(tps)) {
            tp <- tps[i]
            coefs <- rep(0, nrow(grid))
            coefs[grid$timepoint == tp & grid$inf_group == "infected"] <- 1
            coefs[grid$timepoint == "M1" & grid$inf_group == "infected"] <- -1
            coefs[grid$timepoint == tp & grid$inf_group == "uninfected"] <- -1
            coefs[grid$timepoint == "M1" & grid$inf_group == "uninfected"] <- 1
            methods_list[[i]] <- coefs
          }
          res <- contrast(emm, method = methods_list)
          summary(res, infer = TRUE) %>% as.data.frame()
        }
      )
    )
  mod_nest %>%
    dplyr::select(module, dd) %>%
    tidyr::unnest(dd) %>%
    dplyr::rename(
      timepoint = contrast,
      estimate = estimate,
      se = SE,
      df = df,
      t = t.ratio,
      p = p.value
    ) %>%
    dplyr::mutate(site = site_label) %>%
    dplyr::select(site, module, timepoint, estimate, se, df, t, p) %>%
    add_module_fdr(group_vars = c("site", "timepoint"))
}

plot_radar_mts_infstatus <- function() {
  plot_radar_generic(
    radar_obj             = radar_mts_obj,
    group_labels          = c("Uninfected","Infected"),
    colors_border         = cols_infection,
    module_label_colors   = mts_label_style$col,
    module_label_fontface = mts_label_style$font,
    module_label_box      = mts_label_style$box
  )
}

plot_heat_mts_uninf <- function(show_row_names = FALSE,show_legend = FALSE) {
  plot_heatmap_generic(
    mat            = heat_mts_uninf_obj$mat,
    lab_row        = heat_mts_uninf_obj$lab_row,
    legend_title   = "% ΔssGSEA",
    show_row_names = show_row_names,
    show_legend    = show_legend,
    row_fontsize   = 7,
    col_fontsize   = 7,
    main_title     = "Uninfected",
    stars_mat      = stars_mts_uninf,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_mts_inf <- function(show_row_names = TRUE,show_legend = TRUE) {
  plot_heatmap_generic(
    mat            = heat_mts_inf_obj$mat,
    lab_row        = heat_mts_inf_obj$lab_row,
    legend_title   = "% ΔssGSEA",
    show_row_names = show_row_names,
    show_legend    = show_legend,
    row_fontsize   = 7,
    col_fontsize   = 7,
    main_title     = "Infected",
    stars_mat      = stars_mts_inf,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_radar_throat_infstatus <- function() {
  plot_radar_generic(
    radar_obj             = radar_throat_obj,
    group_labels          = c("Uninfected","Infected"),
    colors_border         = cols_infection,
    module_label_colors   = throat_label_style$col,
    module_label_fontface = throat_label_style$font,
    module_label_box      = throat_label_style$box
  )
}

plot_heat_throat_uninf <- function(show_row_names=FALSE,show_legend=FALSE) {
  plot_heatmap_generic(
    mat=heat_throat_uninf_obj$mat,
    lab_row=heat_throat_uninf_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="Uninfected",
    stars_mat=stars_throat_uninf,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

plot_heat_throat_inf <- function(show_row_names=TRUE,show_legend=TRUE) {
  plot_heatmap_generic(
    mat=heat_throat_inf_obj$mat,
    lab_row=heat_throat_inf_obj$lab_row,
    legend_title="% ΔssGSEA",
    show_row_names=show_row_names,
    show_legend=show_legend,
    row_fontsize=7,
    col_fontsize=7,
    main_title="Infected",
    stars_mat=stars_throat_inf,
    heat_min=-100,
    heat_mid=0,
    heat_max=100,
    heat_breaks=seq(-100,100,by=50)
  )
}

extract_baseline_lmm_results <- function(all_module_results, site_label) {
  all_module_results %>%
    dplyr::filter(
      contrast_type == "between_groups",
      timepoint == "M1"
    ) %>%
    dplyr::mutate(
      site = site_label,
      contrast = "Infected - Uninfected",
      timepoint = as.character(timepoint)
    ) %>%
    add_module_fdr(group_vars = c("site", "timepoint", "contrast")) %>%
    dplyr::select(site, module, timepoint, contrast, estimate, se, df, t, p, q) %>%
    dplyr::arrange(q)
}

build_module_label_style <- function(baseline_results) {
  modules <- as.character(module_order)
  label_box <- setNames(rep(FALSE,length(modules)),modules)
  label_cols  <- setNames(rep("black",length(modules)),modules)
  label_font  <- setNames(rep(1L,length(modules)),modules)
  tmp <- baseline_results %>%
    dplyr::mutate(module = as.character(module))
  for (k in seq_len(nrow(tmp))) {
    m <- tmp$module[k]
    q <- tmp$q[k]
    if (is.na(q) || !(m %in% modules)) next
    if (q < 0.10) {
      label_cols[m] <- "black"
      label_font[m] <- 2L
      label_box[m]  <- TRUE
    } else if (q < 0.20) {
      label_cols[m] <- "black"
      label_font[m] <- 2L
    }
  }
  list(col = label_cols,font = label_font,box = label_box)
}

# Build figure objects

combined_mts_uninf <- compute_mean_ssgsea(ssgsea_long_mts,infected=="uninfected")
combined_mts_inf <- compute_mean_ssgsea(ssgsea_long_mts,infected=="infected")

combined_throat_uninf <- compute_mean_ssgsea(ssgsea_long_throat,infected=="uninfected")
combined_throat_inf <- compute_mean_ssgsea(ssgsea_long_throat,infected=="infected")

radar_mts_obj <- build_radar_data(
  combined_mts_uninf,
  combined_mts_inf,
  uninf_lab="Uninfected",
  inf_lab="Infected"
)

radar_throat_obj <- build_radar_data(
  combined_throat_uninf,
  combined_throat_inf,
  uninf_lab="Uninfected",
  inf_lab="Infected"
)

heat_mts_uninf_obj <- build_heatmap_matrix(combined_mts_uninf)
heat_mts_inf_obj <- build_heatmap_matrix(combined_mts_inf)

heat_throat_uninf_obj <- build_heatmap_matrix(combined_throat_uninf)
heat_throat_inf_obj <- build_heatmap_matrix(combined_throat_inf)

# Fit linear mixed models for baseline differences and longitudinal changes
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * infected + (1 | studyID))

all_module_results_mts    <- run_module_mixed_models(ssgsea_long_mts)
all_module_results_throat <- run_module_mixed_models(ssgsea_long_throat)

# Baseline infection-status differences from the longitudinal mixed model
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * infected + (1 | studyID))
# Contrast: Infected - Uninfected at M1
# FDR: BH across modules within the baseline timepoint

baseline_mts_results <- extract_baseline_lmm_results(
  all_module_results_mts,
  site_label = "Nasal"
)

baseline_throat_results <- extract_baseline_lmm_results(
  all_module_results_throat,
  site_label = "Throat"
)

mts_label_style    <- build_module_label_style(baseline_mts_results)
throat_label_style <- build_module_label_style(baseline_throat_results)

heatmap_sig_mts <- all_module_results_mts %>%
  dplyr::filter(contrast_type == "within_vs_baseline") %>%
  dplyr::mutate(
    site = "Nasal",
    timepoint_delta = sub(" - M1", "", contrast)
  ) %>%
  add_module_fdr(group_vars = c("site", "infected", "timepoint_delta"))

heatmap_sig_throat <- all_module_results_throat %>%
  dplyr::filter(contrast_type == "within_vs_baseline") %>%
  dplyr::mutate(
    site = "Throat",
    timepoint_delta = sub(" - M1", "", contrast)
  ) %>%
  add_module_fdr(group_vars = c("site", "infected", "timepoint_delta"))

stars_mts_uninf    <- make_star_matrix(heatmap_sig_mts,heat_mts_uninf_obj,"uninfected")
stars_mts_inf      <- make_star_matrix(heatmap_sig_mts,heat_mts_inf_obj,"infected")
stars_throat_uninf <- make_star_matrix(heatmap_sig_throat,heat_throat_uninf_obj,"uninfected")
stars_throat_inf   <- make_star_matrix(heatmap_sig_throat,heat_throat_inf_obj,"infected")

# Longitudinal module trajectories (absolute ssGSEA)

modules_panel_c <- c("M1.2","M3.4","M4.13","M6.17")
cols_infection_named <- setNames(cols_infection,c("uninfected","infected"))
timepoint_labs <- c("M1"="-1am","1AM"="1am","1PM"="1pm","2AM"="2am","2PM"="2pm","3AM"="3am","3PM"="3pm","5AM"="5am","7AM"="7am","14AM"="14am")

traj_df_mts <- ssgsea_long_mts %>%
  dplyr::mutate(
    module_short=stringr::str_extract(as.character(module),"^M[0-9]+\\.[0-9]+"),
    infected=factor(as.character(infected),levels=c("uninfected","infected")),
    timepoint=factor(as.character(timepoint),levels=names(timepoint_labs))
  )

traj_df_mts %>%
  dplyr::filter(module_short %in% modules_panel_c) %>%
  dplyr::group_by(module_short) %>%
  dplyr::summarise(
    min_ssgsea=min(ssgsea_scaled,na.rm=TRUE),
    max_ssgsea=max(ssgsea_scaled,na.rm=TRUE),
    range_ssgsea=max_ssgsea-min_ssgsea,
    n=dplyr::n(),
    .groups="drop"
  ) %>%
  dplyr::arrange(match(module_short,modules_panel_c)) %>%
  print()

plot_module_trajectory <- function(mod,show_legend=FALSE,show_y_title=FALSE){
  df_plot <- traj_df_mts %>% dplyr::filter(module_short == mod)
  if (nrow(df_plot) == 0) stop("Module not found in traj_df_mts: ",mod)
  df_mean <- df_plot %>%
    dplyr::group_by(timepoint,infected) %>%
    dplyr::summarise(mean_ssgsea=mean(ssgsea_scaled,na.rm=TRUE),.groups="drop")
  ggplot(df_plot,aes(x=timepoint,y=ssgsea_scaled,color=infected)) +
    geom_point(alpha=0.20,size=1.1,position=position_jitter(width=0.35,height=0)) +
    geom_line(data=df_mean,aes(y=mean_ssgsea,group=infected),linewidth=0.4) +
    geom_point(data=df_mean,aes(y=mean_ssgsea),size=1.3) +
    scale_color_manual(values=cols_infection_named,breaks=c("uninfected","infected"),labels=c("Uninfected","Infected")) +
    scale_x_discrete(limits=names(timepoint_labs),labels=unname(timepoint_labs),drop=FALSE) +
    scale_y_continuous(
      breaks=scales::pretty_breaks(n=4),
      expand=expansion(mult=c(0.08,0.08))
    ) +
    labs(title=mod,x=NULL,y=if (show_y_title) "ssGSEA" else NULL,color=NULL) +
    theme_bw(base_size=8) +
    theme(
      text=element_text(color="black"),
      plot.title=element_text(face="bold",hjust=0.5,size=8,margin=margin(b=1)),
      legend.position=if (show_legend) "right" else "none",
      legend.text=element_text(size=7),
      legend.key.height=unit(0.35,"cm"),
      legend.key.width=unit(0.55,"cm"),
      axis.text.x=element_text(angle=45,hjust=1,size=6.5),
      axis.text.y=element_text(size=6.5),
      axis.title.y=if (show_y_title) element_text(size=7.5,face="bold",margin=margin(r=4)) else element_blank(),
      axis.ticks.y=element_line(linewidth=0.3),
      axis.ticks.x=element_line(linewidth=0.3),
      axis.ticks.length=unit(2.5,"pt"),
      panel.grid=element_blank(),
      plot.margin=margin(8,0,1,4)
    )
}

legend_module_c <- modules_panel_c[length(modules_panel_c)]
legend_c <- cowplot::get_legend(plot_module_trajectory(legend_module_c,show_legend=TRUE,show_y_title=FALSE))

panel_c_plots <- purrr::map2(
  modules_panel_c,
  seq_along(modules_panel_c),
  ~plot_module_trajectory(
    .x,
    show_y_title=.y==1
  )
)

fig_1c <- plot_grid(
  plotlist=c(panel_c_plots,list(legend_c)),
  ncol=length(modules_panel_c) + 1,
  rel_widths=c(1.25,rep(1.16,3),0.50),
  align="h"
)

fig_1a_radar_raw <- suppressWarnings(as.ggplot(~ plot_radar_mts_infstatus()))

fig_1a <- ggdraw() +
  draw_plot(
    fig_1a_radar_raw,
    x = 0,
    y = -0.02,
    width = 1,
    height = 1
  )

fig_1b_heat_uninf <- suppressWarnings(as.ggplot(
  ~ plot_heat_mts_uninf(show_row_names=FALSE,show_legend=FALSE)
)) + theme(plot.margin=margin(2,1,2,2))

fig_1b_heat_inf <- suppressWarnings(as.ggplot(
  ~ plot_heat_mts_inf(show_row_names=TRUE,show_legend=TRUE)
)) + theme(plot.margin=margin(2,1,2,0))

fig_1b_heats <- plot_grid(
  fig_1b_heat_uninf,
  fig_1b_heat_inf,
  ncol=2,
  rel_widths=c(0.63,1.28)
)

fig_1b <- ggdraw() +
  draw_plot(fig_1b_heats,x=0.02,y=-0.015,width=0.96,height=1.01)

figure1_top <- plot_grid(
  fig_1a,
  fig_1b,
  ncol=2,
  rel_widths=c(0.9,1.6),
  labels=c("A","B"),
  label_size=8.5
)

figure1 <- plot_grid(
  figure1_top,
  fig_1c,
  ncol=1,
  rel_heights=c(0.95,0.55),
  labels=c("","C"),
  label_size=8.5
)

ggsave(file.path(fig1_dir,"Figure1a.png"),fig_1a,width=3.5,height=3.35,dpi=600)
ggsave(file.path(fig1_dir,"Figure1b.png"),fig_1b,width=6,height=3.35,dpi=600)
ggsave(file.path(fig1_dir,"Figure1c.png"),fig_1c,width=9.5,height=2.0,dpi=600,bg="white")
ggsave(file.path(fig_dir,"Figure1.png"),figure1,width=9.5,height=5.35,dpi=600,bg="white")

# Assemble Figure S2 panels - throat only

fig_s2a_radar_raw <- suppressWarnings(as.ggplot(~ plot_radar_throat_infstatus()))

fig_s2a <- ggdraw() +
  draw_plot(
    fig_s2a_radar_raw,
    x = 0,
    y = -0.02,
    width = 1,
    height = 1
  )

fig_s2b_heat_uninf <- suppressWarnings(as.ggplot(
  ~ plot_heat_throat_uninf(show_row_names = FALSE,show_legend = FALSE)
)) + theme(plot.margin = margin(2,1,2,2))

fig_s2b_heat_inf <- suppressWarnings(as.ggplot(
  ~ plot_heat_throat_inf(show_row_names = TRUE,show_legend = TRUE)
)) + theme(plot.margin = margin(2,1,2,0))

fig_s2b_heats <- plot_grid(
  fig_s2b_heat_uninf,
  fig_s2b_heat_inf,
  ncol = 2,
  rel_widths = c(0.63,1.28)
)

fig_s2b <- ggdraw() +
  draw_plot(
    fig_s2b_heats,
    x = 0.02,
    y = -0.015,
    width = 0.96,
    height = 1.01
  )

# Longitudinal module trajectories in throat

modules_panel_s2c <- c("M1.2","M3.4","M6.17","M3.6")

traj_df_throat <- ssgsea_long_throat %>%
  dplyr::mutate(
    module_short=stringr::str_extract(as.character(module),"^M[0-9]+\\.[0-9]+"),
    infected=factor(as.character(infected),levels=c("uninfected","infected")),
    timepoint=factor(as.character(timepoint),levels=names(timepoint_labs))
  )

plot_module_trajectory_throat <- function(mod,show_legend=FALSE,show_y_title=FALSE){
  df_plot <- traj_df_throat %>% dplyr::filter(module_short == mod)
  if (nrow(df_plot) == 0) stop("Module not found in traj_df_throat: ",mod)
  df_mean <- df_plot %>%
    dplyr::group_by(timepoint,infected) %>%
    dplyr::summarise(mean_ssgsea=mean(ssgsea_scaled,na.rm=TRUE),.groups="drop")
  ggplot(df_plot,aes(x=timepoint,y=ssgsea_scaled,color=infected)) +
    geom_point(alpha=0.15,size=1.1,position=position_jitter(width=0.35,height=0)) +
    geom_line(data=df_mean,aes(y=mean_ssgsea,group=infected),linewidth=0.55) +
    geom_point(data=df_mean,aes(y=mean_ssgsea),size=1.3) +
    scale_color_manual(values=cols_infection_named,breaks=c("uninfected","infected"),labels=c("Uninfected","Infected")) +
    scale_x_discrete(limits=names(timepoint_labs),labels=unname(timepoint_labs),drop=FALSE) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=4),expand=expansion(mult=c(0.08,0.08))) +
    labs(title=mod,x=NULL,y=if (show_y_title) "ssGSEA" else NULL,color=NULL) +
    theme_bw(base_size=8) +
    theme(
      text=element_text(color="black"),
      plot.title=element_text(face="bold",hjust=0.5,size=8,margin=margin(b=1)),
      legend.position=if (show_legend) "right" else "none",
      legend.text=element_text(size=7),
      legend.key.height=unit(0.35,"cm"),
      legend.key.width=unit(0.55,"cm"),
      axis.text.x=element_text(angle=45,hjust=1,size=6.5),
      axis.text.y=element_text(size=6.5),
      axis.title.y=if (show_y_title) element_text(size=7.5,face="bold",margin=margin(r=4)) else element_blank(),
      axis.ticks.y=element_line(linewidth=0.3),
      axis.ticks.x=element_line(linewidth=0.3),
      axis.ticks.length=unit(2.5,"pt"),
      panel.grid=element_blank(),
      plot.margin=margin(8,0,1,4)
    )
}

legend_s2c <- cowplot::get_legend(
  plot_module_trajectory_throat(modules_panel_s2c[length(modules_panel_s2c)],show_legend=TRUE)
)

panel_s2c_plots <- purrr::map2(
  modules_panel_s2c,
  seq_along(modules_panel_s2c),
  ~plot_module_trajectory_throat(.x,show_y_title=.y==1)
)

fig_s2c <- plot_grid(
  plotlist=c(panel_s2c_plots,list(legend_s2c)),
  ncol=length(modules_panel_s2c)+1,
  rel_widths=c(1.25,rep(1.16,3),0.50),
  align="h"
)

figure_s2_top <- plot_grid(
  fig_s2a,
  fig_s2b,
  ncol=2,
  rel_widths=c(0.9,1.6),
  labels=c("A","B"),
  label_size=8.5
)

figure_s2 <- plot_grid(
  figure_s2_top,
  fig_s2c,
  ncol=1,
  rel_heights=c(0.95,0.55),
  labels=c("","C"),
  label_size=8.5
)

ggsave(file.path(figs2_dir,"FigureS2a.png"),fig_s2a,width=3.5,height=3.35,dpi=600)
ggsave(file.path(figs2_dir,"FigureS2b.png"),fig_s2b,width=6,height=3.35,dpi=600)
ggsave(file.path(figs2_dir,"FigureS2c.png"),fig_s2c,width=9.5,height=2.0,dpi=600,bg="white")
ggsave(file.path(fig_dir,"FigureS2.png"),figure_s2,width=9.5,height=5.35,dpi=600,bg="white")

# Sanity checks for radar plot values (ssGSEA at baseline)

get_radar_values <- function(scores_df,label,timepoint_sel="M1") {
  df_m1 <- scores_df %>%
    dplyr::filter(timepoint == timepoint_sel) %>%
    dplyr::group_by(module,infection) %>%
    dplyr::summarise(
      mean_ssgsea = mean(ssgsea_score,na.rm=TRUE),
      sd_ssgsea   = sd(ssgsea_score,na.rm=TRUE),
      n           = dplyr::n(),
      .groups     = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from  = infection,
      values_from = c(mean_ssgsea,sd_ssgsea,n)
    )
  df_m1 <- df_m1 %>%
    dplyr::mutate(
      diff_inf_minus_uninf = mean_ssgsea_Infected - mean_ssgsea_Uninfected
    ) %>%
    dplyr::arrange(factor(module,levels = module_order))
  cat("\n====",label,"(baseline",timepoint_sel,") ====\n")
  print(df_m1)
  invisible(df_m1)
}

mts_scores <- ssgsea_long_mts %>%
  dplyr::transmute(
    module,
    timepoint,
    infection=dplyr::if_else(as.character(infected)=="uninfected","Uninfected","Infected"),
    ssgsea_score=ssgsea_scaled,
    age,
    sex
  )

throat_scores <- ssgsea_long_throat %>%
  dplyr::transmute(
    module,
    timepoint,
    infection=dplyr::if_else(as.character(infected)=="uninfected","Uninfected","Infected"),
    ssgsea_score=ssgsea_scaled,
    age,
    sex
  )

mts_radar_vals    <- get_radar_values(mts_scores,"MTS","M1")
throat_radar_vals <- get_radar_values(throat_scores,"Throat","M1")

# Baseline infection-status differences from the longitudinal mixed model
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * infected + (1 | studyID))
# Contrast: Infected - Uninfected at M1; FDR: BH across modules

print(baseline_mts_results, n = Inf)
print(baseline_throat_results, n = Inf)

# Overall immune profile difference at baseline – nasal (PERMANOVA)

set.seed(123)

permanova_mts_df <- ssgsea_long_mts %>%
  dplyr::filter(timepoint == "M1") %>%
  dplyr::select(studyID,infected,age,sex,module,ssgsea_scaled) %>%
  tidyr::pivot_wider(
    names_from = module,
    values_from = ssgsea_scaled
  )

immune_mat <- permanova_mts_df %>%
  dplyr::select(-studyID,-infected,-age,-sex) %>%
  as.matrix()

permanova_mts_euclidean <- vegan::adonis2(
  immune_mat ~ infected + age + sex,
  data = permanova_mts_df,
  method = "euclidean",
  permutations = 10000,
  by = "margin"
)

print(permanova_mts_euclidean)

# Overall immune profile difference at baseline – throat (PERMANOVA)

set.seed(123)

permanova_throat_df <- ssgsea_long_throat %>%
  dplyr::filter(timepoint == "M1") %>%
  dplyr::select(studyID,infected,age,sex,module,ssgsea_scaled) %>%
  tidyr::pivot_wider(
    names_from = module,
    values_from = ssgsea_scaled
  )

immune_mat_throat <- permanova_throat_df %>%
  dplyr::select(-studyID,-infected,-age,-sex) %>%
  as.matrix()

permanova_throat_euclidean <- vegan::adonis2(
  immune_mat_throat ~ infected + age + sex,
  data = permanova_throat_df,
  method = "euclidean",
  permutations = 10000,
  by = "margin"
)

print(permanova_throat_euclidean)

# Sanity checks for heatmaps (%ΔssGSEA)

check_heatmap_values <- function(heat_obj, label) {
  mat <- heat_obj$mat
  cat("\n====", label, "====\n")
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

mts_uninf_heat_vals <- check_heatmap_values(heat_mts_uninf_obj, "Nasal – Uninfected")
mts_inf_heat_vals <- check_heatmap_values(heat_mts_inf_obj, "Nasal – Infected")
throat_uninf_heat_vals <- check_heatmap_values(heat_throat_uninf_obj, "Throat – Uninfected")
throat_inf_heat_vals <- check_heatmap_values(heat_throat_inf_obj, "Throat – Infected")

# Timepoint-specific infection-status differences in scaled ssGSEA
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * infected + (1 | studyID))
# Contrast: Infected - Uninfected at each timepoint

infection_timepoint_mts <- all_module_results_mts %>%
  dplyr::filter(contrast_type == "between_groups") %>%
  dplyr::mutate(
    site = "Nasal",
    timepoint = as.character(timepoint),
    contrast = "Infected - Uninfected"
  ) %>%
  dplyr::arrange(module, timepoint)

infection_timepoint_throat <- all_module_results_throat %>%
  dplyr::filter(contrast_type == "between_groups") %>%
  dplyr::mutate(
    site = "Throat",
    timepoint = as.character(timepoint),
    contrast = "Infected - Uninfected"
  ) %>%
  dplyr::arrange(module, timepoint)

infection_timepoint_all <- dplyr::bind_rows(
  infection_timepoint_mts,
  infection_timepoint_throat
) %>%
  add_module_fdr(group_vars = c("site", "timepoint")) %>%
  dplyr::select(site, module, timepoint, contrast, estimate, se, df, t, p, q) %>%
  dplyr::arrange(site, timepoint, module)

cat("\n=== Infection-status differences at each timepoint – q < 0.20 ===\n")
print(
  infection_timepoint_all %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(site,timepoint,module),
  n=Inf
)

# Linear mixed effects models comparing module expression at later timepoints to baseline by infection status
# Model: lmer(ssgsea_scaled ~ age + sex + timepoint * infected + (1 | studyID))

cat("\n=== Nasal within-group changes, q < 0.20 ===\n")
print(
  heatmap_sig_mts %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(infected,timepoint_delta,module),
  n=Inf
)

cat("\n=== Throat within-group changes, q < 0.20 ===\n")
print(
  heatmap_sig_throat %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(infected,timepoint_delta,module),
  n=Inf
)

# Timepoint–specific interaction contrasts - differential changes in ssGSEA values from baseline at specific timepoints by infection status

interaction_delta_mts_inf    <- compute_delta_delta_inf(ssgsea_long_mts,"Nasal")
interaction_delta_throat_inf <- compute_delta_delta_inf(ssgsea_long_throat,"Throat")

cat("\n=== Nasal ΔΔ, q < 0.20 ===\n")
print(
  interaction_delta_mts_inf %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(timepoint,module),
  n = Inf
)

cat("\n=== Throat ΔΔ, q < 0.20 ===\n")
print(
  interaction_delta_throat_inf %>%
    dplyr::filter(q < 0.20) %>%
    dplyr::arrange(timepoint,module),
  n = Inf
)