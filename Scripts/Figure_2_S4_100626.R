# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# Figure 2 – nasal microbiome composition by infection status
# Figure S4 – throat microbiome composition by infection status
# Last updated: 10/06/2026

script_summary <- "This script generates Figures 2 and S4 summarizing longitudinal nasal and throat microbiome composition by infection status.
Species-level nasal (V1-V3 16S rRNA) and throat (metagenomic) sequencing data were used to create baseline stacked bar plots and longitudinal
heatmaps. Baseline microbiome composition at the day -1am timepoint was summarized separately for infected and uninfected groups using
relative abundance bar plots. Longitudinal microbiome count data were converted to relative abundances and transformed using a centered log-ratio 
(CLR) transformation following addition of a pseudocount offset. Species-specific linear mixed-effects models of the form CLR_abundance ~ age + 
sex + (day + PM) * infected + (1 | studyID) were fit for each retained microbial species in nasal and throat datasets. Estimated marginal means 
with Satterthwaite denominator degrees of freedom were used to estimate timepoint-specific changes from baseline within infected and uninfected 
groups, PM effects, and infection-associated longitudinal interaction effects. Model estimates and FDR-adjusted significance values were reshaped 
into matrices for visualization as heatmaps showing temporal microbiome changes separately for infected and uninfected participants. Additional 
downstream analyses evaluated baseline differences in microbiome composition by infection status using linear regression models (CLR_abundance ~ 
infected + age + sex) and tested for infection-associated differential longitudinal changes in microbial abundances using timepoint-specific 
delta-delta contrasts derived from mixed-effects models. Baseline community-level microbiome composition was additionally evaluated
using Bray-Curtis PERMANOVA adjusted for age and sex."
cat(script_summary)

remove(list=ls())
version
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)
library(mixOmics)
library(lme4)
library(lmerTest)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplotify)
library(emmeans)
library(scales)
library(vegan)

set.seed(123)

project_dir <- "________________________"
data_dir <- file.path(project_dir,"Data")
fig_dir <- file.path(project_dir,"Figures")
fig2_dir <- file.path(fig_dir,"Figure2")
figs4_dir <- file.path(fig_dir,"FigureS4")

dir.create(fig2_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(figs4_dir,showWarnings=FALSE,recursive=TRUE)

colors_32 <- c(
  "gold1","#6A3D9A","#E31A1C","dodgerblue2","darkturquoise","palegreen2","green4","skyblue2",
  "#FF7F00","#CAB2D6","#FDBF6F","gray50","deeppink1","blue1","steelblue4","maroon",
  "darkorange4","brown","#FB9A99","gray80","#1B9E77","#D95F02","#7570B3","#E7298A",
  "#66A61E","#E6AB02","#A6761D","#666666","#A6CEE3","#B2DF8A","#FB8072","#80B1D3"
)

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each=2),c("am","pm")))[-c(9,11,13)]
timepoint_labs <- c(
  M1="-1am","1AM"="1am","1PM"="1pm","2AM"="2am","2PM"="2pm","3AM"="3am",
  "3PM"="3pm","5AM"="5am","7AM"="7am","14AM"="14am"
)
cols_infection_named <- c(uninfected="#388fcc",infected="#f68c34")

species_panel_2c <- c("Corynebacterium accolens","Corynebacterium macginleyi","Dolosigranulum pigrum","Staphylococcus aureus")
species_panel_2c_labs <- setNames(c("C. accolens","C. macginleyi","D. pigrum","S. aureus"),species_panel_2c)
species_panel_s4c <- c("Neisseria subflava","Prevotella melaninogenica","Streptococcus mitis","Veillonella atypica")
species_panel_s4c_labs <- setNames(c("N. subflava","P. melaninogenica","S. mitis","V. atypica"),species_panel_s4c)

# Define helper functions

force_white <- function(p){
  p + theme(
    plot.background=element_rect(fill="white",colour=NA),
    panel.background=element_rect(fill="white",colour=NA)
  )
}

safe_offset <- function(x){
  nonzero <- x[x>0]
  if(length(nonzero)>0) min(nonzero)/2 else 1e-6
}

make_species_labels_md <- function(species_levels){
  species_levels <- as.character(species_levels)
  labs <- character(length(species_levels))
  names(labs) <- species_levels
  prev_genus <- NA_character_
  for(i in seq_along(species_levels)){
    s <- species_levels[i]
    if(s=="Other"){
      labs[i] <- "Other"
      prev_genus <- NA_character_
      next
    }
    s <- gsub("HMT-([0-9]+)","HMT\\1",s)
    parts <- strsplit(s," ")[[1]]
    genus <- if(length(parts)>=1) parts[1] else NA_character_
    is_hmt_sp <- grepl("\\bsp\\.\\s+HMT[0-9]+$",s)
    genus_display <- if(!is_hmt_sp && !is.na(prev_genus) && genus==prev_genus) paste0(substr(genus,1,1),".") else genus
    prev_genus <- genus
    if(grepl("\\bsp\\.$",s)){
      labs[i] <- paste0("<i>",genus_display,"</i> sp.")
    }else if(grepl("\\bsp\\.\\s+HMT[0-9]+$",s)){
      strain <- sub("^.*sp\\.\\s+","",s)
      labs[i] <- paste0("<i>",genus_display,"</i> sp. ",strain)
    }else if(length(parts)>=2){
      labs[i] <- paste0("<i>",genus_display," ",parts[2],"</i>")
    }else{
      labs[i] <- s
    }
  }
  labs
}

md_to_expr <- function(x){
  if(!grepl("<i>",x,fixed=TRUE)) return(shQuote(x))
  italic_part <- sub(".*<i>([^<]+)</i>.*","\\1",x)
  rest <- gsub(".*</i>\\s*","",x)
  rest <- trimws(rest)
  if(rest=="" || rest==italic_part){
    sprintf('italic("%s")',italic_part)
  }else{
    sprintf('italic("%s")~" %s"',italic_part,rest)
  }
}

make_heatmap_labels <- function(md_labels){
  expr_strings <- vapply(md_labels,md_to_expr,character(1))
  out <- parse(text=expr_strings)
  names(out) <- names(md_labels)
  out
}

fix_species_labels <- function(labels_md){
  if("Peptoniphilus sp.HMT187" %in% names(labels_md)){
    labels_md["Peptoniphilus sp.HMT187"] <- "<i>Peptoniphilus</i> sp. HMT187"
  }
  if("Fusobacterium necrophorum" %in% names(labels_md)){
    labels_md["Fusobacterium necrophorum"] <- "<i>Fusobacterium necrophorum</i>"
  }
  labels_md
}

add_day_pm <- function(phy){
  meta <- data.frame(sample_data(phy))
  meta$day <- NA
  meta$day[meta$timepoint=="M1"] <- "D0"
  meta$day[meta$timepoint %in% c("1AM","1PM")] <- "D1"
  meta$day[meta$timepoint %in% c("2AM","2PM")] <- "D2"
  meta$day[meta$timepoint %in% c("3AM","3PM")] <- "D3"
  meta$day[meta$timepoint=="5AM"] <- "D5"
  meta$day[meta$timepoint=="7AM"] <- "D7"
  meta$day[meta$timepoint=="14AM"] <- "D14"
  meta$day <- factor(meta$day,levels=paste0("D",c(0:3,5,7,14)))
  meta$PM <- "N"
  meta$PM[grep("PM",meta$timepoint)] <- "Y"
  meta$PM <- factor(meta$PM,levels=c("N","Y"))
  sample_data(phy) <- meta
  phy
}

prep_baseline_bar_data <- function(phy){
  phy_m1 <- subset_samples(phy,timepoint=="M1")
  phy_m1_merged <- merge_samples(phy_m1,group="infected")
  phy_m1_rel <- transform_sample_counts(phy_m1_merged,function(x) x/sum(x))
  tax_df <- as.data.frame(tax_table(phy_m1_rel))
  tax_df$Species_full <- trimws(as.character(tax_df$Species))
  tax_df$Species_full[is.na(tax_df$Species_full) | tax_df$Species_full %in% c("NA","")] <- "Other"
  tax_table(phy_m1_rel) <- tax_table(as.matrix(tax_df))
  df <- psmelt(phy_m1_rel) %>%
    dplyr::mutate(
      Species_full=tax_df[OTU,"Species_full"],
      infection_status=ifelse(as.character(Sample)=="N","uninfected","infected")
    )
  vec_spec <- sort(unique(df$Species_full))
  vec_spec <- c(setdiff(vec_spec,"Other"),"Other")
  labels_md <- fix_species_labels(make_species_labels_md(vec_spec))
  labels_hm <- make_heatmap_labels(labels_md)
  labels_bar <- vapply(labels_md,md_to_expr,character(1))
  colors <- setNames(rep(colors_32,length.out=length(vec_spec)),vec_spec)
  if("Other" %in% names(colors)) colors["Other"] <- "grey20"
  bar_df <- df %>%
    dplyr::mutate(
      infection_status=factor(infection_status,levels=c("uninfected","infected")),
      Species_plot=factor(Species_full,levels=vec_spec)
    ) %>%
    dplyr::group_by(infection_status,Species_plot) %>%
    dplyr::summarise(Abundance=sum(Abundance),.groups="drop")
  list(df=df,vec_spec=vec_spec,labels_md=labels_md,labels_hm=labels_hm,labels_bar=labels_bar,colors=colors,bar_df=bar_df)
}

prep_species_clr <- function(phy){
  phy <- add_day_pm(phy)
  count <- as(otu_table(phy),"matrix")
  if(taxa_are_rows(phy)) count <- t(count)
  tax_df <- as.data.frame(tax_table(phy))
  tax_df$Species <- trimws(as.character(tax_df$Species))
  tax_df$Species[is.na(tax_df$Species) | tax_df$Species %in% c("NA","")] <- "Other"
  prop <- count/rowSums(count)
  prop <- as.matrix(as.data.frame(prop))
  offset <- safe_offset(prop)
  clr <- logratio.transfo(prop,logratio="CLR",offset=offset)
  meta <- data.frame(sample_data(phy))
  dat <- cbind(meta,clr[rownames(meta),,drop=FALSE])
  list(phy=phy,count=count,tax_df=tax_df,clr=clr,meta=meta,dat=dat,species_cols=colnames(clr))
}

prep_relative_species_df <- function(phy,species_keep){
  transform_sample_counts(phy,function(x) x/sum(x)) %>%
    psmelt() %>%
    dplyr::mutate(
      Species=trimws(as.character(Species)),
      Species=ifelse(is.na(Species) | Species %in% c("NA","") ,"Other",Species),
      infected=ifelse(infected=="N","uninfected","infected"),
      infected=factor(infected,levels=c("uninfected","infected")),
      timepoint=factor(timepoint,levels=names(timepoint_labs))
    ) %>%
    dplyr::filter(Species %in% species_keep)
}

run_eachtime_models <- function(dat,species_cols,fmlr0){
  tp_map <- tibble::tribble(
    ~timepoint,~day,~PM,
    "M1","D0","N",
    "1am","D1","N",
    "1pm","D1","Y",
    "2am","D2","N",
    "2pm","D2","Y",
    "3am","D3","N",
    "3pm","D3","Y",
    "5am","D5","N",
    "7am","D7","N",
    "14am","D14","N"
  )
  build_contrasts <- function(emm_grid,infected_level){
    grid <- as.data.frame(emm_grid@grid)
    make_coef <- function(day_a,pm_a,inf_a,day_b,pm_b,inf_b){
      coefs <- rep(0,nrow(grid))
      coefs[grid$day==day_a & grid$PM==pm_a & grid$infected==inf_a] <- 1
      coefs[grid$day==day_b & grid$PM==pm_b & grid$infected==inf_b] <- -1
      coefs
    }
    out <- vector("list",length(all_timepoint[-1])+1)
    names(out) <- c(paste0(all_timepoint[-1],"_vs_baseline_",ifelse(infected_level=="Y","inf","uninf")),
                    paste0("PM_",ifelse(infected_level=="Y","infected","uninfected")))
    for(i in seq_along(all_timepoint[-1])){
      tp <- all_timepoint[-1][i]
      row <- tp_map[tp_map$timepoint==tp,,drop=FALSE]
      out[[i]] <- make_coef(row$day,row$PM,infected_level,"D0","N",infected_level)
    }
    out[[length(out)]] <- make_coef("D1","Y",infected_level,"D1","N",infected_level)
    out
  }
  res_eachtime <- purrr::map_dfr(species_cols,function(sp){
    fmlr_tmp <- as.formula(paste0("`",sp,"`",fmlr0))
    fit <- lmerTest::lmer(fmlr_tmp,data=dat,control=lmerControl(check.conv.singular="ignore"))
    emm <- emmeans::emmeans(fit,~day * PM * infected,lmer.df="satterthwaite")
    ct_uninf <- emmeans::contrast(emm,method=build_contrasts(emm,"N"))
    ct_inf <- emmeans::contrast(emm,method=build_contrasts(emm,"Y"))
    dplyr::bind_rows(
      summary(ct_uninf,infer=TRUE) %>% as.data.frame() %>% dplyr::mutate(infected="N"),
      summary(ct_inf,infer=TRUE) %>% as.data.frame() %>% dplyr::mutate(infected="Y")
    ) %>%
      dplyr::mutate(
        bacteria=sp,
        timepoint=factor(c(all_timepoint[-1],"PM",all_timepoint[-1],"PM"),levels=c(all_timepoint[-1],"PM"))
      ) %>%
      dplyr::select(bacteria,timepoint,infected,comparison=contrast,estimate,std.err=SE,df,t.ratio,pvalue=p.value)
  })
  pred_names <- unique(res_eachtime$comparison)
  res_eachtime$padj <- NA_real_
  for(ii in seq_along(pred_names)){
    idx <- res_eachtime$comparison==pred_names[ii]
    res_eachtime$padj[idx] <- p.adjust(res_eachtime$pvalue[idx],method="fdr")
  }
  res_eachtime %>%
    dplyr::mutate(log10qval=-sign(estimate)*log10(padj))
}

build_heat_df <- function(res_eachtime,id,species_order=NULL,label_lookup=NULL){
  coef_wide <- res_eachtime %>%
    dplyr::filter(infected==id) %>%
    dplyr::select(timepoint,bacteria,estimate) %>%
    tidyr::pivot_wider(names_from="bacteria",values_from="estimate") %>%
    dplyr::arrange(timepoint) %>%
    as.data.frame()
  rownames(coef_wide) <- coef_wide$timepoint
  coef_wide <- coef_wide[,-1,drop=FALSE]
  qval_wide <- res_eachtime %>%
    dplyr::filter(infected==id) %>%
    dplyr::select(timepoint,bacteria,log10qval) %>%
    tidyr::pivot_wider(names_from="bacteria",values_from="log10qval") %>%
    dplyr::arrange(timepoint) %>%
    as.data.frame()
  rownames(qval_wide) <- qval_wide$timepoint
  qval_wide <- qval_wide[,-1,drop=FALSE]
  heat_df <- dplyr::left_join(
    coef_wide %>% as.data.frame() %>% tibble::rownames_to_column("timepoint") %>% tidyr::pivot_longer(-timepoint,names_to="bacteria",values_to="estimate"),
    qval_wide %>% as.data.frame() %>% tibble::rownames_to_column("timepoint") %>% tidyr::pivot_longer(-timepoint,names_to="bacteria",values_to="log10qval"),
    by=c("timepoint","bacteria")
  ) %>%
    dplyr::mutate(
      star=dplyr::case_when(
        abs(log10qval)>=3 ~ "***",
        abs(log10qval)>=2 ~ "**",
        abs(log10qval)>=-log10(0.10) ~ "*",
        TRUE ~ ""
      ),
      timepoint=factor(timepoint,levels=c(all_timepoint[-1],"PM"))
    )
  if(!is.null(species_order)){
    bact_levels <- species_order[species_order %in% unique(heat_df$bacteria)]
  }else{
    bact_levels <- unique(heat_df$bacteria)
  }
  bact_levels <- rev(bact_levels)
  heat_df$bacteria <- factor(heat_df$bacteria,levels=bact_levels)
  lab_vec <- if(!is.null(label_lookup)) label_lookup[bact_levels] else setNames(bact_levels,bact_levels)
  list(heat_df=heat_df,bact_vec=bact_levels,lab_vec=lab_vec)
}

make_tc_mats <- function(h){
  df <- h$heat_df %>% dplyr::filter(!is.na(bacteria),!is.na(timepoint))
  est <- df %>%
    dplyr::group_by(bacteria,timepoint) %>%
    dplyr::summarise(estimate=mean(estimate,na.rm=TRUE),.groups="drop") %>%
    tidyr::pivot_wider(names_from=timepoint,values_from=estimate) %>%
    tibble::column_to_rownames("bacteria") %>%
    as.matrix()
  star <- df %>%
    dplyr::group_by(bacteria,timepoint) %>%
    dplyr::summarise(star=paste0(unique(star[star!=""]),collapse=""),.groups="drop") %>%
    tidyr::pivot_wider(names_from=timepoint,values_from=star) %>%
    tibble::column_to_rownames("bacteria") %>%
    as.matrix()
  list(est=est,star=star)
}

plot_micro_heat <- function(mat,stars,lab_vec,title,col_fun,legend_at,
                            show_row_names=TRUE,show_legend=TRUE,
                            draw_padding=unit(c(1,5,3,3),"pt"),
                            column_names_centered=FALSE, 
                            row_fontsize=6.5,col_fontsize=6.5,title_fontsize=7.2,
                            star_size=8.0,grid_col="grey60",grid_lwd=0.6,
                            legend_label_size=6.0,legend_title_size=6.8,
                            legend_height=unit(6,"mm"),legend_width=unit(0.3,"mm"),
                            custom_legend=FALSE,legend_x=0.88,legend_y=0.52,
                            custom_legend_height=unit(12,"mm")){
  col_order <- c(all_timepoint[-1],"PM")
  col_order <- col_order[col_order %in% colnames(mat)]
  mat <- mat[,col_order,drop=FALSE]
  stars <- stars[,col_order,drop=FALSE]
  lab_vec <- lab_vec[rownames(mat)]
  ord <- rev(seq_len(nrow(mat)))
  mat <- mat[ord,,drop=FALSE]
  stars <- stars[ord,,drop=FALSE]
  lab_vec <- lab_vec[ord]
  ht <- Heatmap(
    mat,
    name="Estimate",
    col=col_fun,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    rect_gp=gpar(col=grid_col,lwd=grid_lwd),
    na_col="gray90",
    show_row_names=show_row_names,
    row_labels=lab_vec,
    row_names_side="right",
    row_names_gp=gpar(fontsize=row_fontsize),
    show_column_names=TRUE,
    column_names_gp=gpar(fontsize=col_fontsize),
    column_names_rot=45,
    column_names_centered=column_names_centered,
    column_title=title,
    column_title_gp=gpar(fontsize=title_fontsize,fontface="bold"),
    show_heatmap_legend=show_legend && !custom_legend,
    heatmap_legend_param=list(
      legend_height=legend_height,
      legend_width=legend_width,
      title="Estimate",
      at=legend_at,
      labels_gp=gpar(fontsize=legend_label_size),
      title_gp=gpar(fontsize=legend_title_size,fontface="bold")
    ),
    cell_fun=function(j,i,x,y,w,h,fill){
      lab <- stars[i,j]
      if(!is.na(lab) && lab!="") grid.text(lab,x,y,gp=gpar(col="grey20",fontsize=star_size,fontface="bold"))
    }
  )
  draw(ht,heatmap_legend_side="right",show_heatmap_legend=show_legend && !custom_legend,padding=draw_padding)
  if(show_legend && custom_legend){
    lgd <- Legend(
      title="Estimate",
      col_fun=col_fun,
      at=legend_at,
      legend_height=custom_legend_height,
      labels_gp=gpar(fontsize=legend_label_size),
      title_gp=gpar(fontsize=legend_title_size,fontface="bold")
    )
    draw(lgd,x=unit(legend_x,"npc"),y=unit(legend_y,"npc"),just=c("left","center"))
  }
}

build_bar_plot <- function(bar_df,vec_spec,species_colors,species_labels,
                           legend_size=6.7,legend_key_height=0.32,legend_spacing_y=0.03,
                           plot_margin=margin(t=10,r=-5,b=-5,l=5)){
  ggplot(bar_df,aes(x=infection_status,y=Abundance,fill=Species_plot)) +
    geom_col(width=0.95,color="grey70",linewidth=0.1) +
    theme_bw(base_size=8) +
    scale_fill_manual(values=species_colors,breaks=vec_spec,labels=parse(text=species_labels[vec_spec])) +
    coord_cartesian(ylim=c(0,1),expand=0) +
    scale_y_continuous(labels=scales::percent_format(accuracy=1),breaks=seq(0,1,0.2),expand=c(0,0)) +
    labs(x="",y="Relative abundance") +
    guides(fill=guide_legend(title=NULL,ncol=1,byrow=TRUE)) +
    scale_x_discrete(labels=c(uninfected="Uninfected",infected="Infected")) +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.border=element_blank(),
      axis.line=element_line(linewidth=0.3),
      axis.ticks=element_line(linewidth=0.3),
      axis.title.y=element_text(size=7.2,face="bold",margin=margin(r=4)),
      axis.text.x=element_text(size=7.0,face="bold",hjust=0.5,color="black"),
      axis.text.y=element_text(size=6.5,color="black"),
      legend.title=element_blank(),
      legend.text=element_text(size=legend_size),
      legend.key.width=unit(0.25,"cm"),
      legend.key.height=unit(legend_key_height,"cm"),
      legend.spacing.y=unit(legend_spacing_y,"cm"),
      legend.spacing.x=unit(0.02,"cm"),
      legend.box.spacing=unit(0.35,"cm"),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,4,0,0),
      legend.position="right",
      plot.margin=plot_margin
    )
}

plot_species_trajectory <- function(
    sp,
    df_rel,
    panel_labs,
    panel_ylim,
    sigma=0.0001,
    show_legend=FALSE,
    show_y_title=FALSE
){
  df_plot <- df_rel %>% dplyr::filter(Species==sp)
  if(nrow(df_plot)==0) stop("Species not found: ",sp)
  df_mean <- df_plot %>%
    dplyr::group_by(timepoint,infected) %>%
    dplyr::summarise(mean_abundance=mean(Abundance,na.rm=TRUE),.groups="drop")
  ylim_sp <- panel_ylim[[sp]]
  breaks_sp <- c(0,0.0002,0.001,0.01,0.05,0.10,0.25,1)
  if(sp=="Corynebacterium accolens"){
    breaks_sp <- c(0,0.05,0.10,0.25,0.50,1)
  }
  if(sp=="Prevotella melaninogenica"){
    breaks_sp <- c(0,0.01,0.05,0.10,0.25,0.50)
  }
  if(sp %in% c("Neisseria subflava","Streptococcus mitis")){
    breaks_sp <- c(0,0.01,0.05,0.10,0.25,0.50)
  }
  if(sp %in% "Veillonella atypica"){
    breaks_sp <- c(0,0.005,0.01,0.02,0.05,0.10,0.25)
  }
  breaks_sp <- breaks_sp[
    breaks_sp >= ylim_sp[1] &
      breaks_sp <= ylim_sp[2]
  ]
  ggplot(df_plot,aes(x=timepoint,y=Abundance,color=infected)) +
    geom_point(alpha=0.15,size=1.1,position=position_jitter(width=0.35,height=0)) +
    geom_line(data=df_mean,aes(y=mean_abundance,group=infected),linewidth=0.55) +
    geom_point(data=df_mean,aes(y=mean_abundance),size=1.5) +
    scale_color_manual(values=cols_infection_named,breaks=c("uninfected","infected"),labels=c("Uninfected","Infected")) +
    scale_x_discrete(limits=names(timepoint_labs),labels=unname(timepoint_labs),drop=FALSE) +
    scale_y_continuous(
      trans=scales::pseudo_log_trans(
        base=10,
        sigma=sigma
      ),
      breaks=breaks_sp,
      labels=function(x){
        sapply(x,function(v){
          if(v==0) return("0%")
          if(abs(v-0.0002)<1e-12) return("0.02%")
          if(abs(v-0.001)<1e-12) return("0.1%")
          if(abs(v-0.005)<1e-12) return("0.5%")
          if(abs(v-0.01)<1e-12) return("1%")
          if(abs(v-0.02)<1e-12) return("2%")
          if(abs(v-0.05)<1e-12) return("5%")
          if(abs(v-0.10)<1e-12) return("10%")
          if(abs(v-0.25)<1e-12) return("25%")
          if(abs(v-0.50)<1e-12) return("50%")
          if(abs(v-1)<1e-12) return("100%")
          paste0(round(v*100,1),"%")
        })
      },
      limits=ylim_sp,
      expand=expansion(mult=c(0.02,0.04))
    ) +
    labs(title=panel_labs[sp],x=NULL,y=if(show_y_title) "Relative abundance" else NULL,color=NULL) +
    theme_bw(base_size=8) +
    theme(
      plot.title=element_text(face="bold.italic",hjust=0.5,size=7,margin=margin(b=1)),
      legend.position=if(show_legend) "right" else "none",
      legend.text=element_text(size=7),
      legend.key.height=unit(0.35,"cm"),
      legend.key.width=unit(0.65,"cm"),
      axis.text.x=element_text(angle=45,hjust=1,size=6.2),
      axis.text.y=element_text(size=6.0),
      axis.title.y=if(show_y_title) element_text(size=7.2,face="bold",margin=margin(r=4)) else element_blank(),
      axis.ticks.y=element_line(linewidth=0.3),
      axis.ticks.x=element_line(linewidth=0.3),
      axis.ticks.length=unit(2.5,"pt"),
      panel.grid=element_blank(),
      plot.margin=margin(8,1,1,4)
    )
}

build_species_panel <- function(species_vec,df_rel,panel_labs,panel_ylim,panel_sigma=NULL,legend_width=0.50){
  if(is.null(panel_sigma)){
    panel_sigma <- setNames(rep(0.0001,length(species_vec)),species_vec)
  }
  legend_obj <- cowplot::get_legend(
    plot_species_trajectory(species_vec[length(species_vec)],df_rel,panel_labs,panel_ylim,sigma=unname(panel_sigma[species_vec[length(species_vec)]]),
                            show_legend=TRUE,show_y_title=FALSE)
  )
  species_plots <- purrr::map2(
    species_vec,
    seq_along(species_vec),
    ~plot_species_trajectory(.x,df_rel,panel_labs,panel_ylim,sigma=unname(panel_sigma[.x]),show_y_title=.y==1)
  )
  plot_grid(
    plotlist=c(species_plots,list(legend_obj)),
    ncol=length(species_vec)+1,
    rel_widths=c(1.25,rep(1.17,length(species_vec)-1),legend_width),
    align="h"
  )
}

run_baseline_regression <- function(meta,count_mat,vec_spec){
  idx_m1 <- meta$timepoint=="M1"
  meta_m1 <- meta[idx_m1,,drop=FALSE]
  count_m1 <- count_mat[rownames(meta_m1),vec_spec,drop=FALSE]
  prop_m1 <- count_m1/rowSums(count_m1)
  offset <- safe_offset(prop_m1)
  clr_m1 <- logratio.transfo(prop_m1,logratio="CLR",offset=offset)
  dat <- cbind(meta_m1,clr_m1[rownames(meta_m1),,drop=FALSE])
  res <- lapply(colnames(clr_m1),function(sp){
    fml <- as.formula(paste0("`",sp,"` ~ infected + age + sex"))
    fit <- lm(fml,data=dat)
    cf <- coef(summary(fit))
    infected_row <- grep("^infected",rownames(cf),value=TRUE)[1]
    if(is.na(infected_row)) return(data.frame(Species=sp,term=NA_character_,estimate=NA_real_,pvalue=NA_real_))
    data.frame(Species=sp,term=infected_row,estimate=cf[infected_row,"Estimate"],pvalue=cf[infected_row,"Pr(>|t|)"])
  }) %>%
    dplyr::bind_rows()
  res$padj <- p.adjust(res$pvalue,method="fdr")
  res[order(res$pvalue),]
}

extract_sig <- function(df,label_map){
  df %>%
    dplyr::mutate(
      bacteria_chr=as.character(bacteria),
      Species=unname(label_map[bacteria_chr]),
      qval=padj
    ) %>%
    dplyr::filter(qval<0.10) %>%
    dplyr::arrange(qval) %>%
    dplyr::select(Species,timepoint,infected,estimate,pvalue,qval,log10qval)
}

compute_delta_delta_species <- function(df_long,species_cols,site_label){
  df_inf <- df_long %>%
    dplyr::filter(!is.na(infected)) %>%
    dplyr::mutate(inf_group=if(is.factor(infected)) infected else factor(infected))
  
  lvl <- levels(df_inf$inf_group)
  if(length(lvl)!=2) stop("infected must have exactly 2 levels in df_long")
  sp_nest <- df_inf %>%
    tidyr::pivot_longer(
      dplyr::all_of(species_cols),
      names_to="species",
      values_to="clr"
    ) %>%
    dplyr::group_by(species) %>%
    tidyr::nest()
  sp_nest <- sp_nest %>%
    dplyr::mutate(
      model=purrr::map(
        data,
        ~lmer(clr ~ age + sex + timepoint * inf_group + (1 | studyID),data=.x)
      ),
      dd=purrr::map(model,~{
        emm <- emmeans(.x,~timepoint * inf_group,lmer.df="satterthwaite")
        grid <- as.data.frame(emm@grid)
        tps <- setdiff(levels(grid$timepoint),"M1")
        lvl_here <- levels(grid$inf_group)
        ref <- lvl_here[1]
        trg <- lvl_here[2]
        methods_list <- vector("list",length(tps))
        names(methods_list) <- tps
        for(i in seq_along(tps)){
          tp <- tps[i]
          coefs <- rep(0,nrow(grid))
          coefs[grid$timepoint==tp & grid$inf_group==trg] <- 1
          coefs[grid$timepoint=="M1" & grid$inf_group==trg] <- -1
          coefs[grid$timepoint==tp & grid$inf_group==ref] <- -1
          coefs[grid$timepoint=="M1" & grid$inf_group==ref] <- 1
          methods_list[[i]] <- coefs
        }
        
        res <- contrast(emm,method=methods_list)
        out <- summary(res,infer=TRUE) %>% as.data.frame()
        out$timepoint <- names(methods_list)
        out
      })
    )
  sp_nest %>%
    dplyr::select(species,dd) %>%
    tidyr::unnest(dd) %>%
    dplyr::transmute(
      site=site_label,
      species=species,
      timepoint=timepoint,
      estimate=estimate,
      se=SE,
      df=df,
      t=t.ratio,
      p=p.value
    ) %>%
    dplyr::group_by(site,timepoint) %>%
    dplyr::mutate(
      q_timepoint=p.adjust(p,method="BH")
    ) %>%
    dplyr::ungroup()
}

# Data preprocessing

path_mts <- file.path(data_dir,"phy.mts.16s.pruned.species.rds")
path_throat <- file.path(data_dir,"phy.throat.ngs.pruned.species.rds")

phy_mts <- readRDS(path_mts)
phy_throat <- readRDS(path_throat)

mts_base <- prep_baseline_bar_data(phy_mts)
throat_base <- prep_baseline_bar_data(phy_throat)

vec_spec_mts <- mts_base$vec_spec
species_labels_mts <- mts_base$labels_md
species_labels_mts_hm <- mts_base$labels_hm
species_labels_mts_bar <- mts_base$labels_bar
species_colors_mts <- mts_base$colors
mts_bar_df <- mts_base$bar_df

vec_spec_throat <- throat_base$vec_spec
species_labels_throat <- throat_base$labels_md
species_labels_throat_hm <- throat_base$labels_hm
species_labels_throat_bar <- throat_base$labels_bar
species_colors_throat <- throat_base$colors
throat_bar_df <- throat_base$bar_df

fmlr0 <- "~ age + sex + (day + PM) * infected + (1 | studyID)"

mts_long <- prep_species_clr(phy_mts)
throat_long <- prep_species_clr(phy_throat)

res_mts_eachtime <- run_eachtime_models(
  dat=mts_long$dat,
  species_cols=mts_long$species_cols,
  fmlr0=fmlr0
)

res_throat_eachtime <- run_eachtime_models(
  dat=throat_long$dat,
  species_cols=throat_long$species_cols,
  fmlr0=fmlr0
)

mts_rel_df <- prep_relative_species_df(phy_mts,species_panel_2c)
throat_rel_df <- prep_relative_species_df(phy_throat,species_panel_s4c)

# Barplots: Figure 2a and Figure S4a

fig_2a <- force_white(
  build_bar_plot(
    bar_df=mts_bar_df,
    vec_spec=vec_spec_mts,
    species_colors=species_colors_mts,
    species_labels=species_labels_mts_bar,
    legend_size=6.5
  )
)

fig_s4a <- force_white(
  build_bar_plot(
    bar_df=throat_bar_df,
    vec_spec=vec_spec_throat,
    species_colors=species_colors_throat,
    species_labels=species_labels_throat_bar,
    legend_size=6.5,
    legend_key_height=0.15,
    legend_spacing_y=0.025,
    plot_margin=margin(t=10,r=0,b=-5,l=5)
  )
)

# Heatmaps: Figure 2b and Figure S4b

mts_uninf_heat <- build_heat_df(res_eachtime=res_mts_eachtime,id="N",species_order=vec_spec_mts,label_lookup=species_labels_mts)
mts_inf_heat <- build_heat_df(res_eachtime=res_mts_eachtime,id="Y",species_order=vec_spec_mts,label_lookup=species_labels_mts)
throat_uninf_heat <- build_heat_df(res_eachtime=res_throat_eachtime,id="N",species_order=vec_spec_throat,label_lookup=species_labels_throat)
throat_inf_heat <- build_heat_df(res_eachtime=res_throat_eachtime,id="Y",species_order=vec_spec_throat,label_lookup=species_labels_throat)

mts_u <- make_tc_mats(mts_uninf_heat)
mts_i <- make_tc_mats(mts_inf_heat)
thr_u <- make_tc_mats(throat_uninf_heat)
thr_i <- make_tc_mats(throat_inf_heat)

col_fun_mts <- circlize::colorRamp2(
  c(-3,-1.5,0,1.5,3),
  c("#4A1486","#B9A3D5","white","#F2BE8A","#C65A11")
)

col_fun_thr <- circlize::colorRamp2(
  c(-3,-1.5,0,1.5,3),
  c("#4A1486","#B9A3D5","white","#F2BE8A","#C65A11")
)

fig_2b_heat_uninf <- suppressWarnings(
  as.ggplot(~plot_micro_heat(
    mts_u$est,mts_u$star,species_labels_mts_hm,
    title="Uninfected",col_fun=col_fun_mts,legend_at=c(-3,-1.5,0,1.5,3),
    show_row_names=FALSE,show_legend=FALSE
  ))
) + theme(plot.margin=margin(2,1,2,0))

fig_2b_heat_inf <- suppressWarnings(
  as.ggplot(~plot_micro_heat(
    mts_i$est,mts_i$star,species_labels_mts_hm,
    title="Infected",col_fun=col_fun_mts,legend_at=c(-3,-1.5,0,1.5,3),
    show_row_names=TRUE,show_legend=TRUE
  ))
) + theme(plot.margin=margin(2,1,2,0))

fig_2b_heats <- plot_grid(
  fig_2b_heat_uninf,
  fig_2b_heat_inf,
  ncol=2,
  rel_widths=c(0.58,1.20),
  align="hv"
)

fig_2b <- ggdraw() +
  draw_plot(fig_2b_heats,x=0.02,y=-0.025,width=0.96,height=1.025)

fig_s4b_heat_uninf <- suppressWarnings(
  as.ggplot(~plot_micro_heat(
    thr_u$est,thr_u$star,species_labels_throat_hm,
    title="Uninfected",col_fun=col_fun_thr,legend_at=c(-3,-1.5,0,1.5,3),
    show_row_names=FALSE,show_legend=FALSE,
    draw_padding=unit(c(1,1,3,5),"pt"),
    column_names_centered=TRUE,
    star_size=7.0
  ))
) + theme(plot.margin=margin(2,1,2,0))

fig_s4b_heat_inf <- suppressWarnings(
  as.ggplot(~plot_micro_heat(
    thr_i$est,thr_i$star,species_labels_throat_hm,
    title="Infected",col_fun=col_fun_thr,legend_at=c(-3,-1.5,0,1.5,3),
    show_row_names=TRUE,show_legend=TRUE,
    draw_padding=unit(c(1,1,3,5),"pt"),
    column_names_centered=TRUE,
    star_size=7.0,
    custom_legend=TRUE,
    legend_x=0.86,
    legend_y=0.52,
    legend_label_size=6.5,
    legend_title_size=7.0,
    custom_legend_height=unit(12,"mm")
  ))
) + theme(plot.margin=margin(2,1,2,0))

fig_s4b_heats <- plot_grid(
  fig_s4b_heat_uninf,
  fig_s4b_heat_inf,
  ncol=2,
  rel_widths=c(0.65,1.20),
  align="hv"
)

fig_s4b <- ggdraw() +
  draw_plot(fig_s4b_heats,x=0.02,y=-0.025,width=0.96,height=1.025)

# Species trajectories: Figure 2c and Figure S4c

species_panel_2c_ylim <- list(
  "Corynebacterium accolens"=c(0,1),
  "Corynebacterium macginleyi"=c(0,1),
  "Dolosigranulum pigrum"=c(0,1),
  "Staphylococcus aureus"=c(0,1)
)

species_panel_2c_sigma <- c(
  "Corynebacterium accolens"=0.05,
  "Corynebacterium macginleyi"=0.0001,
  "Dolosigranulum pigrum"=0.0001,
  "Staphylococcus aureus"=0.0001
)

fig_2c <- build_species_panel(
  species_vec=species_panel_2c,
  df_rel=mts_rel_df,
  panel_labs=species_panel_2c_labs,
  panel_ylim=species_panel_2c_ylim,
  panel_sigma=species_panel_2c_sigma,
  legend_width=0.50
)

species_panel_s4c_ylim <- list(
  "Prevotella melaninogenica"=c(0,0.5),
  "Streptococcus mitis"=c(0,0.5),
  "Neisseria subflava"=c(0,0.5),
  "Veillonella atypica"=c(0,0.1)
)

species_panel_s4c_sigma <- c(
  "Prevotella melaninogenica" = 0.01,
  "Streptococcus mitis" = 0.01,
  "Neisseria subflava" = 0.01,
  "Veillonella atypica" = 0.002
)

fig_s4c <- build_species_panel(
  species_vec=species_panel_s4c,
  df_rel=throat_rel_df,
  panel_labs=species_panel_s4c_labs,
  panel_ylim=species_panel_s4c_ylim,
  panel_sigma=species_panel_s4c_sigma,
  legend_width=0.50
)

# Assemble final figures

figure2_top <- plot_grid(
  fig_2a,
  fig_2b,
  ncol=2,
  rel_widths=c(1.2,1.6),
  labels=c("A","B"),
  label_size=7.8
)

figure2 <- plot_grid(
  figure2_top,
  fig_2c,
  ncol=1,
  rel_heights=c(0.9,0.55),
  labels=c("","C"),
  label_size=7.8
)

ggsave(file.path(fig2_dir,"Figure2a.png"),fig_2a,width=3.5,height=3.25,dpi=600,bg="white")
ggsave(file.path(fig2_dir,"Figure2b.png"),fig_2b,width=6.0,height=3.25,dpi=600,bg="white")
ggsave(file.path(fig2_dir,"Figure2c.png"),fig_2c,width=9.5,height=1.75,dpi=600,bg="white")
ggsave(file.path(fig_dir,"Figure2.png"),figure2,width=9.5,height=5.0,dpi=600,bg="white")

figureS4_top <- plot_grid(
  fig_s4a,
  fig_s4b,
  ncol=2,
  rel_widths=c(1.25,1.6),
  labels=c("A","B"),
  label_size=8.5
)

figureS4 <- plot_grid(
  figureS4_top,
  fig_s4c,
  ncol=1,
  rel_heights=c(1.05,0.55),
  labels=c("","C"),
  label_size=8.5
)

ggsave(file.path(figs4_dir,"FigureS4a.png"),fig_s4a,width=3.5,height=4.25,dpi=600,bg="white")
ggsave(file.path(figs4_dir,"FigureS4b.png"),fig_s4b,width=6.0,height=4.25,dpi=600,bg="white")
ggsave(file.path(figs4_dir,"FigureS4c.png"),fig_s4c,width=9.5,height=1.65,dpi=600,bg="white")
ggsave(file.path(fig_dir,"FigureS4.png"),figureS4,width=9.5,height=5.90,dpi=600,bg="white")

# PERMANOVA comparing microbiome composition at baseline by infection status

set.seed(123)

mts_m1 <- subset_samples(phy_mts,timepoint=="M1")
mts_rel <- transform_sample_counts(
  mts_m1,
  function(x) x/sum(x)
)
mts_mat <- as(otu_table(mts_rel),"matrix")
if(taxa_are_rows(mts_rel)) mts_mat <- t(mts_mat)
mts_meta <- data.frame(sample_data(mts_rel))
mts_bray <- vegan::vegdist(
  mts_mat,
  method="bray"
)
mts_permanova <- vegan::adonis2(
  mts_bray ~ infected + age + sex,
  data=mts_meta,
  permutations=10000,
  by="margin"
)

cat("\n=== Baseline nasal microbiome composition by infection status (PERMANOVA) ===\n")
print(mts_permanova)

set.seed(123)

throat_m1 <- subset_samples(phy_throat,timepoint=="M1")
throat_rel <- transform_sample_counts(
  throat_m1,
  function(x) x/sum(x)
)
throat_mat <- as(otu_table(throat_rel),"matrix")
if(taxa_are_rows(throat_rel)) throat_mat <- t(throat_mat)
throat_meta <- data.frame(sample_data(throat_rel))
throat_bray <- vegan::vegdist(
  throat_mat,
  method="bray"
)
throat_permanova <- vegan::adonis2(
  throat_bray ~ infected + age + sex,
  data=throat_meta,
  permutations=10000,
  by="margin"
)

cat("\n=== Baseline throat microbiome composition by infection status (PERMANOVA) ===\n")
print(throat_permanova)

# Baseline differences in CLR-transformed species abundances by infection status
# Model: CLR_abundance_species ~ infected + age + sex; FDR: BH

res_mts_baseline <- run_baseline_regression(
  meta=mts_long$meta,
  count_mat=mts_long$count,
  vec_spec=vec_spec_mts
)

res_mts_baseline %>% dplyr::filter(padj < 0.20)
# No differences in relative abundances of retained species in the nasal microbiome at baseline

res_throat_baseline <- run_baseline_regression(
  meta=throat_long$meta,
  count_mat=throat_long$count,
  vec_spec=vec_spec_throat
)

res_throat_baseline %>% dplyr::filter(padj < 0.20)
# No differences in relative abundances of retained species in the throat microbiome at baseline

# Timepoint-specific interaction contrasts - differential changes in CLR-transformed abundances from baseline by infection status

interaction_delta_mts_species <- compute_delta_delta_species(
  df_long=mts_long$dat,
  species_cols=mts_long$species_cols,
  site_label="Nasal"
)

cat("\n=== Nasal species ΔΔ, q_timepoint < 0.20 ===\n")
interaction_delta_mts_species %>%
  dplyr::filter(q_timepoint < 0.20) %>%
  dplyr::arrange(timepoint,species) %>%
  print(n=Inf)

interaction_delta_throat_species <- compute_delta_delta_species(
  df_long=throat_long$dat,
  species_cols=throat_long$species_cols,
  site_label="Throat"
)

cat("\n=== Throat species ΔΔ, q_timepoint < 0.20 ===\n")
interaction_delta_throat_species %>%
  dplyr::filter(q_timepoint < 0.20) %>%
  dplyr::arrange(timepoint,species) %>%
  print(n=Inf)