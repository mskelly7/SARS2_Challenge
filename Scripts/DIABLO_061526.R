# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH
# DIABLO multi-omic analyses
# Last updated: 06/15/2026

script_summary <- "This script prepares and analyzes paired microbiome and host transcriptomic module datasets using DIABLO (mixOmics) to identify
multi-omic signatures associated with infection status, symptom severity, and viral load across longitudinal timepoints. Previously computed Chaussabel
ssGSEA module scores, scaled to a 0-1 range within each anatomical site, are reshaped into timepoint-specific matrices for DIABLO modeling. The  
ssGSEA dataset provides both immune-module matrices and outcome metadata, while species-level microbiome phyloseq objects are split into
by-timepoint datasets for nasal and throat samples. Microbiome data are transformed to relative abundance by timepoint, the aggregated Other
category is removed prior to modeling, and the remaining species are CLR-transformed before integration. DIABLO models are fit separately for site,
outcome, and timepoint using a design matrix strength of 0.25. Models are initially fit with up to 3 latent components, and the final
number of components is selected using repeated M-fold cross-validation based on balanced error rate with centroid-distance classification. Sparse
feature selection is tuned using candidate keepX values from 5 to 10 features per block, with 100 repeated tuning iterations. Model tuning and
performance assessment use 5-fold cross-validation with 100 repeated iterations. Final model performance is assessed using leave-one-out validation.
Samples lacking paired microbiome and ssGSEA data are excluded."
cat(script_summary)

remove(list=ls())
version
set.seed(123)

library(mixOmics)
library(phyloseq)
library(tidyverse)

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each=2),c("AM","PM")))[-c(9,11,13)]

project_dir <- "_________________"
data_dir <- file.path(project_dir,"Data")
diablo_data_dir <- file.path(project_dir,"DIABLO","data")
diablo_results_base_dir <- file.path(project_dir,"DIABLO")

dir.create(diablo_data_dir,showWarnings=FALSE,recursive=TRUE)

ssgsea <- readRDS(file.path(data_dir,"ssgsea.rds"))

ssgsea <- ssgsea %>%
  dplyr::mutate(
    timepoint=factor(as.character(timepoint),levels=all_timepoint),
    timepoint_order=as.numeric(timepoint),
    infected=factor(as.character(infected),levels=c("uninfected","infected")),
    sex=factor(sex),
    studyID=factor(studyID)
  )

ssgsea_mts <- ssgsea %>%
  dplyr::filter(site=="MTS")

ssgsea_throat <- ssgsea %>%
  dplyr::filter(site=="Throat")

####################################
# MICROBIOME BY-TIMEPOINT DATASETS #
####################################

phy_mts_16s_pruned <- readRDS(file.path(data_dir,"phy.mts.16s.pruned.species.rds"))
phy_throat_ngs_pruned <- readRDS(file.path(data_dir,"phy.throat.ngs.pruned.species.rds"))

make_microbiome_timepoint_list <- function(phy,infected_level=NULL){
  meta <- data.frame(sample_data(phy))
  if(!is.null(infected_level)){
    if(!"infected"%in%colnames(meta)) stop("infected not found in sample_data")
    keep_inf <- rownames(meta)[meta$infected==infected_level]
    phy <- prune_samples(keep_inf,phy)
  }
  meta <- data.frame(sample_data(phy))
  if(!"timepoint"%in%colnames(meta)) stop("timepoint not found in sample_data")
  out_list <- vector("list",length(all_timepoint))
  names(out_list) <- all_timepoint
  for(time in all_timepoint){
    samps_time <- rownames(meta)[meta$timepoint==time]
    phy_sub <- prune_samples(samps_time,phy)
    if(nsamples(phy_sub)==0){
      out_list[[time]] <- phy_sub
    }else{
      phy.prop.species <- transform_sample_counts(phy_sub,function(x) x/sum(x))
      out_list[[time]] <- phy.prop.species
    }
  }
  return(out_list)
}

micro_mts_all <- make_microbiome_timepoint_list(phy_mts_16s_pruned,infected_level=NULL)
saveRDS(micro_mts_all,file.path(diablo_data_dir,"microbiome_mts_species_by_timepoint.rds"))
micro_throat_all <- make_microbiome_timepoint_list(phy_throat_ngs_pruned,infected_level=NULL)
saveRDS(micro_throat_all,file.path(diablo_data_dir,"microbiome_throat_species_by_timepoint.rds"))
micro_mts_infonly <- make_microbiome_timepoint_list(phy_mts_16s_pruned,infected_level="Y")
saveRDS(micro_mts_infonly,file.path(diablo_data_dir,"microbiome_mts_species_by_timepoint_infected_only.rds"))
micro_throat_infonly <- make_microbiome_timepoint_list(phy_throat_ngs_pruned,infected_level="Y")
saveRDS(micro_throat_infonly,file.path(diablo_data_dir,"microbiome_throat_species_by_timepoint_infected_only.rds"))

#################################
# DIABLO MODELLING BY TIMEPOINT #
#################################

make_ssgsea_timepoint_list <- function(ssgsea_df,infected_only=FALSE){
  df <- ssgsea_df
  if(infected_only){
    df <- df %>% dplyr::filter(infected=="infected")
  }
  out <- vector("list",length(all_timepoint))
  names(out) <- all_timepoint
  for(time in all_timepoint){
    tmp <- df %>%
      dplyr::filter(timepoint==time)
    if(nrow(tmp)==0){
      out[[time]] <- NULL
    }else{
      tmp_metadata <- tmp %>%
        dplyr::select(
          sample,
          timepoint,
          timepoint_order,
          infected,
          sx_score_cat,
          vl_cat
        ) %>%
        dplyr::distinct() %>%
        tibble::column_to_rownames("sample")
      tmp_mat <- tmp %>%
        dplyr::select(sample,module,ssgsea_scaled) %>%
        tidyr::pivot_wider(names_from=module,values_from=ssgsea_scaled) %>%
        tibble::column_to_rownames("sample") %>%
        as.matrix()
      out[[time]] <- list(metadata=tmp_metadata,mat=tmp_mat)
    }
  }
  return(out)
}

run_diablo_for_timepoint <- function(time,
                                     micro_list,
                                     ssgsea_list,
                                     outcome_type = c("infection","symptoms","viral_load"),
                                     site_label = c("MTS","Throat"),
                                     res_dir,
                                     design_strength = 0.1,
                                     use_transposed_otu = FALSE) {
  outcome_type <- match.arg(outcome_type)
  site_label <- match.arg(site_label)
  filt.phy.prop <- micro_list[[time]]
  if (is.null(filt.phy.prop) || nsamples(filt.phy.prop) == 0){
    message("Skipping ",site_label," ",time," ",outcome_type," design=",design_strength,": no microbiome samples")
    return(invisible(NULL))
  }
  # Remove aggregated "Other" category prior to CLR transformation
  if("Species" %in% colnames(tax_table(filt.phy.prop))){
    filt.phy.prop <- subset_taxa(filt.phy.prop,Species!="Other")
  }
  if (ntaxa(filt.phy.prop) == 0){
    message("Skipping ",site_label," ",time," ",outcome_type," design=",design_strength,": no microbiome taxa after removing Other")
    return(invisible(NULL))
  }
  df_asv <- data.frame(tax_table(filt.phy.prop))
  df_asv$asv <- rownames(df_asv)
  otu_raw <- otu_table(filt.phy.prop)
  otu <- if (use_transposed_otu) t(otu_raw) else otu_raw
  if (outcome_type == "infection") {
    clr_otu <- logratio.transfo(otu,logratio="CLR",offset=1)
    meta_cols <- c("timepoint","timepoint_order","infected")
    outcome_col <- "infected"
    outcome_levels <- c("uninfected","infected")
  } else if (outcome_type == "symptoms") {
    prop_otu <- as.matrix(as.data.frame(otu))
    my_offset <- min(prop_otu[prop_otu > 0]) / 2
    clr_otu <- logratio.transfo(otu,logratio="CLR",offset=my_offset)
    meta_cols <- c("timepoint","timepoint_order","sx_score_cat")
    outcome_col <- "sx_score_cat"
    outcome_levels <- c("low","high")
  } else {
    prop_otu <- as.matrix(as.data.frame(otu))
    my_offset <- min(prop_otu[prop_otu > 0]) / 2
    clr_otu <- logratio.transfo(otu,logratio="CLR",offset=my_offset)
    meta_cols <- c("timepoint","timepoint_order","vl_cat")
    outcome_col <- "vl_cat"
    outcome_levels <- c("low","high")
  }
  lookup_vector <- setNames(df_asv$Species,df_asv$asv)
  updated_vector <- unname(lookup_vector[colnames(clr_otu)])
  colnames(clr_otu) <- updated_vector
  class(clr_otu) <- "matrix"
  df_otu_metadata <- data.frame(sample_data(filt.phy.prop))
  keep_cols_otu <- intersect(meta_cols,colnames(df_otu_metadata))
  df_otu_metadata <- df_otu_metadata[,keep_cols_otu,drop=FALSE]
  df_otu_metadata$sample <- rownames(df_otu_metadata)
  this_ssgsea <- ssgsea_list[[time]]
  if(is.null(this_ssgsea)){
    message("Skipping ",site_label," ",time," ",outcome_type," design=",design_strength,": no precomputed ssGSEA data")
    return(invisible(NULL))
  }
  df_ssgsea_metadata <- this_ssgsea$metadata
  b <- this_ssgsea$mat
  if(nrow(df_ssgsea_metadata) == 0 || nrow(b) == 0){
    message("Skipping ",site_label," ",time," ",outcome_type," design=",design_strength,": no precomputed ssGSEA samples")
    return(invisible(NULL))
  }
  a <- clr_otu
  comm_sample_id <- intersect(rownames(a),rownames(b))
  if (length(comm_sample_id) < 3){
    message("Skipping ",site_label," ",time," ",outcome_type," design=",design_strength,": fewer than 3 overlapping samples")
    return(invisible(NULL))
  }
  if (outcome_type == "infection") {
    inf_raw <- df_ssgsea_metadata[comm_sample_id,outcome_col]
    if(all(na.omit(as.character(inf_raw)) %in% c("N","Y"))){
      infection_status <- ifelse(inf_raw == "N","uninfected","infected")
    }else{
      infection_status <- as.character(inf_raw)
    }
    Y <- factor(infection_status,levels=outcome_levels)
  } else {
    var_status <- df_ssgsea_metadata[comm_sample_id,outcome_col]
    Y <- factor(var_status,levels=outcome_levels)
  }
  if(any(table(Y) < 3)){
    message("Skipping ",site_label," ",time," ",outcome_type," design=",design_strength,": fewer than 3 samples in one outcome group")
    return(invisible(NULL))
  }
  a <- a[comm_sample_id,,drop=FALSE]
  b <- b[comm_sample_id,,drop=FALSE]
  data_block <- list(rnaseq=b,microbiome=a)  
  design <- matrix(design_strength,
                   ncol=length(data_block),
                   nrow=length(data_block),
                   dimnames=list(names(data_block),names(data_block)))
  diag(design) <- 0
  design_label <- paste0("design_",gsub("\\.","p",as.character(design_strength)))
  folds_perf <- 5
  nrepeat_perf <- 100
  folds_tune <- 5
  nrepeat_tune <- 100
  basic.diablo.model <- block.splsda(X=data_block,Y=Y,ncomp=3,design=design)
  perf.diablo <- perf(basic.diablo.model,
                      validation="Mfold",
                      folds=folds_perf,
                      nrepeat=nrepeat_perf)
  ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER","centroids.dist"]
  if (ncomp == 1) ncomp <- 2
  n_workers <- 4
  BPPARAM <- BiocParallel::SnowParam(workers=n_workers)
  test.keepX <- list(rnaseq=5:10,
                     microbiome=5:10)
  tune.feature <- tune.block.splsda(X=data_block,
                                    Y=Y,
                                    ncomp=ncomp,
                                    test.keepX=test.keepX,
                                    design=design,
                                    validation="Mfold",
                                    folds=folds_tune,
                                    nrepeat=nrepeat_tune,
                                    dist="centroids.dist",
                                    BPPARAM=BPPARAM)
  list.keepX <- tune.feature$choice.keepX
  final.diablo.model <- block.splsda(X=data_block,
                                     Y=Y,
                                     ncomp=ncomp,
                                     keepX=list.keepX,
                                     design=design)
  model_prefix <- if (site_label == "MTS") "MTS_" else "Throat_"
  model_path <- file.path(res_dir,paste0(model_prefix,time,"_",design_label,"_final_model.rds"))
  saveRDS(final.diablo.model,file=model_path)
  perf_final <- perf(final.diablo.model,validation="loo")
  perf_path <- file.path(res_dir,paste0(model_prefix,time,"_",design_label,"_perf.rds"))
  saveRDS(perf_final,file=perf_path)
  invisible(list(model=final.diablo.model,perf=perf_final,design=design,design_strength=design_strength))
}

# Use final selected DIABLO design strength

design_strength <- 0.25
design_label <- paste0("design_",gsub("\\.","p",as.character(design_strength)))
module_set_name <- "chaussabel"

# Infection susceptibility datasets

micro_mts_inf_all <- readRDS(file.path(diablo_data_dir,"microbiome_mts_species_by_timepoint.rds"))
micro_throat_inf_all <- readRDS(file.path(diablo_data_dir,"microbiome_throat_species_by_timepoint.rds"))

# Symptom/viral load datasets

micro_mts_sx <- readRDS(file.path(diablo_data_dir,"microbiome_mts_species_by_timepoint_infected_only.rds"))
micro_throat_sx <- readRDS(file.path(diablo_data_dir,"microbiome_throat_species_by_timepoint_infected_only.rds"))

micro_mts_vl <- micro_mts_sx
micro_throat_vl <- micro_throat_sx

ssgsea_mts_inf_all <- make_ssgsea_timepoint_list(ssgsea_mts,infected_only=FALSE)
ssgsea_throat_inf_all <- make_ssgsea_timepoint_list(ssgsea_throat,infected_only=FALSE)
ssgsea_mts_sx <- make_ssgsea_timepoint_list(ssgsea_mts,infected_only=TRUE)
ssgsea_throat_sx <- make_ssgsea_timepoint_list(ssgsea_throat,infected_only=TRUE)
ssgsea_mts_vl <- ssgsea_mts_sx
ssgsea_throat_vl <- ssgsea_throat_sx

diablo_results_dir <- file.path(diablo_results_base_dir,paste0("results_",module_set_name))
res_infection_mts_dir <- file.path(diablo_results_dir,"infection","mts")
res_infection_throat_dir <- file.path(diablo_results_dir,"infection","throat")
res_symptoms_mts_dir <- file.path(diablo_results_dir,"symptoms","mts")
res_symptoms_throat_dir <- file.path(diablo_results_dir,"symptoms","throat")
res_viral_mts_dir <- file.path(diablo_results_dir,"viral_load","mts")
res_viral_throat_dir <- file.path(diablo_results_dir,"viral_load","throat")

dir.create(res_infection_mts_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(res_infection_throat_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(res_symptoms_mts_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(res_symptoms_throat_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(res_viral_mts_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(res_viral_throat_dir,showWarnings=FALSE,recursive=TRUE)

message("========================================")
message("Starting DIABLO runs for ",module_set_name," | ",design_label)
message("========================================")

for(time in all_timepoint){
  message("Running infection MTS: ",time," | ",module_set_name," | ",design_label)
  run_diablo_for_timepoint(time,micro_mts_inf_all,ssgsea_mts_inf_all,"infection","MTS",res_infection_mts_dir,design_strength,FALSE)
  message("Running infection Throat: ",time," | ",module_set_name," | ",design_label)
  run_diablo_for_timepoint(time,micro_throat_inf_all,ssgsea_throat_inf_all,"infection","Throat",res_infection_throat_dir,design_strength,TRUE)
}

for(time in all_timepoint){
  message("Running symptoms MTS: ",time," | ",module_set_name," | ",design_label)
  run_diablo_for_timepoint(time,micro_mts_sx,ssgsea_mts_sx,"symptoms","MTS",res_symptoms_mts_dir,design_strength,FALSE)
  message("Running symptoms Throat: ",time," | ",module_set_name," | ",design_label)
  run_diablo_for_timepoint(time,micro_throat_sx,ssgsea_throat_sx,"symptoms","Throat",res_symptoms_throat_dir,design_strength,TRUE)
}

for(time in all_timepoint){
  message("Running viral load MTS: ",time," | ",module_set_name," | ",design_label)
  run_diablo_for_timepoint(time,micro_mts_vl,ssgsea_mts_vl,"viral_load","MTS",res_viral_mts_dir,design_strength,FALSE)
  message("Running viral load Throat: ",time," | ",module_set_name," | ",design_label)
  run_diablo_for_timepoint(time,micro_throat_vl,ssgsea_throat_vl,"viral_load","Throat",res_viral_throat_dir,design_strength,TRUE)
}

message("Completed DIABLO runs for ",module_set_name," | ",design_label)