# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH 
# RNASeq data preprocessing
# Last updated: 10/06/2026

script_summary <- "This script functionally annotates Chaussabel immune transcriptional modules using ReactomePA and Gene Ontology Biological Process
enrichment analyses to identify dominant biological pathways and themes within each module. Final module annotations are informed by pathway enrichment
results and manual review of module genes, and the curated 17-module Chaussabel gene set is saved for downstream analyses. The script then processes
the STAR RNA-seq count table by annotating genes via Ensembl/biomaRt and combining count data, gene annotations, and sample metadata into a phyloseq
object. After removal of negative controls (<300,000 reads), the resulting object is saved as phy.rnaseq.rds. Gene filtering removes non-protein-coding
genes, genes lacking external gene symbols, and mitochondrial and ribosomal genes. All remaining samples have >5 million reads, so no samples are
pruned, and the filtered object is saved as phy.rnaseq.filter.rds. Finally, single-sample gene set enrichment analysis (ssGSEA) is performed separately
for nasal and throat RNA-seq samples using the curated Chaussabel modules. Raw ssGSEA scores are retained, and module enrichment scores are additionally
normalized to a 0-1 range separately within each anatomical site. The resulting longitudinal ssGSEA dataset, together with curated Chaussabel module 
definitions, pathway annotation tables, and processed phyloseq objects, is saved for downstream statistical analyses and figure generation."
cat(script_summary)

remove(list=ls())
version
set.seed(123)
setwd("________________________") 

library(biomaRt)
library(readxl)
library(openxlsx)
library(tidyverse)
library(phyloseq)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(GSVA)
library(reshape2)
library(edgeR)

data_dir <- "________________________"
dir.create(data_dir,showWarnings=FALSE,recursive=TRUE)

loaded_obj <- load("orig_chaussabel.rda")
chaussabel <- get(loaded_obj[1])
names(chaussabel) <- sub("_.*","",names(chaussabel))
chaussabel <- lapply(chaussabel,function(x){
  x <- as.character(x)
  # remove NA / blanks
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[x!=""]
  # remove duplicates
  x <- unique(x)
  return(x)
})

chaussabel_df <- as.data.frame(
  lapply(chaussabel,function(x){
    length(x) <- max(lengths(chaussabel))
    x
  }),
  stringsAsFactors=FALSE
)

mods <- names(chaussabel)
length(mods)
head(mods)

# Run Reactome pathway enrichment analysis 
react_results <- list()
react_tbl <- list()
for(mod in mods){
  cat("\n========================\n")
  cat(mod,"\n")
  cat("========================\n")
  genes <- chaussabel[[mod]]
  ids <- bitr(genes,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)
  if(is.null(ids) || nrow(ids)==0){
    cat("No genes mapped for",mod,"\n")
    react_results[[mod]] <- NULL
    react_tbl[[mod]] <- data.frame(module=mod,status="No genes mapped")
    next
  }
  react <- enrichPathway(gene=unique(ids$ENTREZID),organism="human",pAdjustMethod="BH",qvalueCutoff=0.05,readable=TRUE)
  react_results[[mod]] <- react
  df <- as.data.frame(react)
  if(is.null(df) || nrow(df)==0){
    cat("No enriched Reactome pathways for",mod,"\n")
    react_tbl[[mod]] <- data.frame(module=mod,status="No enriched Reactome pathways")
    next
  }
  top_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n=10) %>%
    mutate(
      module=mod,
      status="Enriched",
      module_gene_count=length(unique(genes)),
      mapped_entrez_count=length(unique(ids$ENTREZID)),
      genes_in_enrichment=as.integer(sub("/.*","",GeneRatio)),
      gene_ratio_denominator=as.integer(sub(".*/","",GeneRatio))
    ) %>%
    dplyr::select(
      module,
      status,
      Description,
      genes_in_enrichment,
      module_gene_count,
      mapped_entrez_count,
      gene_ratio_denominator,
      p.adjust,
      geneID
    )
  print(top_df)
  react_tbl[[mod]] <- top_df
}
react_tbl <- bind_rows(react_tbl)
write.csv(react_tbl,file.path(data_dir,"Chaussabel_Reactome_top10_by_module.csv"),row.names=FALSE)

# Run GO Biological Process enrichment analysis
go_results <- list()
go_tbl <- list()
for(mod in mods){
  cat("\n========================\n")
  cat(mod,"GO BP\n")
  cat("========================\n")
  genes <- chaussabel[[mod]]
  ids <- bitr(genes,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)
  if(is.null(ids) || nrow(ids)==0){
    cat("No genes mapped for",mod,"\n")
    go_results[[mod]] <- NULL
    go_tbl[[mod]] <- data.frame(module=mod,status="No genes mapped")
    next
  }
  go <- enrichGO(gene=unique(ids$ENTREZID),OrgDb=org.Hs.eg.db,ont="BP",pAdjustMethod="BH",qvalueCutoff=0.05,readable=TRUE)
  go_results[[mod]] <- go
  df <- as.data.frame(go)
  if(is.null(df) || nrow(df)==0){
    cat("No enriched GO BP terms for",mod,"\n")
    go_tbl[[mod]] <- data.frame(module=mod,status="No enriched GO BP terms")
    next
  }
  top_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n=10) %>%
    mutate(
      module=mod,
      status="Enriched",
      module_gene_count=length(unique(genes)),
      mapped_entrez_count=length(unique(ids$ENTREZID)),
      genes_in_enrichment=as.integer(sub("/.*","",GeneRatio)),
      gene_ratio_denominator=as.integer(sub(".*/","",GeneRatio))
    ) %>%
    dplyr::select(
      module,
      status,
      Description,
      genes_in_enrichment,
      module_gene_count,
      mapped_entrez_count,
      gene_ratio_denominator,
      p.adjust,
      geneID
    )
  print(top_df)
  go_tbl[[mod]] <- top_df
}
go_tbl <- bind_rows(go_tbl)
write.csv(go_tbl,file.path(data_dir,"Chaussabel_GO_BP_top10_by_module.csv"),row.names=FALSE)

final_module_map <- c(
  "M1.2"  = "M1.2_Type_1_interferon",
  "M3.4"  = "M3.4_ISGs",
  "M5.4"  = "M5.4_Antiviral_sensing",
  "M9.51" = "M9.51_Innate_signaling",
  "M4.13" = "M4.13_Innate_inflammation",
  "M4.6"  = "M4.6_Monocyte/APC_activation",
  "M5.1"  = "M5.1_Phagocytosis",
  "M7.1"  = "M7.1_Neutrophil_migration",
  "M4.9"  = "M4.9_Neutrophil_inflammation",
  "M7.2"  = "M7.2_Neutrophil_activation",
  "M5.7"  = "M5.7_Neutrophil_degranulation",
  "M3.2"  = "M3.2_TLR_signaling",
  "M6.17" = "M6.17_NFkB_signaling",
  "M3.6"  = "M3.6_NK_cells",
  "M4.2"  = "M4.2_B_cells",
  "M4.1"  = "M4.1_T_cells",
  "M4.15" = "M4.15_T_cell_activation"
)

missing_mods <- setdiff(names(final_module_map),names(chaussabel))
if(length(missing_mods)>0){
  stop(paste("Missing modules:",paste(missing_mods,collapse=", ")))
}

chaussabel_final <- chaussabel[names(final_module_map)]
names(chaussabel_final) <- unname(final_module_map)
max_len_final <- max(lengths(chaussabel_final))

chaussabel_final_df <- as.data.frame(
  lapply(chaussabel_final,function(x){
    length(x) <- max_len_final
    x
  }),
  stringsAsFactors=FALSE
)

saveRDS(chaussabel_final,file.path(data_dir,"chaussabel_modules_17.rds"))
openxlsx::write.xlsx(chaussabel_final_df,file.path(data_dir,"chaussabel_modules_17.xlsx"),overwrite=TRUE)

###########################################################
# GENERATE INITIAL PHYLOSEQ OBJECT CONTAINING ALL SAMPLES #
###########################################################

# Create RNASeq count table
rnaseq <- read.csv("SARS2_Challenge_STAR_output.csv",sep=",",stringsAsFactors=FALSE,check.names=FALSE)
rnaseq$target_id <- gsub("\\..*","",rnaseq$target_id)
count_cols <- setdiff(colnames(rnaseq),"target_id")
rnaseq_counts <- rnaseq %>%
  group_by(target_id) %>%
  summarise(across(all_of(count_cols),\(x) sum(x,na.rm=TRUE)),.groups="drop") %>%
  column_to_rownames("target_id")
summary(colSums(rnaseq_counts))
ensembl_ids <- rownames(rnaseq_counts)
mart <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl",version = 114)
gene_ids <- getBM(attributes=c("external_gene_name","ensembl_gene_id","ensembl_gene_id_version","gene_biotype","chromosome_name","start_position",
                               "end_position"),filters="ensembl_gene_id",values=ensembl_ids,mart=mart)
gene_ids$Length <- gene_ids$end_position - gene_ids$start_position + 1
summary(gene_ids$Length)
rnaseq_counts <- rnaseq_counts %>% rownames_to_column("ensembl_gene_id")
rnaseq_final <- merge(gene_ids,rnaseq_counts,by="ensembl_gene_id")
rnaseq_final$ensembl_gene_id <- NULL
rnaseq_final$start_position <- NULL
rnaseq_final$end_position <- NULL
names(rnaseq_final)[names(rnaseq_final)=="ensembl_gene_id_version"] <- "EnsemblID"
names(rnaseq_final)[names(rnaseq_final)=="external_gene_name"] <- "Gene"
names(rnaseq_final)[names(rnaseq_final)=="chromosome_name"] <- "Chromosome"
names(rnaseq_final)[names(rnaseq_final)=="gene_biotype"] <- "Biotype"
colnames(rnaseq_final) <- gsub("-","_",colnames(rnaseq_final))

# Create gene count table (otu_table)
counts_raw <- rnaseq_final[,-c(1,3,4,5)] # dropping Gene, Biotype, Chromosome, Length
counts_raw <- counts_raw %>% remove_rownames %>% column_to_rownames(var="EnsemblID")
counts_raw_otu <- otu_table(counts_raw, taxa_are_rows = TRUE)

# Create gene annotation table (tax_table)
genes <- rnaseq_final[,c(1:5)] # Keeping Gene, EnsemblID, Biotype, Chromosome, Length
genes <- genes %>% remove_rownames %>% column_to_rownames(var="EnsemblID")
genes_tax <- tax_table(as.matrix(genes))

# Upload metadata file
samples <- read_excel("________________________/SARS2_challenge_metadata.xlsx", sheet="samples", na="NA")
controls <- read_excel("________________________/SARS2_challenge_RNASeq_controls.xlsx", sheet="controls", na="NA")
metadata <- rbind(samples, controls)
metadata <- metadata %>% column_to_rownames(var = "sampleID")

# Assemble initial phyloseq object
phy.rnaseq <- phyloseq(counts_raw_otu, genes_tax)
sample_data(phy.rnaseq) <- metadata
nsamples(phy.rnaseq)
ntaxa(phy.rnaseq)
# Negative extraction and library prep controls have <300,000 reads
# Remove negative control samples
phy.rnaseq <- subset_samples(phy.rnaseq, neg_control==FALSE)
nsamples(phy.rnaseq)
saveRDS(phy.rnaseq,file.path(data_dir,"phy.rnaseq.rds"))
remove(mart, rnaseq, ensembl_ids, gene_ids, counts_raw, genes, rnaseq_final)

##################
# GENE FILTERING #
##################

phy.rnaseq <- readRDS(file.path(data_dir,"phy.rnaseq.rds"))
genes <- data.frame(tax_table(phy.rnaseq))
summary(sample_sums(phy.rnaseq))
nsamples(phy.rnaseq)
ntaxa(phy.rnaseq)

# Remove transcripts for non-protein-coding genes
phy.rnaseq.filter <- subset_taxa(phy.rnaseq, Biotype=="protein_coding")
ntaxa(phy.rnaseq.filter)

# Remove transcripts without gene identifier
phy.rnaseq.filter <- subset_taxa(phy.rnaseq.filter, Gene!="")
ntaxa(phy.rnaseq.filter)

# Remove mitochondrial genes
MT_genes <- genes[grep("^MT", genes$Gene, perl=T),]
MT_genes <- MT_genes$Gene
phy.rnaseq.filter <- subset_taxa(phy.rnaseq.filter, !(Gene %in% MT_genes))
ntaxa(phy.rnaseq.filter)
# Remove large subunit mitochondrial ribosomal proteins
MRPL_genes <- genes[grep("^MRPL", genes$Gene, perl=T),]
MRPL_genes <- MRPL_genes$Gene
phy.rnaseq.filter <- subset_taxa(phy.rnaseq.filter, !(Gene %in% MRPL_genes))
ntaxa(phy.rnaseq.filter)
# Remove small subunit mitochondrial ribosomal proteins
MRPS_genes <- genes[grep("^MRPS", genes$Gene, perl=T),]
MRPS_genes <- MRPS_genes$Gene
phy.rnaseq.filter <- subset_taxa(phy.rnaseq.filter, !(Gene %in% MRPS_genes))
ntaxa(phy.rnaseq.filter)
# Remove L ribosomal genes
RPL_genes <- genes[grep("^RPL", genes$Gene, perl=T),]
RPL_genes <- RPL_genes$Gene
phy.rnaseq.filter <- subset_taxa(phy.rnaseq.filter, !(Gene %in% RPL_genes))
ntaxa(phy.rnaseq.filter)
# Remove S ribosomal genes
RPS_genes <- genes[grep("^RPS", genes$Gene, perl=T),]
RPS_genes <- RPS_genes$Gene
phy.rnaseq.filter <- subset_taxa(phy.rnaseq.filter, !(Gene %in% RPS_genes))
ntaxa(phy.rnaseq.filter)

# Remove samples with <5 million reads
nsamples(phy.rnaseq.filter)
summary(sample_sums(phy.rnaseq.filter))
phy.rnaseq.filter <- subset_samples(phy.rnaseq.filter, sample_sums(phy.rnaseq.filter)>=5000000)
nsamples(phy.rnaseq.filter)
# No samples removed as all have >=5 million reads
nsamples(subset_samples(phy.rnaseq.filter, sampletype=="MTS"))
nsamples(subset_samples(phy.rnaseq.filter, sampletype=="Throat"))
ntaxa(phy.rnaseq.filter)
summary(sample_sums(phy.rnaseq.filter))
summary(sample_sums(subset_samples(phy.rnaseq.filter, sampletype=="MTS")))
summary(sample_sums(subset_samples(phy.rnaseq.filter, sampletype=="Throat")))
saveRDS(phy.rnaseq.filter,file.path(data_dir,"phy.rnaseq.filter.rds"))
remove(genes)

# Generate immune module expression scores using ssGSEA

phy.rnaseq.filter <- readRDS(file.path(data_dir,"phy.rnaseq.filter.rds"))
gene_sets <- readRDS(file.path(data_dir,"chaussabel_modules_17.rds"))

all_timepoint <- c("M1",paste0(rep(c(1,2,3,5,7,14),each=2),c("AM","PM")))[-c(9,11,13)]

build_ssgsea_long_from_phyloseq <- function(phy,gene_sets,site_label){
  meta_all <- data.frame(sample_data(phy))
  keep_samples <- rownames(meta_all)[meta_all$sampletype == site_label]
  phy_site <- prune_samples(keep_samples,phy)
  geneinfo <- data.frame(tax_table(phy_site))
  countdata <- data.frame(otu_table(phy_site))
  metadata <- data.frame(sample_data(phy_site))
  geneinfo$Gene <- make.unique(geneinfo$Gene)
  rownames(countdata) <- geneinfo$Gene
  metadata$timepoint <- factor(metadata$timepoint,levels=all_timepoint)
  metadata$infection_status <- ifelse(metadata$infected == "N","uninfected","infected")
  dge <- edgeR::DGEList(counts=countdata)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  expr_mat <- edgeR::cpm(dge,log=TRUE,prior.count=1)
  param <- GSVA::ssgseaParam(exprData=expr_mat,geneSets=gene_sets,normalize=TRUE)
  ssgsea <- t(as.matrix(GSVA::gsva(param,verbose=FALSE)))
  df_ssgsea <- reshape2::melt(ssgsea,varnames=c("sample","module"),value.name="ssgsea")
  meta_with_id <- metadata %>% tibble::rownames_to_column("sample")
  df_ssgsea %>%
    dplyr::left_join(meta_with_id,by="sample") %>%
    dplyr::mutate(site=site_label,timepoint=factor(timepoint,levels=all_timepoint),
                  infected=factor(infection_status,levels=c("uninfected","infected")),sex=factor(sex))
}

add_site_scaled_ssgsea <- function(ssgsea_long_df){
  ssgsea_min <- min(ssgsea_long_df$ssgsea,na.rm=TRUE)
  ssgsea_max <- max(ssgsea_long_df$ssgsea,na.rm=TRUE)
  ssgsea_long_df %>%
    dplyr::mutate(
      ssgsea_scaled=(ssgsea-ssgsea_min)/(ssgsea_max-ssgsea_min),
      ssgsea_scaled=pmin(pmax(ssgsea_scaled,0),1)
    )
}

summarize_ssgsea_distribution <- function(ssgsea_df,label){
  cat("\n==============================\n")
  cat(label,"\n")
  cat("==============================\n")
  cat("\nRaw ssGSEA summary:\n")
  print(summary(ssgsea_df$ssgsea))
  cat("Raw ssGSEA range:\n")
  print(range(ssgsea_df$ssgsea,na.rm=TRUE))
  cat("\nScaled ssGSEA summary:\n")
  print(summary(ssgsea_df$ssgsea_scaled))
  cat("Scaled ssGSEA range:\n")
  print(range(ssgsea_df$ssgsea_scaled,na.rm=TRUE))
  cat("\nScaled ssGSEA by module:\n")
  print(
    ssgsea_df %>%
      dplyr::group_by(module) %>%
      dplyr::summarise(
        min=min(ssgsea_scaled,na.rm=TRUE),
        q25=quantile(ssgsea_scaled,0.25,na.rm=TRUE),
        median=median(ssgsea_scaled,na.rm=TRUE),
        mean=mean(ssgsea_scaled,na.rm=TRUE),
        q75=quantile(ssgsea_scaled,0.75,na.rm=TRUE),
        max=max(ssgsea_scaled,na.rm=TRUE),
        .groups="drop"
      ) %>%
      dplyr::arrange(module),
    n=Inf
  )
}

ssgsea_mts <- build_ssgsea_long_from_phyloseq(
  phy=phy.rnaseq.filter,
  gene_sets=gene_sets,
  site_label="MTS"
) %>%
  add_site_scaled_ssgsea()

ssgsea_throat <- build_ssgsea_long_from_phyloseq(
  phy=phy.rnaseq.filter,
  gene_sets=gene_sets,
  site_label="Throat"
) %>%
  add_site_scaled_ssgsea()

summarize_ssgsea_distribution(ssgsea_mts,"MTS ssGSEA")
summarize_ssgsea_distribution(ssgsea_throat,"Throat ssGSEA")

p_mts_ssgsea <- ggplot(ssgsea_mts,aes(x=ssgsea_scaled)) +
  geom_histogram(bins=40,color="black",fill="grey80") +
  theme_bw() +
  labs(title="MTS scaled ssGSEA distribution",x="Scaled ssGSEA",y="Count")

p_throat_ssgsea <- ggplot(ssgsea_throat,aes(x=ssgsea_scaled)) +
  geom_histogram(bins=40,color="black",fill="grey80") +
  theme_bw() +
  labs(title="Throat scaled ssGSEA distribution",x="Scaled ssGSEA",y="Count")

print(p_mts_ssgsea)
print(p_throat_ssgsea)

ssgsea <- dplyr::bind_rows(ssgsea_mts,ssgsea_throat)
saveRDS(ssgsea,file.path(data_dir,"ssgsea.rds"))