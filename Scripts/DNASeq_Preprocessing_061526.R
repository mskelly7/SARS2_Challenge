# SARS-CoV-2 Challenge Study
# Matthew Kelly, MD, MPH 
# Microbiome sequencing data preprocessing
# Last updated: 06/15/2026

script_summary <- "This script preprocesses microbiome sequencing data from throat and nasal swab samples and creates finalized phyloseq objects for 
downstream analyses. It defines helper functions to collapse taxa to the species level and to create qPCR-scaled abundance objects by multiplying 
relative abundances by 16S qPCR bacterial load. For throat shotgun metagenomic data, the script reads KneadData read-count summaries, Bracken taxonomic 
profiles, and sample metadata. The script then examines negative controls, removes host reads, and saves the initial object as phy.throat.ngs.rds. 
It then uses decontam2 to identify reagent contaminants using negative-control prevalence (threshold 0.1) and DNA-concentration frequency (threshold 0.1) 
approaches, removes contaminant taxa, and saves phy.throat.ngs.decontam.rds. Control samples are removed and diversity-read-depth plots are used to 
inform a 10,000-read cutoff for sample pruning, although no samples are removed based on this cutoff. Species were retained if they demonstrated mean 
relative abundance >0.5% and were detected in >10% of samples; remaining taxa were collapsed into an 'Other' category. The final objects are saved as 
phy.throat.ngs.pruned.species.rds and the qPCR-scaled phy.throat.ngs.pruned.species.qpcr.rds. For nasal MTS V1-V3 16S rRNA sequencing data, the script 
reads DADA2/eHOMD ASV taxonomy and sequence tables and sample metadata and saves phy.mts.16s.rds. It then uses decontam2 to identify contaminants by 
negative-control prevalence (0.3) and DNA-concentration frequency (0.3) approaches, manually retains several Staphylococcus/Lawsonella ASVs after 
plot review, and saves phy.mts.16s.decontam.rds. Control samples are removed and read-depth diversity plots are used to inform a 25,000-read cutoff 
for nasal sample pruning. Species were retained if they demonstrated mean relative abundance >0.5% and were detected in >10% of samples; remaining taxa 
were collapsed into an 'Other' category. Final objects are saved as phy.mts.16s.pruned.species.rds and the qPCR-scaled phy.mts.16s.pruned.species.qpcr.rds."
cat(script_summary)

remove(list=ls())
version
set.seed(123)
setwd("________________________") 

library(phyloseq)
library(tidyverse)
library(readxl)
library(decontam)
library(data.table)
library(Biostrings)

data_dir <- "________________________"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

make_qpcr_phyloseq <- function(ps, qpcr_col="qpcr_log10"){
  meta <- data.frame(sample_data(ps))
  if(!(qpcr_col %in% colnames(meta))) stop(paste(qpcr_col,"not found in sample_data"))
  if(any(is.na(meta[[qpcr_col]]))) stop("Missing qpcr_log10 values in pruned phyloseq object")
  ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
  otu_mat <- as(otu_table(ps_rel), "matrix")
  qpcr_counts <- 10^meta[sample_names(ps_rel), qpcr_col]
  if(taxa_are_rows(ps_rel)){
    otu_mat <- sweep(otu_mat, 2, qpcr_counts, "*")
  } else {
    otu_mat <- sweep(otu_mat, 1, qpcr_counts, "*")
  }
  otu_table(ps_rel) <- otu_table(otu_mat, taxa_are_rows=taxa_are_rows(ps_rel))
  sample_data(ps_rel) <- sample_data(ps)
  ps_rel
}

collapse_phyloseq_to_species <- function(ps){
  tax_df <- as.data.frame(tax_table(ps))
  tax_df$Species_full <- paste(tax_df$Genus,tax_df$Species)
  tax_df$Species_full[tax_df$Species_full %in% c("NA NA","Other Other")] <- "Other"
  tax_df$Species_full <- trimws(tax_df$Species_full)
  otu_mat <- as(otu_table(ps),"matrix")
  if(taxa_are_rows(ps)){
    otu_species <- rowsum(otu_mat,group=tax_df$Species_full)
    taxa_are_rows_out <- TRUE
  } else {
    otu_species <- t(rowsum(t(otu_mat),group=tax_df$Species_full))
    taxa_are_rows_out <- FALSE
  }
  tax_species <- tax_df %>%
    dplyr::mutate(Species=Species_full) %>%
    dplyr::group_by(Species) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Species_full) %>%
    as.data.frame()
  rownames(tax_species) <- tax_species$Species
  phyloseq(
    otu_table(otu_species,taxa_are_rows=taxa_are_rows_out),
    tax_table(as.matrix(tax_species)),
    sample_data(ps)
  )
}

######################################################
# THROAT SWAB SAMPLES - SHOTGUN METAGENOMIC SEQUENCING
######################################################

#########################
# CREATE PHYLOSEQ OBJECT
#########################

counts <- read.csv("Pipeline_Outputs/SARS2_Challenge_NGS_kneaddata_counts.csv", sep=",", stringsAsFactors = FALSE, na.strings=("NA"))
counts <- counts[,c("Sample","raw.pair1","trimmed.pair1","final.pair1")]
counts$pct_metagenomic <- ((counts$final.pair1)/(counts$trimmed.pair1))*100
names(counts)[names(counts)=="Sample"] <- "sampleID"
names(counts)[names(counts)=="raw.pair1"] <- "raw_reads"
names(counts)[names(counts)=="trimmed.pair1"] <- "trimmed_reads"
names(counts)[names(counts)=="final.pair1"] <- "microbial_reads_kneaddata"
counts$sampleID <- gsub("_NGS_R1_kneaddata","",as.factor(counts$sampleID))
counts$sampleID <- gsub("-","_",as.factor(counts$sampleID))
nrow(counts)
summary(counts$raw_reads)
summary(counts$microbial_reads_kneaddata)
summary(counts$pct_metagenomic)

# Upload metadata and bracken output files
samples <- read_excel("________________________/SARS2_challenge_metadata.xlsx", sheet="samples", na="NA")
controls <- read_excel("________________________/SARS2_challenge_DNASeq_controls.xlsx", sheet="controls", na="NA")
metadata <- rbind(samples, controls)
bracken <- read.csv("Pipeline_Outputs/SARS2_Challenge_NGS_bracken.csv", sep=",", stringsAsFactors = FALSE, na.strings=("NA"))
colnames(bracken) <- gsub("_NGS","", colnames(bracken))

# Create initial OTU table file
otutable <- bracken[,c(8,10:911)] # tax_id followed by sample IDs
row.names(otutable) <- otutable[,1]
otutable$tax_id <- NULL
otutable <- otutable %>% dplyr:: select(contains("_num"))
colnames(otutable) <- gsub("_num","", colnames(otutable))
sampleID <- colnames(otutable)
microbial_reads_bracken <- colSums(otutable)
bracken_counts <- data.frame(sampleID, microbial_reads_bracken)
metadata_counts <- merge(metadata, counts, by="sampleID", all.x=F)
metadata_counts <- merge(metadata_counts, bracken_counts, by="sampleID", all.x=F)
metadata_counts <- metadata_counts %>% column_to_rownames(var = "sampleID")

# Create initial taxonomy table
taxtable <- bracken[,c(8,1:7)] # tax_id followed by taxonomic classification
row.names(taxtable) <- taxtable[,1]
taxtable$tax_id <- NULL
names(taxtable)[names(taxtable)=="domain"] <- "Kingdom"
names(taxtable)[names(taxtable)=="phylum"] <- "Phylum"
names(taxtable)[names(taxtable)=="class"] <- "Class"
names(taxtable)[names(taxtable)=="order"] <- "Order"
names(taxtable)[names(taxtable)=="family"] <- "Family"
names(taxtable)[names(taxtable)=="species"] <- "Species"
# Replace NA with data from higher taxonomic rank (identifies taxonomic rank with letter)
rplc <- which(taxtable$Phylum=="" | is.na(taxtable$Phylum))
taxtable$Phylum[rplc] <- paste(taxtable$Kingdom[rplc],";p",sep="")
rplc <- which(taxtable$Class=="" | is.na(taxtable$Class))
taxtable$Class[rplc] <- paste(taxtable$Phylum[rplc],";c",sep="")
rplc <- which(taxtable$Order=="" | is.na(taxtable$Order))
taxtable$Order[rplc] <- paste(taxtable$Class[rplc],";o",sep="")
rplc <- which(taxtable$Family=="" | is.na(taxtable$Family))
taxtable$Family[rplc] <- paste(taxtable$Order[rplc],";f",sep="")
rplc <- which(taxtable$Species=="" | is.na(taxtable$Species) | taxtable$Species=="")
taxtable$Species[rplc] <- paste(taxtable$Family[rplc],";g",sep="")
rplc <- which(taxtable$Species=="" | is.na(taxtable$Species) | taxtable$Species=="")
taxtable$Species[rplc] <- paste(taxtable$Species[rplc],";s",sep="")
remove(rplc)

# Assemble phyloseq objects with taxa read counts
taxtable <- tax_table(as.matrix(taxtable))
otutable <- otu_table(otutable, taxa_are_rows = TRUE)
phy.throat.ngs <- phyloseq(otutable, taxtable)
sample_data(phy.throat.ngs) <- metadata_counts
nsamples(phy.throat.ngs)
ntaxa(phy.throat.ngs)
summary(sample_sums(phy.throat.ngs))

# Examine negative control samples
neg_controls <- str_subset(sampleID, "neg")
neg_controls <- append(neg_controls,"Species",after=0)
bracken_neg_controls <- bracken %>% dplyr::select(7, matches(neg_controls))
bracken_neg_controls <- dplyr::select(bracken_neg_controls, -contains("frac"))
bracken_neg_controls <- aggregate(.~species, bracken_neg_controls, sum)
bracken_neg_controls <- bracken_neg_controls %>% column_to_rownames(var = "species")
summary(colSums(bracken_neg_controls))
# Calculate mean # of counts for taxa in negative control samples
mean_abundance <- rowMeans(bracken_neg_controls)
sorted_abundance <- sort(mean_abundance, decreasing = TRUE)
top10 <- head(sorted_abundance, 10)
top10_df <- data.frame(Taxon = names(top10),
                       mean_counts = as.numeric(top10),
                       row.names = NULL)
print(top10_df)
# Most reads from negative control samples are assigned as common reagent contaminants

# Remove human reads and microbial reads from kingdoms other than Archaea, Bacteria, or Eukaryota
ntaxa(phy.throat.ngs)
phy.throat.ngs <- subset_taxa(phy.throat.ngs, Kingdom=="Archaea" | Kingdom=="Bacteria" | Kingdom=="Eukaryota")
ntaxa(phy.throat.ngs)
phy.throat.ngs <- subset_taxa(phy.throat.ngs, Species!="Homo sapiens")
ntaxa(phy.throat.ngs)
summary(sample_sums(phy.throat.ngs))
saveRDS(phy.throat.ngs,file.path(data_dir, "phy.throat.ngs.rds"))

remove(bracken_counts, counts, metadata, samples, controls, otutable, microbial_reads_bracken, taxtable, neg_controls)

##########################################################################
# REMOVAL OF REAGENT CONTAMINANTS USING DECONTAM (USE RELATIVE ABUNDANCES)
##########################################################################

set.seed(1234)
library(decontam); packageVersion("decontam")

phy.throat.ngs <- readRDS(file.path(data_dir, "phy.throat.ngs.rds"))
metadata_counts <- data.frame(sample_data(phy.throat.ngs))
metadata_counts$sampleID <- row.names(metadata_counts)

# First create a figure that shows the distribution of sample read counts
ggplot(data=metadata_counts, aes(x=sampleID, y=microbial_reads_bracken, color=control)) + geom_point()
ggplot(data=metadata_counts, aes(x=sampleID, y=microbial_reads_bracken, color=neg_control)) + geom_point()

# Identify contaminants based on prevalence in negative controls (use relative abundances)
phy.throat.ngs.relative <- transform_sample_counts(phy.throat.ngs, function(Abundance) Abundance/sum(Abundance))
taxtable <- data.frame(tax_table(phy.throat.ngs.relative))
suppressWarnings(prev <- isContaminant(phy.throat.ngs.relative, method="prevalence", neg="neg_control", batch="batch_ngs", threshold=0.1))
table(prev$contaminant)
contam_prev <- subset(prev, contaminant==TRUE)
contam_prev <- merge(taxtable, contam_prev, by="row.names")
ntaxa(phy.throat.ngs)
phy.throat.ngs.decontam.prev <- prune_taxa(!prev$contaminant, phy.throat.ngs)
ntaxa(phy.throat.ngs.decontam.prev)

# Identify contaminants based on frequency (use relative abundances)
suppressWarnings(freq <- isContaminant(phy.throat.ngs.decontam.prev, method="frequency", conc="dna_ngs", batch="batch_ngs", threshold=0.1))
table(freq$contaminant)
contam_freq <- subset(freq, contaminant==TRUE)
contam_freq <- merge(taxtable, contam_freq, by="row.names")
names(contam_freq)[names(contam_freq)=="Row.names"] <- "Taxon"
contam_freq <- contam_freq[,-c(2:7)]
freq_abundant <- freq
freq_abundant$contaminant[freq_abundant$freq<0.001] <- FALSE
contaminant_plot <- plot_frequency(phy.throat.ngs.decontam.prev, taxa_names(phy.throat.ngs.decontam.prev)[sample(which(freq_abundant$contaminant))], conc="dna_ngs") + xlab("DNA Concentration")
png(file="R_Plots/contaminant_plot_throat_ngs.png", width = 10, height = 10, units = 'in', res = 600)
print(contaminant_plot)
dev.off()
ntaxa(phy.throat.ngs.decontam.prev)
phy.throat.ngs.decontam <- prune_taxa(!freq$contaminant, phy.throat.ngs.decontam.prev)
ntaxa(phy.throat.ngs.decontam)

# Create phyloseq object after reagent contaminant removal
nsamples(phy.throat.ngs.decontam)
sample_sums(phy.throat.ngs.decontam)
decontam_reads <- sample_sums(phy.throat.ngs.decontam)
decontam_counts <- data.frame(sampleID, decontam_reads)
metadata_counts <- merge(metadata_counts, decontam_counts, by="sampleID", all.x=F)
rownames(metadata_counts) <- metadata_counts[,1]
sample_data(phy.throat.ngs.decontam) <- metadata_counts
saveRDS(phy.throat.ngs.decontam, file.path(data_dir, "phy.throat.ngs.decontam.rds"))

remove(decontam_counts, freq, decontam_reads, phy.throat.ngs.relative, phy.throat.ngs.decontam.prev, prev, taxtable, contaminant_plot)

#################
# SAMPLE PRUNING
#################

phy.throat.ngs.decontam <- readRDS(file.path(data_dir, "phy.throat.ngs.decontam.rds"))
# Remove control samples 
nsamples(phy.throat.ngs.decontam)
phy.throat.ngs.decontam <- subset_samples(phy.throat.ngs.decontam, control==FALSE)
nsamples(phy.throat.ngs.decontam)
metadata <- data.frame(sample_data(phy.throat.ngs.decontam))
decontam_reads <- sample_sums(phy.throat.ngs.decontam)

# Create dataframe with Shannon index and observed species for each specimen
suppressWarnings(diversity <- estimate_richness(phy.throat.ngs.decontam, measures = c("Shannon", "Observed")))
setDT(diversity, keep.rownames = TRUE)[]
diversity$sampleID <- diversity$rn
diversity$rn <- NULL
metadata <- merge(metadata, diversity, by="sampleID")
row.names(metadata) <- metadata$sampleID
sample_data(phy.throat.ngs.decontam) <- metadata
remove(diversity)

# Plots of alpha diversity measures by sample read counts
brks <- c(10000,25000,50000,100000,500000,1000000,5000000,25000000,100000000)
main_brks <- c(10000,100000,1000000,100000000)
shannon <- ggplot(metadata, aes(decontam_reads, Shannon)) +
  geom_point(alpha=0.7, size=1.8) +
  geom_vline(xintercept=main_brks, linetype="dashed", linewidth=0.3) +
  scale_x_log10(breaks=brks, limits=c(10000,100000000), expand=c(0.01,0)) +
  theme_bw()
png(file="R_Plots/shannon_throat_ngs.png", width=12, height=6, units='in', res=600)
print(shannon)
dev.off()

observed <- ggplot(metadata, aes(decontam_reads, Observed)) +
  geom_point(alpha=0.7, size=1.8) +
  geom_vline(xintercept=main_brks, linetype="dashed", linewidth=0.3) +
  scale_x_log10(breaks=brks, limits=c(10000,100000000), expand=c(0.01,0)) +
  coord_cartesian(ylim=c(0,10000)) +
  theme_bw()
png(file="R_Plots/observed_throat_ngs.png", width=12, height=6, units='in', res=600)
print(observed)
dev.off()
remove(shannon, observed)

# Based on this plot, and prior studies on shallow metagenomics, use a cutoff of 10,000 reads for sample pruning
# https://www.nature.com/articles/s41598-023-33489-1
nsamples(phy.throat.ngs.decontam)
summary(sample_sums(phy.throat.ngs.decontam))
phy.throat.ngs.pruned <- prune_samples(sample_sums(phy.throat.ngs.decontam)>=10000, phy.throat.ngs.decontam)
nsamples(phy.throat.ngs.pruned)
# No samples were pruned as the minimum microbial read count is 98,427 reads

# Remove taxa with no reads in any remaining samples
ntaxa(phy.throat.ngs.pruned)
phy.throat.ngs.pruned <- filter_taxa(phy.throat.ngs.pruned, function(x) sum(x) > 0, TRUE)
nsamples(phy.throat.ngs.pruned)
ntaxa(phy.throat.ngs.pruned)
metadata_pruned <- data.frame(sample_data(phy.throat.ngs.pruned))
sum(metadata_pruned$raw_reads)
summary(metadata_pruned$raw_reads)
summary(metadata_pruned$pct_metagenomic)
summary(metadata_pruned$decontam_reads)
summary(metadata_pruned$Shannon)
summary(metadata_pruned$Observed)
saveRDS(phy.throat.ngs.pruned,file.path(data_dir, "phy.throat.ngs.pruned.rds"))

# Final dataset aggregated at species level with read counts
phy.throat.ngs.pruned <- readRDS(file.path(data_dir, "phy.throat.ngs.pruned.rds"))
ntaxa(phy.throat.ngs.pruned)
phy.throat.ngs.pruned.species <- collapse_phyloseq_to_species(phy.throat.ngs.pruned)
tax.df <- as.data.frame(tax_table(phy.throat.ngs.pruned.species))
tax.df[] <- lapply(tax.df,as.character)
no.species <- is.na(tax.df$Species) | tax.df$Species == "" | tax.df$Species == "NA" | tax.df$Species == "Other" | grepl("uncultured|unclassified|unidentified|metagenome|Candidatus",tax.df$Species,ignore.case=TRUE) | grepl(";s$",tax.df$Species) | grepl("^([A-Za-z]+) \\1;s$",tax.df$Species) | (grepl(" sp\\.",tax.df$Species) & !grepl("^Schaalia sp\\. HMT-172$",tax.df$Species,ignore.case=TRUE))
tax.df$Species[no.species] <- "Other"
tax.df[no.species,colnames(tax.df)!="Species"] <- NA
tax_table(phy.throat.ngs.pruned.species) <- tax_table(as.matrix(tax.df))
phy.rel <- transform_sample_counts(phy.throat.ngs.pruned.species,function(x) x/sum(x))
otu.rel <- as(otu_table(phy.rel),"matrix")
if(!taxa_are_rows(phy.rel)) otu.rel <- t(otu.rel)
mean.rel <- rowMeans(otu.rel)
prev <- rowMeans(otu.rel>0)
keep.taxa <- names(mean.rel)[mean.rel>0.005 & prev>0.10]
tax.df <- as.data.frame(tax_table(phy.throat.ngs.pruned.species))
tax.df[] <- lapply(tax.df,as.character)
tax.df$Species[!(rownames(tax.df) %in% keep.taxa)] <- "Other"
tax.df[!(rownames(tax.df) %in% keep.taxa),colnames(tax.df)!="Species"] <- NA
tax_table(phy.throat.ngs.pruned.species) <- tax_table(as.matrix(tax.df))
phy.throat.ngs.pruned.species <- tax_glom(phy.throat.ngs.pruned.species,taxrank="Species",NArm=FALSE)
taxa_names(phy.throat.ngs.pruned.species) <- make.unique(as.character(tax_table(phy.throat.ngs.pruned.species)[,"Species"]))
nsamples(phy.throat.ngs.pruned.species)
ntaxa(phy.throat.ngs.pruned.species)
summary(sample_sums(phy.throat.ngs.pruned.species))
phy.rel.final <- transform_sample_counts(phy.throat.ngs.pruned.species,function(x) x/sum(x))
otu.rel.final <- as(otu_table(phy.rel.final),"matrix")
if(!taxa_are_rows(phy.rel.final)) otu.rel.final <- t(otu.rel.final)
mean.rel.final <- rowMeans(otu.rel.final)
sort(mean.rel.final,decreasing=TRUE)
names(sort(mean.rel.final,decreasing=TRUE))
saveRDS(phy.throat.ngs.pruned.species,file.path(data_dir,"phy.throat.ngs.pruned.species.rds"))

# Create a dataset with "absolute" species abundances obtained through relative abundances x qPCR
nsamples(phy.throat.ngs.pruned.species)
phy.throat.ngs.pruned.species.qpcr <- prune_samples(!is.na(sample_data(phy.throat.ngs.pruned.species)$qpcr_log10), phy.throat.ngs.pruned.species)
phy.throat.ngs.pruned.species.qpcr <- make_qpcr_phyloseq(phy.throat.ngs.pruned.species.qpcr)
nsamples(phy.throat.ngs.pruned.species.qpcr)
summary(sample_sums(phy.throat.ngs.pruned.species.qpcr))
saveRDS(phy.throat.ngs.pruned.species.qpcr, file.path(data_dir, "phy.throat.ngs.pruned.species.qpcr.rds"))

#####################################################
# NASAL SWAB SAMPLES - V1-V3 16S RRNA GENE SEQUENCING
#####################################################

#########################
# CREATE PHYLOSEQ OBJECT
#########################

# Upload sequencing file with assigned taxonomy from eHOMD
mts.taxa.ehomd <- readRDS("16S_Sequencing_Files/mts.taxa.ehomd.rds")
mts.seqtab.nochim <- readRDS("16S_Sequencing_Files/mts.seqtab.nochim.rds")
rownames(mts.seqtab.nochim) <- gsub("_F_filt.fastq.gz","",rownames(mts.seqtab.nochim))
nrow(mts.seqtab.nochim)
row.names(mts.seqtab.nochim)

# Create sample metadata dataframe
samples <- read_excel("________________________/SARS2_challenge_metadata.xlsx", sheet="samples", na="NA")
controls <- read_excel("________________________/SARS2_challenge_DNASeq_controls.xlsx", sheet="controls", na="NA")
metadata_16s <- rbind(samples, controls)
library(tibble)
metadata_16s <- data.frame(column_to_rownames(metadata_16s, var = "sampleID"))

taxtable <- as.data.frame(mts.taxa.ehomd)
taxtable$Species <- NULL
names(taxtable)[7] <- "Species"
# Replace NA in taxtable
rplc <- which(taxtable$Phylum=="" | is.na(taxtable$Phylum))
taxtable$Phylum[rplc] <- paste(taxtable$Kingdom[rplc],";p",sep="")
rplc <- which(taxtable$Class=="" | is.na(taxtable$Class))
taxtable$Class[rplc] <- paste(taxtable$Phylum[rplc],";c",sep="")
rplc <- which(taxtable$Order=="" | is.na(taxtable$Order))
taxtable$Order[rplc] <- paste(taxtable$Class[rplc],";o",sep="")
rplc <- which(taxtable$Family=="" | is.na(taxtable$Family))
taxtable$Family[rplc] <- paste(taxtable$Order[rplc],";f",sep="")
rplc <- which(taxtable$Genus=="" | is.na(taxtable$Genus) | taxtable$Genus=="")
taxtable$Genus[rplc] <- paste(taxtable$Family[rplc],";g",sep="")
rplc <- which(taxtable$Genus=="NA;c;o;f;g")
taxtable$Genus[rplc] <- "Unassigned"
rplc <- which(taxtable$Species=="" | is.na(taxtable$Species) | taxtable$Species=="")
taxtable$Species[rplc] <- paste(taxtable$Genus[rplc],";s",sep="")
remove(rplc)
taxtable <- as.matrix(taxtable)

# Generate phyloseq object
phy.mts.16s <- phyloseq(otu_table(mts.seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxtable), sample_data(metadata_16s))
nsamples(phy.mts.16s)
ntaxa(phy.mts.16s)

# Remove samples with 0 reads (7 negative extraction controls with 0 reads)
phy.mts.16s <- prune_samples(sample_sums(phy.mts.16s)>0, phy.mts.16s)
nsamples(phy.mts.16s)
# Remove 2 throat samples that underwent 16S rRNA sequencing
phy.mts.16s <- subset_samples(phy.mts.16s, sampletype!="Throat")
nsamples(phy.mts.16s)
# Remove taxa with no reads after sample pruning
phy.mts.16s <- filter_taxa(phy.mts.16s, function(x) sum(x) > 0, TRUE)
ntaxa(phy.mts.16s)

# Number ASVs based on relative abundance and store the exact sequences in the reference sequence slot - can access with refseq(ps)
dna <- Biostrings::DNAStringSet(taxa_names(phy.mts.16s))
names(dna) <- taxa_names(phy.mts.16s)
phy.mts.16s <- merge_phyloseq(phy.mts.16s, dna)
taxa_names(phy.mts.16s) <- paste0("ASV", seq(ntaxa(phy.mts.16s)))
phy.mts.16s
sampleID <- sample_names(phy.mts.16s)
raw_reads <- sample_sums(phy.mts.16s)
# Total number of bacteria sequencing reads prior to contaminant removal
saveRDS(phy.mts.16s,file.path(data_dir, "phy.mts.16s.rds"))

remove(controls, dna, metadata_16s, mts.seqtab.nochim, mts.taxa.ehomd, samples, taxtable)

##################################
# EXAMINE NEGATIVE CONTROL SAMPLES
##################################

phy.controls <- subset_samples(phy.mts.16s, control==TRUE)
phy.controls <- filter_taxa(phy.controls, function(x) sum(x) > 0, TRUE)
sample_sums(phy.controls)
ntaxa(phy.controls)
control_asvs <- data.frame(t(otu_table(phy.controls)))
control_asvs$ASV <- rownames((control_asvs))
taxtable <- data.frame(tax_table(phy.mts.16s))
taxtable$ASV <- rownames((taxtable))
control_asvs <- merge(taxtable, control_asvs, by="ASV", all.x=FALSE)
rownames(control_asvs) <- control_asvs[,1]
control_asvs <- control_asvs[,-1]
phy.controls.agglom <- tax_glom(phy.controls, taxrank = 'Genus')
ntaxa(phy.controls.agglom)
phy.controls.relative <- transform_sample_counts(phy.controls.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.controls.relative))  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
relative_controls <- psmelt(phy.controls.relative)
relative_controls$Genus <- as.character(relative_controls$Genus)
control_genera <- aggregate(relative_controls$Abundance, by=list(Phylum=relative_controls$Phylum, Genus=relative_controls$Genus,
                                                                 OTU=relative_controls$OTU), FUN=mean)
control_genera <- aggregate(control_genera$x, by=list(Phylum=control_genera$Phylum,
                                                      Genus=control_genera$Genus), FUN=sum)
colnames(control_genera)[3] <- "genus_Ab"
sum(control_genera$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(control_genera)            # Corresponds to # of unique genera
control_genera <- control_genera[order(-control_genera$genus_Ab),]
head(control_genera)

remove(phy.controls, phy.controls.agglom, phy.controls.relative, relative_controls)

#########################################
# REMOVAL OF CONTAMINANT ASVS IN DECONTAM
#########################################

set.seed(1234)
library(decontam); packageVersion("decontam")

# Identify contaminants based on prevalence in negative controls (use relative abundances)
phy.mts.16s.relative <- transform_sample_counts(phy.mts.16s, function(Abundance) Abundance/sum(Abundance))
taxtable <- data.frame(tax_table(phy.mts.16s.relative))
suppressWarnings(prev <- isContaminant(phy.mts.16s.relative, method="prevalence", neg="neg_control", batch="batch_16s", threshold=0.3))
table(prev$contaminant)
contam_prev <- subset(prev, contaminant==TRUE)
contam_prev <- merge(taxtable, contam_prev, by="row.names")
ntaxa(phy.mts.16s)
phy.mts.16s.decontam.prev <- prune_taxa(!prev$contaminant, phy.mts.16s)
ntaxa(phy.mts.16s.decontam.prev)

# Identify contaminants based on frequency (use relative abundances)
suppressWarnings(freq <- isContaminant(phy.mts.16s.decontam.prev, method="frequency", conc="dna_16s", batch="batch_16s", threshold=0.3))
table(freq$contaminant)
contam_freq <- subset(freq, contaminant==TRUE)
contam_freq <- merge(taxtable, contam_freq, by="row.names")
names(contam_freq)[names(contam_freq)=="Row.names"] <- "Taxon"
contam_freq <- contam_freq[,-c(2:6)]
freq_abundant <- freq
freq_abundant$contaminant[freq_abundant$prev<20] <- FALSE
contaminant_plot <- plot_frequency(phy.mts.16s.decontam.prev, taxa_names(phy.mts.16s.decontam.prev)[sample(which(freq_abundant$contaminant))], conc="dna_16s") + 
  xlab("DNA Concentration") + aes(color=batch_16s)
png(file="R_Plots/contaminant_plot_mts_16s.png", width = 10, height = 10, units = 'in', res = 600)
print(contaminant_plot)
dev.off()
# Remove the following ASVs from the list of contaminants after review of contaminant plots: 
# ASV118 (Lawsonella clevelandensis), ASV31 (Staphylococcus epidermidis), ASV34 (Staphylococcus epidermidis), 
# ASV132 (Staphylococcus epidermidis), ASV44 (Staphylococcus capitis/caprae), ASV62 (Staphylococcus epidermidis), 
# ASV65 (Staphylococcus epidermidis), ASV71 (Staphylococcus epidermidis)
table(freq$contaminant)
freq$contaminant[row.names(freq)=="ASV118" | row.names(freq)=="ASV31" | row.names(freq)=="ASV34" | 
                   row.names(freq)=="ASV132" | row.names(freq)=="ASV44" | row.names(freq)=="ASV62" | 
                   row.names(freq)=="ASV65" | row.names(freq)=="ASV71"] <- FALSE
table(freq$contaminant)
ntaxa(phy.mts.16s.decontam.prev)
phy.mts.16s.decontam <- prune_taxa(!freq$contaminant, phy.mts.16s.decontam.prev)
ntaxa(phy.mts.16s.decontam)

# Create phyloseq object after reagent contaminant removal
nsamples(phy.mts.16s.decontam)
sample_sums(phy.mts.16s.decontam)
metadata_16s_decontam <- data.frame(sample_data(phy.mts.16s.decontam))
metadata_16s_decontam$sampleID <- row.names(metadata_16s_decontam)
decontam_reads <- sample_sums(phy.mts.16s.decontam)
raw_reads <- data.frame(sampleID, raw_reads)
decontam_reads <- data.frame(sampleID, decontam_reads)
read_counts <- merge(raw_reads, decontam_reads, by="sampleID", all.x=F)
metadata_16s_decontam <- merge(metadata_16s_decontam, read_counts, by="sampleID", all.x=F)
rownames(metadata_16s_decontam) <- metadata_16s_decontam[,1]
sample_data(phy.mts.16s.decontam) <- metadata_16s_decontam
saveRDS(phy.mts.16s.decontam,file.path(data_dir, "phy.mts.16s.decontam.rds"))

remove(contaminant_plot, decontam_reads, freq, freq_abundant, metadata_16s_decontam, phy.mts.16s, phy.mts.16s.decontam.prev, 
       phy.mts.16s.relative, prev, raw_reads, read_counts, taxtable)

#################
# SAMPLE PRUNING
#################

phy.mts.16s.decontam <- readRDS(file.path(data_dir, "phy.mts.16s.decontam.rds"))

# Remove control samples 
nsamples(phy.mts.16s.decontam)
phy.mts.16s.decontam <- subset_samples(phy.mts.16s.decontam, control==FALSE)
nsamples(phy.mts.16s.decontam)
metadata <- data.frame(sample_data(phy.mts.16s.decontam))
nrow(metadata)
decontam_reads <- sample_sums(phy.mts.16s.decontam)
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

# Create dataframe with Shannon index and observed species for each specimen
suppressWarnings(diversity <- estimate_richness(phy.mts.16s.decontam, measures = c("Shannon", "Observed")))
setDT(diversity, keep.rownames = TRUE)[]
diversity$sampleID <- diversity$rn
diversity$rn <- NULL
metadata <- merge(metadata, diversity, by="sampleID")
row.names(metadata) <- metadata$sampleID
sample_data(phy.mts.16s.decontam) <- metadata
remove(diversity)

# Plots of alpha diversity measures by sample read counts
brks <- c(1000,2500,5000,7500,10000,15000,20000,25000,30000,40000,50000,75000,100000,150000,200000,250000)
main_brks <- c(5000,10000,25000,50000,100000,250000)
shannon <- ggplot(metadata, aes(decontam_reads, Shannon)) +
  geom_point(alpha=0.7, size=1.8) +
  geom_vline(xintercept=main_brks, linetype="dashed", linewidth=0.3) +
  scale_x_log10(breaks=brks, limits=c(1000,250000), expand=c(0.01,0)) +
  theme_bw()
png(file="R_Plots/shannon_mts_16s.png", width=12, height=6, units='in', res=600)
print(shannon)
dev.off()

observed <- ggplot(metadata, aes(decontam_reads, Observed)) +
  geom_point(alpha=0.7, size=1.8) +
  geom_vline(xintercept=main_brks, linetype="dashed", linewidth=0.3) +
  scale_x_log10(breaks=brks, limits=c(1000,250000), expand=c(0.01,0)) +
  coord_cartesian(ylim=c(0,100)) +
  theme_bw()
png(file="R_Plots/observed_mts_16s.png", width=12, height=6, units='in', res=600)
print(observed)
dev.off()
remove(shannon, observed)

# Based on these plots, use a cutoff of 25,000 reads for assessing nasal microbiome diversity
nsamples(phy.mts.16s.decontam)
summary(sample_sums(phy.mts.16s.decontam))
phy.mts.16s.pruned <- prune_samples(sample_sums(phy.mts.16s.decontam)>=25000, phy.mts.16s.decontam)
nsamples(phy.mts.16s.pruned)

# Remove taxa with no reads in any remaining samples
ntaxa(phy.mts.16s.pruned)
phy.mts.16s.pruned <- filter_taxa(phy.mts.16s.pruned, function(x) sum(x) > 0, TRUE)
nsamples(phy.mts.16s.pruned)
ntaxa(phy.mts.16s.pruned)
metadata_pruned <- data.frame(sample_data(phy.mts.16s.pruned))
sum(metadata_pruned$raw_reads)
summary(metadata_pruned$raw_reads)
summary(metadata_pruned$decontam_reads)
summary(metadata_pruned$Shannon)
summary(metadata_pruned$Observed)
saveRDS(phy.mts.16s.pruned,file.path(data_dir, "phy.mts.16s.pruned.rds"))

# Final dataset aggregated at species level with read counts
phy.mts.16s.pruned <- readRDS(file.path(data_dir, "phy.mts.16s.pruned.rds"))
ntaxa(phy.mts.16s.pruned)
phy.mts.16s.pruned.species <- collapse_phyloseq_to_species(phy.mts.16s.pruned)
tax.df <- as.data.frame(tax_table(phy.mts.16s.pruned.species))
tax.df[] <- lapply(tax.df,as.character)
no.species <- is.na(tax.df$Species) |
  tax.df$Species == "" |
  tax.df$Species == "NA" |
  tax.df$Species == "Other" |
  grepl("uncultured|unclassified|unidentified|metagenome",tax.df$Species,ignore.case=TRUE) |
  grepl(";s$",tax.df$Species) |
  grepl("^([A-Za-z]+) \\1;s$",tax.df$Species)
tax.df$Species[no.species] <- "Other"
tax.df[no.species,colnames(tax.df)!="Species"] <- NA
tax_table(phy.mts.16s.pruned.species) <- tax_table(as.matrix(tax.df))
phy.rel <- transform_sample_counts(phy.mts.16s.pruned.species,function(x) x/sum(x))
otu.rel <- as(otu_table(phy.rel),"matrix")
if(!taxa_are_rows(phy.rel)) otu.rel <- t(otu.rel)
mean.rel <- rowMeans(otu.rel)
prev <- rowMeans(otu.rel>0)
keep.taxa <- names(mean.rel)[mean.rel>0.005 & prev>0.10]
tax.df <- as.data.frame(tax_table(phy.mts.16s.pruned.species))
tax.df[] <- lapply(tax.df,as.character)
tax.df$Species[!(rownames(tax.df) %in% keep.taxa)] <- "Other"
tax.df[!(rownames(tax.df) %in% keep.taxa),colnames(tax.df)!="Species"] <- NA
tax_table(phy.mts.16s.pruned.species) <- tax_table(as.matrix(tax.df))
phy.mts.16s.pruned.species <- tax_glom(phy.mts.16s.pruned.species,taxrank="Species",NArm=FALSE)
taxa_names(phy.mts.16s.pruned.species) <- make.unique(as.character(tax_table(phy.mts.16s.pruned.species)[,"Species"]))
nsamples(phy.mts.16s.pruned.species)
ntaxa(phy.mts.16s.pruned.species)
summary(sample_sums(phy.mts.16s.pruned.species))
phy.rel.final <- transform_sample_counts(phy.mts.16s.pruned.species,function(x) x/sum(x))
otu.rel.final <- as(otu_table(phy.rel.final),"matrix")
if(!taxa_are_rows(phy.rel.final)) otu.rel.final <- t(otu.rel.final)
mean.rel.final <- rowMeans(otu.rel.final)
sort(mean.rel.final,decreasing=TRUE)
names(sort(mean.rel.final,decreasing=TRUE))
saveRDS(phy.mts.16s.pruned.species,file.path(data_dir,"phy.mts.16s.pruned.species.rds"))

# Create a dataset with "absolute" abundances obtained through relative abundances x qPCR
nsamples(phy.mts.16s.pruned.species)
phy.mts.16s.pruned.species.qpcr <- prune_samples(!is.na(sample_data(phy.mts.16s.pruned.species)$qpcr_log10),phy.mts.16s.pruned.species)
phy.mts.16s.pruned.species.qpcr <- make_qpcr_phyloseq(phy.mts.16s.pruned.species.qpcr)
nsamples(phy.mts.16s.pruned.species.qpcr)
summary(sample_sums(phy.mts.16s.pruned.species.qpcr))
saveRDS(phy.mts.16s.pruned.species.qpcr,file.path(data_dir,"phy.mts.16s.pruned.species.qpcr.rds"))