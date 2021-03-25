#************************************
# Main data preparation script for publication:
# "The upper respiratory tract microbiome of Australian Aboriginal and Torres 
# Strait Islander children in ear and nose health and disease; a prospective 
# cohort study"
# 
# Data is formatted for use in companion scripts.
#
# NOTE - Although the below script refers to each representative sequence as 
# an OTU, these sequences are actually amplicon sequence variants (ASVs).
# For more information, see : https://www.nature.com/articles/ismej2017119
#************************************

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

# Uncomment and run to install the libraries that might be needed 
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("ComplexHeatmap")
# BiocManager::install("phyloseq")
# BiocManager::install("decontam")
# BiocManager::install("DESeq2")

# install.packages("ggplot2")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("RColorBrewer")
# install.packages("vegan")
# install.packages("reshape2")
# install.packages("gplots")
# install.packages("heatmap3")
# install.packages("ggfortify")
# install.packages("seqinr")
# install.packages("psych")
# install.packages("ggnewscale")
# install.packages("cowplot")

# install.packages("devtools")
# devtools::install_github('thomasp85/ggraph')
# devtools::install_github("donaldtmcknight/microDecon")

library(decontam)
library(microDecon)
library(DESeq2)
library(ggplot2)
library(ggnewscale)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(gplots)
library(phyloseq)
library(seqinr)


####################################
# Define various colour palettes
colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
colour_palette_10 <- c("#4bafd0","#d14c38","#52a874","#b05cc6","#77b241","#7178ca","#d19b42","#c75a8e","#837f38","#c06a5a")
colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
colour_palette_25_distinct <- c("#005e27","#34adff","#352138","#f10039","#92bb00","#99e0ff","#373300","#680015","#fca000","#d0dbb1","#507400","#9900b7","#3a32cd","#ff6f21","#ff9794","#e700a6","#54032d","#ffc74c","#dfda84","#8d2700","#271472","#eecbf4","#013e77","#d59bff","#ff81ce")
colour_palette_25_distinct_b <- c("#82ff7a","#f1beff","#26324c","#c82be3","#006f91","#7a003e","#efe5ff","#ff5bbc","#9bff42","#004e1c","#ce8c00","#3a246d","#c6b700","#7c1100","#2affbc","#92008e","#9a9900","#c60024","#ff753e","#009b8c","#ffd68a","#009008","#489aff","#6e19cc","#ecffc8")
colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
# colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#a85bd2","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")

####################################

common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
  legend.key=element_blank(),
  legend.direction="vertical",
  legend.background = element_rect(colour ="white", size = .3),
  legend.text.align = 0,
  legend.title = element_text(size=10, face="bold"),
  legend.title.align = 0.5,
  legend.margin = margin(c(2,2,2,2)),
  legend.key.height=unit(.4,"cm"),
  legend.text = element_text(size = 8),
  axis.text = element_text(size = 9, colour = "black"),
  axis.title = element_text(size = 10,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))


################################################################################################
# Set the working directory

source("code/helper_functions.R")


###############################################################
# Create result directories if they are missing
dir.create(file.path(".", "Result_figures"), showWarnings = FALSE)
dir.create(file.path(".", "Result_tables"), showWarnings = FALSE)
# dir.create(file.path(".", "Result_objects"), showWarnings = FALSE)

dir.create(file.path("./Result_tables", "other"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "count_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "relative_abundance_tables"), showWarnings = FALSE)

dir.create(file.path("./Result_tables/diversity_analysis/otu"),showWarnings = FALSE, recursive = T)
dir.create(file.path("./Result_tables/diversity_analysis/genus"),showWarnings = FALSE, recursive = T)


dir.create(file.path("./Result_tables/fastspar_inputs/Nose"),showWarnings = FALSE, recursive = T)
dir.create(file.path("./Result_tables/fastspar_inputs/Community"),showWarnings = FALSE, recursive = T)
dir.create(file.path("./Result_tables/fastspar_inputs/Otitis_Status"),showWarnings = FALSE, recursive = T)
dir.create(file.path("./Result_tables/fastspar_inputs/Community__Gold_Star"),showWarnings = FALSE, recursive = T)
dir.create(file.path("./Result_tables/fastspar_inputs/Community__Otitis_Status"),showWarnings = FALSE, recursive = T)

dir.create(file.path("./Result_tables", "DESeq_results"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "abundance_analysis_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "combined_counts_abundances_and_metadata_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "stats_various"), showWarnings = FALSE)

dir.create(file.path("./Result_figures", "abundance_analysis_plots"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "pca_plots"), showWarnings = FALSE)

dir.create(file.path("./Result_figures", "diversity_analysis"), showWarnings = FALSE)

dir.create(file.path("./Result_figures/diversity_analysis", "genus"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/diversity_analysis", "genus_within_community"), showWarnings = FALSE,recursive = T)


dir.create(file.path("./Result_figures", "heatmaps"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "exploratory_analysis"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "DESeq_plots"), showWarnings = FALSE)

dir.create(file.path("./Result_figures/pca_plots", "otu"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/pca_plots", "genus"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/pca_plots", "otu_within_community"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/pca_plots", "genus_within_community"), showWarnings = FALSE,recursive = T)

dir.create(file.path("./Result_other", "sequences"), showWarnings = FALSE,recursive = T)

dir.create(file.path("./Result_figures/correlation_analysis", "networks"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/correlation_analysis/networks","genus"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/correlation_analysis/networks","otu"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/correlation_analysis", "corrplots"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_figures/correlation_analysis", "by_feature"), showWarnings = FALSE,recursive = T)

dir.create(file.path("./Result_tables/correlation_analysis", "networks"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_tables/correlation_analysis", "corrplots"), showWarnings = FALSE,recursive = T)
dir.create(file.path("./Result_tables/correlation_analysis", "by_feature"), showWarnings = FALSE,recursive = T)

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
create_combined_dataframe_no_rare <- function(counts.df, abundances.df, mymetadata, mylevel = "OTU", otu_map.df = NULL){
  counts <- clean_dataframe(counts.df)
  rel_abundances <- clean_dataframe(abundances.df)
  # Ensure ordering is the same
  rel_abundances <- rel_abundances[rownames(counts),,drop = F]
  # Combine the datasets. Passing as.matrix(counts) captures the rownames as a column. This can be renamed after
  combined_data <- cbind(melt(as.matrix(counts), variable.name = "sample", value.name = "Read_count"),
                         melt(rel_abundances, value.name = "Relative_abundance")[,2, drop = F])
  # Remove samples with a read count of zero
  combined_data <- combined_data[combined_data$Read_count > 0,]
  # Calculate logged read counts
  combined_data$Read_count_logged <- log(combined_data$Read_count, 10)
  # Fix the Var2 column
  names(combined_data)[2] <- "Sample"
  # Merge with metadata. Assumes an Index column matching Sample
  combined_data <- merge(combined_data, mymetadata, by.x = "Sample", by.y = "Index")
  if (mylevel == "OTU.ID"){
    names(combined_data)[names(combined_data) == "Var1"] <- "OTU.ID"
    combined_data <- merge(combined_data, otu_map.df, by.x = "OTU.ID", by.y = "OTU.ID")
  }
  else if (mylevel == "Species"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_species"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species", "taxonomy_species")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_species", by.y = "taxonomy_species")
  }
  else if (mylevel == "Genus"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_genus"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "taxonomy_genus")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_genus", by.y = "taxonomy_genus")
  }
  else if (mylevel == "Family"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_family"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "taxonomy_family")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_family", by.y = "taxonomy_family")
  }
  else if (mylevel == "Order"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_order"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "Order", "taxonomy_order")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_order", by.y = "taxonomy_order")
  }
  else if (mylevel == "Class"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_class"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "Class", "taxonomy_class")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_class", by.y = "taxonomy_class")
  }
  else if (mylevel == "Phylum"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_phylum"
    otu_map_reduced.df <- unique(otu_map.df[,c("Domain","Phylum", "taxonomy_phylum")])
    combined_data <- merge(combined_data, otu_map_reduced.df, by.x = "taxonomy_phylum", by.y = "taxonomy_phylum")
  }
  return(combined_data)
}
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Load the metadata
metadata.df <- read.table("data/metadata.tsv", header = T, sep = "\t")

# Change metadata name for clean sample id to index
metadata.df$Index <- metadata.df$Sequence_file_ID_clean

# Make the index the rowname
rownames(metadata.df) <- metadata.df$Index

# Make empty cells NA
metadata.df[metadata.df == ''] <- NA

# Remove unwanted columns
metadata.df <- metadata.df[,!names(metadata.df) %in% c("N_Adeno","N_WUKI","N_BOCA","N_COV_OC43","N_COV_NL63","N_HKU_1","N_ENT","N_hMPV","N_PARA_1","N_PARA_2","N_RSV_A","N_RSV_B","N_FLU_B","N_FLU_A","Virus_any")]
names(metadata.df)
# ------------------------

# Create custom variable combinations
metadata.df$Community__Season <- with(metadata.df, paste0(Community, "__", Season))
metadata.df$Community__Season[grepl("NA", metadata.df$Community__Season)] <- NA
metadata.df$Community__Season <- factor(metadata.df$Community__Season)

metadata.df$Community__Gold_Star <- with(metadata.df, paste0(Community, "__", Gold_Star))
metadata.df$Community__Gold_Star[grepl("NA", metadata.df$Community__Gold_Star)] <- NA
metadata.df$Community__Gold_Star <- factor(metadata.df$Community__Gold_Star)

metadata.df$Community__Otitis_Status <- with(metadata.df, paste0(Community, "__", Otitis_Status))
metadata.df$Community__Otitis_Status[grepl("NA", metadata.df$Community__Otitis_Status)] <- NA
metadata.df$Community__Otitis_Status <- factor(metadata.df$Community__Otitis_Status)

metadata.df$Tympanic_membrane__Gold_Star <- with(metadata.df, paste0(Tympanic_membrane, "__", Gold_Star))
metadata.df$Tympanic_membrane__Gold_Star[grepl("NA", metadata.df$Tympanic_membrane__Gold_Star)] <- NA
metadata.df$Tympanic_membrane__Gold_Star <- factor(metadata.df$Tympanic_membrane__Gold_Star)

metadata.df$Otitis_Status__Gold_Star <- with(metadata.df, paste0(Otitis_Status, "__", Gold_Star))
metadata.df$Otitis_Status__Gold_Star[grepl("NA", metadata.df$Otitis_Status__Gold_Star)] <- NA
metadata.df$Otitis_Status__Gold_Star <- factor(metadata.df$Otitis_Status__Gold_Star)


# Collapse the qPCR results into a single column
metadata.df$H.Influenzae_qPCR <- NA
metadata.df[which(metadata.df$H.Influenzae_ND == 1),]$H.Influenzae_qPCR <- "H.Influenzae_ND"
metadata.df[which(metadata.df$H.Influenzae_1st_IQR == 1),]$H.Influenzae_qPCR <- "H.Influenzae_1st_IQR"
metadata.df[which(metadata.df$H.Influenzae_2nd_to_3rd_IQR == 1),]$H.Influenzae_qPCR <- "H.Influenzae_2nd_to_3rd_IQR"
metadata.df[which(metadata.df$H.Influenzae_more_than_3rd_IQR == 1),]$H.Influenzae_qPCR <- "H.Influenzae_more_than_3rd_IQR"

metadata.df$M.catarrhalis_qPCR <- NA
metadata.df[which(metadata.df$M.catarrhalis_ND == 1),]$M.catarrhalis_qPCR <- "M.catarrhalis_ND"
metadata.df[which(metadata.df$M.catarrhalis_1st_IQR == 1),]$M.catarrhalis_qPCR <- "M.catarrhalis_1st_IQR"
metadata.df[which(metadata.df$M.catarrhalis_2nd_to_3rd_IQR == 1),]$M.catarrhalis_qPCR <- "M.catarrhalis_2nd_to_3rd_IQR"
metadata.df[which(metadata.df$M.catarrhalis_more_than_3rd_IQR == 1),]$M.catarrhalis_qPCR <- "M.catarrhalis_more_than_3rd_IQR"

metadata.df$S.pneumoniae_qPCR <- NA
metadata.df[which(metadata.df$S.pneumoniae_ND == 1),]$S.pneumoniae_qPCR <- "S.pneumoniae_ND"
metadata.df[which(metadata.df$S.pneumoniae_1st_IQR == 1),]$S.pneumoniae_qPCR <- "S.pneumoniae_1st_IQR"
metadata.df[which(metadata.df$S.pneumoniae_2nd_to_3rd_IQR == 1),]$S.pneumoniae_qPCR <- "S.pneumoniae_2nd_to_3rd_IQR"
metadata.df[which(metadata.df$S.pneumoniae_more_than_3rd_IQR == 1),]$S.pneumoniae_qPCR <- "S.pneumoniae_more_than_3rd_IQR"


# Set order of Otitis Status
metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Never OM", "HxOM", "Acute Otitis Media", "Effusion", "Perforation"))
# Set order of Season
metadata.df$Season <- factor(metadata.df$Season, levels = c("Spring", "Winter", "Autumn"))
# Set order of Nose
metadata.df$Nose <- factor(metadata.df$Nose, levels = c("Normal", "Serous", "Purulent"))
# Set order of No_peop_res_discrete
metadata.df$No_peop_res_discrete <- factor(metadata.df$No_peop_res_discrete, levels = c("2 to 3", "4 to 6", "7 to 12", "Unknown"))

# Retain copy of metadata prior to filtering samples
metadata_unfiltered.df <- metadata.df
# ------------------------------------------------------------------------------------------------------------

# Load and process the OTU table
project_otu_table.df <- read.csv("data/features_statistics_paired_noqc_250_230.csv")

# Remove spaces around " ; " in Taxon
project_otu_table.df$Taxon <- gsub("; ", ";", project_otu_table.df$Taxon)

# Fix name of first column
names(project_otu_table.df)[1] <- "OTU.ID"

# Sort by hash
project_otu_table.df <- project_otu_table.df[order(project_otu_table.df$OTU.ID),]

# Assign short form ID
project_otu_table.df$ASV_ID <- paste0("ASV_", 1:length(project_otu_table.df$OTU.ID))

# Remove J001 from sample names
names(project_otu_table.df) <- gsub("_J001", "", names(project_otu_table.df))

# Remove PGDCPos sample
project_otu_table.df[,"PGDCPos"] <- NULL

# Remove negative samples with zero counts
bad_negs <- grep("Neg",colnames(project_otu_table.df), value= T)[colSums(project_otu_table.df[grep("Neg",colnames(project_otu_table.df), value= T)]) == 0]
project_otu_table.df[,bad_negs] <- NULL

# Get the sample ids from the OTU table
sample_ids <- names(project_otu_table.df)[!names(project_otu_table.df) %in% c("OTU.ID","ASV_ID","Frequency", "Taxon", "Confidence", "RepSeq") ]
sample_ids_original <- sample_ids

# Get the negative sample IDs
negative_sample_ids <- grep("Neg",sample_ids, value= T)
not_negative_sample_ids <- sample_ids[!sample_ids %in% negative_sample_ids]

# Ensure rows with no counts are removed
project_otu_table.df <- project_otu_table.df[apply(project_otu_table.df[sample_ids],MARGIN = 1,function(z) any(z!=0)),]

# Split the Taxon column into Domain, Phylum...Species 
project_otu_table.df <- separate(project_otu_table.df, "Taxon", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"), remove =F, sep = ";")

# Splitting taxa strings that are not specified at certain taxa levels will produce NA entries at those levels. 
# NA entries should be changed to "Unassigned"
project_otu_table.df[is.na(project_otu_table.df)] <- "Unassigned"

# Create taxonomy strings for phylum, class, order, family, genus and specie levels
project_otu_table.df$taxonomy_phylum <- with(project_otu_table.df, paste(Domain, Phylum, sep =";"))
project_otu_table.df$taxonomy_class <- with(project_otu_table.df, paste(Domain, Phylum, Class, sep =";"))
project_otu_table.df$taxonomy_order <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, sep =";"))
project_otu_table.df$taxonomy_family <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, Family, sep =";"))
project_otu_table.df$taxonomy_genus <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, Family, Genus, sep =";"))
project_otu_table.df$taxonomy_species <- with(project_otu_table.df, paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep =";"))

# Store a version of the unfiltered project table
project_otu_table_unfiltered.df <- project_otu_table.df

# ------------------------------------------------
# ------------------------------------------------
# Assign unique colours for each discrete state
discrete_variables <- c("Nose","Tympanic_membrane", "Otitis_Status",
                        "Season","Community","Gold_Star","No_peop_res_discrete",
                        "H.influenzae_culture","M.catarrhalis_culture","S.pneumoniae_culture",
                        "Otitis_Status__Gold_Star", "Tympanic_membrane__Gold_Star",
                        "Community__Season","Community__Gold_Star","Community__Otitis_Status",
                        "H.Influenzae_qPCR", "M.catarrhalis_qPCR", "S.pneumoniae_qPCR",
                        "Corynebacterium_pseudodiphtheriticum","Dolosigranulum_pigrum","N_HRV")

# never - hx- "Acute Otitis Media" -effusion
otitis_status_palette <- setNames(c("#4395bf","#43bf70", "#ffc400","#ff7b00", "#b31313"),c("Never OM", "HxOM", "Acute Otitis Media", "Effusion", "Perforation"))
all_variable_colours <- as.character(lapply(as.character(metadata.df[,"Otitis_Status"]), function(x) otitis_status_palette[x]))
metadata.df[,"Otitis_Status_colour"] <- all_variable_colours

#c49d6c
# c43e00
community_colours <- c("Remote" = "#c49d6c", "Rural" = "#8ea100")
all_community_colours <- as.character(lapply(as.character(metadata.df$Community), function(x) community_colours[x]))
metadata.df$Community_colour <- all_community_colours

nose_colours <- c("Normal" = "#79aec9", "Serous" = "#ffec45", "Purulent" = "#8cbd24")
all_nose_colours <- as.character(lapply(as.character(metadata.df$Nose), function(x) nose_colours[x]))
metadata.df$Nose_colour <- all_nose_colours

season_colours <- c("Spring" = "#549627", "Winter" = "#bae1ff", "Autumn" = "#de8900")
all_season_colours <- as.character(lapply(as.character(metadata.df$Season), function(x) season_colours[x]))
metadata.df$Season_colour <- all_season_colours


for (myvar in discrete_variables){
  if (!myvar %in% c("Community", "Season", "Otitis_Status", "Nose", "Season", "No_peop_res_discrete")){
    myvar_values <- factor(as.character(sort(unique(metadata.df[,myvar]))))  
  }else{
    myvar_values <- unique(metadata.df[,myvar])
  }
  if (length(myvar_values) <= 8){
    myvar_colours <- setNames(brewer.pal(n = 8, name = "Dark2")[1:length(myvar_values)], myvar_values)  
  } else{
    myvar_colours <- setNames(colour_palette_10[1:length(myvar_values)], myvar_values)  
  }
  
  if (!paste0(myvar,"_colour") %in% colnames(metadata.df)){
    all_variable_colours <- as.character(lapply(as.character(metadata.df[,myvar]), function(x) myvar_colours[x]))
    metadata.df[,paste0(myvar,"_colour")] <- all_variable_colours    
  }
}
metadata.df[with(metadata.df, !is.na(No_peop_res_discrete) & No_peop_res_discrete == "Unknown"),]$No_peop_res_discrete_colour <- "grey"

# ------------------------------------------------
#         Remove unwanted lineages
fungal_phyla <- c("Basidiomycota","Ascomycota","Mucoromycota","Cryptomycota","Peronosporomycetes","Blastocladiomycota","Chytridiomycota","Zoopagomycota","Neocallimastigomycota")
#"Aphelidea","LKM15" - unclear if these are fungal
fungal_pattern <- paste0(fungal_phyla,collapse = "|p__")
fungal_pattern <- gsub("^", "p__", fungal_pattern)
#  ----------------------
# Rough stats on the discarded features
project_otu_table_discarded.df <- project_otu_table.df[!grepl(paste0("d__Bacteria|d__Archaea|",fungal_pattern), project_otu_table.df$Taxon) | 
                                                         project_otu_table.df$Phylum == "Unassigned" | 
                                                         grepl("o__Chloroplast", project_otu_table.df$Taxon) |
                                                         grepl("f__Mitochondria", project_otu_table.df$Taxon) ,]

# Abundance of discarded features in each sample based on unfiltered data
temp <- melt(round(colSums(project_otu_table_discarded.df[sample_ids]) / colSums(project_otu_table_unfiltered.df[sample_ids]) * 100,4), value.name = "Relative_abundance")
temp <- m2df(temp,column_name = "Index")
temp <- merge(temp,metadata.df,by = "Index") %>% select(Index, Community, Gold_Star, Relative_abundance)
temp$Relative_abundance[is.na(temp$Relative_abundance)] <- 0
temp <- temp[order(temp$Relative_abundance,decreasing = T),]
max(temp$Relative_abundance)
min(temp$Relative_abundance)
mean(temp$Relative_abundance)
median(temp$Relative_abundance)

# Read counts of discarded features in each sample based on unfiltered data
temp <- melt(colSums(project_otu_table_discarded.df[sample_ids]), value.name = "Read_count")
temp <- m2df(temp,column_name = "Index")
temp <- merge(temp,metadata.df,by = "Index") %>% select(Index, Community, Gold_Star, Read_count)
temp <- temp[order(temp$Read_count,decreasing = T),]
max(temp$Read_count)
min(temp$Read_count)
mean(temp$Read_count)
median(temp$Read_count)
# Number of features discarded
length(project_otu_table_discarded.df$OTU.ID) / length(project_otu_table_unfiltered.df$OTU.ID) *100

# Abundance of Unassigned features in each sample based on unfiltered data
dim(project_otu_table_discarded.df[project_otu_table_discarded.df$Taxon == "Unassigned",])
temp <- melt(round(colSums(project_otu_table_discarded.df[project_otu_table_discarded.df$Taxon == "Unassigned",sample_ids])/colSums(project_otu_table_unfiltered.df[sample_ids]) * 100,4), value.name = "Relative_abundance")
temp <- m2df(temp,column_name = "Index")
temp <- merge(temp,metadata.df,by = "Index") %>% select(Index, Community, Gold_Star, Relative_abundance)
temp$Relative_abundance[is.na(temp$Relative_abundance)] <- 0
temp <- temp[order(temp$Relative_abundance,decreasing = T),]
max(temp$Relative_abundance)
min(temp$Relative_abundance)
mean(temp$Relative_abundance)
median(temp$Relative_abundance)
#  ----------------------

# Remove OTUs that are Unassigned
# project_otu_table.df <- project_otu_table.df[project_otu_table.df$Taxon != "Unassigned",]

# Discard anything not Bacterial or Archaeal or fungal
project_otu_table.df <- project_otu_table.df[grepl(paste0("d__Bacteria|d__Archaea|",fungal_pattern), project_otu_table.df$Taxon),]

# Discard anything that is Unassigned at the Phylum level
project_otu_table.df <- project_otu_table.df[!project_otu_table.df$Phylum == "Unassigned",]

# Discard chloroplast features
project_otu_table.df <- project_otu_table.df[!grepl("o__Chloroplast", project_otu_table.df$Taxon,ignore.case = T),]

# Discard mitochondria features
project_otu_table.df <- project_otu_table.df[!grepl("f__Mitochondria", project_otu_table.df$Taxon,ignore.case = T),]

# dim(project_otu_table.df)
# dim(project_otu_table_unfiltered.df)
#  ----------------------------------------------------------------------------------------------------------------
#  ----------------------------------------------------------------------------------------------------------------

# Remove old Taxon column
project_otu_table.df$Taxon <- NULL
project_otu_table_unfiltered.df$Taxon <- NULL

# Store the OTUs and corresponding taxonomy information in a separate dataframe
otu_taxonomy_map.df <- project_otu_table.df[c("OTU.ID",
                                              "ASV_ID",
                                           "Domain", 
                                           "Phylum", 
                                           "Class", 
                                           "Order", 
                                           "Family",
                                           "Genus",
                                           "Species",
                                           "taxonomy_phylum",
                                           "taxonomy_class",
                                           "taxonomy_order",
                                           "taxonomy_family",
                                           "taxonomy_genus",
                                           "taxonomy_species", 
                                           "RepSeq")]

rownames(otu_taxonomy_map.df) <- otu_taxonomy_map.df$OTU.ID

# Version without the RepSeq
reduced_tax_map <- otu_taxonomy_map.df
reduced_tax_map$RepSeq <- NULL

# Save this OTU taxonomy map for later use
write.table(otu_taxonomy_map.df, file = "Result_tables/other/otu_taxonomy_map.csv", sep = ",", quote = F, row.names = F)

# Also save the unfiltered table, to avoid processing the original data table again 
write.table(project_otu_table_unfiltered.df, file = "Result_tables/other/project_otu_table_unfiltered.csv", sep = ",", quote = F, row.names = F)

# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# Now we can generate the tables that we will need for different analyses at both the OTU and various taxa levels

# -----------------------------
# -----------OTU LEVEL---------
# Dataframe containing the counts for each OTU, for each sample
otu.df <- project_otu_table.df[c("OTU.ID", sample_ids)]
otu_unfiltered.df <- project_otu_table_unfiltered.df[c("OTU.ID", sample_ids)]

## Create a matrix version for ease of processing
# For filtered OTU matrix
otu.m <- df2matrix(otu.df)

temp <- read_counts_and_unique_features(otu.m, sort(c(sample_ids_original, "PNeg1P1")))
df2matrix(temp$sample_read_counts)
df2matrix(temp$sample_feature_counts)

# And unfiltered OTU matrix
otu_unfiltered.m <- df2matrix(otu_unfiltered.df)

temp <- read_counts_and_unique_features(otu_unfiltered.m, sample_ids_original)
df2matrix(temp$sample_read_counts)
df2matrix(temp$sample_feature_counts)

# Create relative abundance matrix from counts matrix
otu_rel.m <- t(t(otu.m)/ colSums(otu.m))
otu_unfiltered_rel.m <- t(t(otu_unfiltered.m) / colSums(otu_unfiltered.m))

# Change nans to 0. Occurs when a sample has no hits at this point.
otu_rel.m[is.nan(otu_rel.m)] <- 0
otu_unfiltered_rel.m[is.nan(otu_unfiltered_rel.m)] <- 0
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# 
# Create feature table for negative samples
negative_otu.m <- otu.m[,negative_sample_ids]
negative_otu_rel.m <- otu_rel.m[,negative_sample_ids]

# Create dataframes for negative feature tables
negative_otu.df <- m2df(negative_otu.m, "OTU.ID")
negative_otu_rel.df <- m2df(negative_otu_rel.m, "OTU.ID")

# Generate and write the final counts and abundances for each taxonomy level to file
# Merge the otu counts back with the taxonomy data
otu_metadata_merged.df <- merge(negative_otu.df, otu_taxonomy_map.df, by.x = "OTU.ID", by.y = "OTU.ID")

tax_data.l <- generate_tax_level_data(otu_metadata_merged.df, 
                                      tax_string_levels = c("taxonomy_species", "taxonomy_genus", "taxonomy_family", "taxonomy_class", "taxonomy_order", "taxonomy_phylum"),
                                      sample_ids = negative_sample_ids)

species.df <- m2df(tax_data.l$taxonomy_species$counts,"taxonomy_species")
species_rel.df <-  m2df(tax_data.l$taxonomy_species$abundances,"taxonomy_species")
genus.df <- m2df(tax_data.l$taxonomy_genus$counts,"taxonomy_genus")
genus_rel.df <-  m2df(tax_data.l$taxonomy_genus$abundances,"taxonomy_genus")
family.df <- m2df(tax_data.l$taxonomy_family$counts,"taxonomy_family")
family_rel.df <-  m2df(tax_data.l$taxonomy_family$abundances,"taxonomy_family")
order.df <- m2df(tax_data.l$taxonomy_order$counts,"taxonomy_order")
order_rel.df <-  m2df(tax_data.l$taxonomy_order$abundances,"taxonomy_order")
class.df <- m2df(tax_data.l$taxonomy_class$counts,"taxonomy_class")
class_rel.df <-  m2df(tax_data.l$taxonomy_class$abundances,"taxonomy_class")
phylum.df <- m2df(tax_data.l$taxonomy_phylum$counts,"taxonomy_phylum")
phylum_rel.df <-  m2df(tax_data.l$taxonomy_phylum$abundances,"taxonomy_phylum")

write.table(species.df, file = "Result_tables/count_tables/Negative_samples_Specie_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(species_rel.df, file = "Result_tables/relative_abundance_tables/Negative_samples_Specie_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus.df, file = "Result_tables/count_tables/Negative_samples_Genus_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus_rel.df, file = "Result_tables/relative_abundance_tables/Negative_samples_Genus_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family.df, file = "Result_tables/count_tables/Negative_samples_Family_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family_rel.df, file = "Result_tables/relative_abundance_tables/Negative_samples_Family_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order.df, file = "Result_tables/count_tables/Negative_samples_Order_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order_rel.df, file = "Result_tables/relative_abundance_tables/Negative_samples_Order_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class.df, file = "Result_tables/count_tables/Negative_samples_Class_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class_rel.df, file = "Result_tables/relative_abundance_tables/Negative_samples_Class_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum.df, file = "Result_tables/count_tables/Negative_samples_Phylum_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum_rel.df, file = "Result_tables/relative_abundance_tables/Negative_samples_Phylum_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)


# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                     CONTAMINANT IDENTIFICATION AND REMOVAL

# One simple approach to removing contaminants, or 'noise', is to simply remove low abundance features or
# remove features that are present / more prevalent in the negative controls.
# Alternatively (or in addition), packages like microdecon and decontam can be used to identity putative contaminants
# We are going to try multiple approaches.

# 1. Use Decontam to identify contaminants
decontam_contaminants.df <- isContaminant(t(otu.m), 
                                          method = "prevalence", 
                                          neg = rownames(t(otu.m)) %in% negative_sample_ids, 
                                          threshold = 0.5)
hist(decontam_contaminants.df$p,100, ylim = c(0,400), xlim = c(0,1))
abline(v = 0.5, col = 'red')
decontam_contaminants_features <- rownames(subset(decontam_contaminants.df, contaminant == T))

# 2. Use microdecon
# OTU Neg1 ....Sample 1....taxa (optional)
microdecon_data.df <- m2df(otu.m[,colnames(otu.m) %in% negative_sample_ids], column_name = "OTU.ID")
microdecon_data.df <- cbind(microdecon_data.df, otu.m[,!colnames(otu.m) %in% negative_sample_ids])
rownames(microdecon_data.df) <- NULL
microdecon_data.df$taxonomy_species <- subset(otu_taxonomy_map.df, OTU.ID %in% microdecon_data.df$OTU.ID)$taxonomy_species

microdecon_contaminants.df <- decon(microdecon_data.df, 
                                    numb.blanks = length(negative_sample_ids),
                                    numb.ind = length(sample_ids) - length(negative_sample_ids),
                                    taxa = T,
                                    runs = 2)

microdecon_contaminants_features <- microdecon_contaminants.df$OTUs.removed$OTU.ID
# summary(microdecon_contaminants_features %in% decontam_contaminants_features)
# summary(decontam_contaminants_features %in% microdecon_contaminants_features)

# 3. Identify features that are present in negative controls (extreme approach) 
neg_present_features <- rownames(otu.m[,colnames(otu.m) %in% negative_sample_ids][rowSums(otu.m[,colnames(otu.m) %in% negative_sample_ids]) > 0,])

# subset(otu_taxonomy_map.df, OTU.ID %in% neg_present_features)$taxonomy_species

# 4. Identify features that are more prevalent in negative controls. This is not that useful for this study 
# as we only have a small number of negative samples.
otu_rare_count.m <- t(rrarefy(x = t(otu.m), sample=2000))
negative_otu_rare_count.m <- t(rrarefy(x = t(negative_otu.m), sample=2000))
otu_negative_sample_prevalences <- apply(negative_otu_rare_count.m, 1, function(x) {length(which(x > 0))}) /length(negative_sample_ids)
not_negative_sample_ids <- sample_ids[!sample_ids %in% negative_sample_ids]
otu_not_negative_sample_prevalences <- apply(otu_rare_count.m[,not_negative_sample_ids], 1, function(x) {length(which(x > 0))}) /length(not_negative_sample_ids)
contaminating_otus_from_prevalences <- names(otu_negative_sample_prevalences[otu_negative_sample_prevalences > otu_not_negative_sample_prevalences])

# Final contaminating features to remove
contaminating_features <- unique(c(decontam_contaminants_features, microdecon_contaminants_features))

# -----------------------
# Write contaminant taxonomy data to file
contaminating_feature_taxonomy_data.df <- otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% contaminating_features,]
contaminating_feature_taxonomy_data.df <- unique(contaminating_feature_taxonomy_data.df)
contaminating_feature_taxonomy_data.df$Decontam_contaminant <- contaminating_feature_taxonomy_data.df$OTU.ID %in% decontam_contaminants_features
contaminating_feature_taxonomy_data.df$Microdecon_contaminant <- contaminating_feature_taxonomy_data.df$OTU.ID %in% microdecon_contaminants_features
temp <- contaminating_feature_taxonomy_data.df$RepSeq
contaminating_feature_taxonomy_data.df$RepSeq <- NULL
contaminating_feature_taxonomy_data.df$RepSeq <- temp
write.csv(x = contaminating_feature_taxonomy_data.df, file = "Result_tables/other/contaminate_taxonomy_data.csv", quote =F, row.names = F)

# Write putative contaminant sequences to file
write.fasta(sequences = as.list(contaminating_feature_taxonomy_data.df$RepSeq),open = "w",
            names = as.character(contaminating_feature_taxonomy_data.df$OTU.ID),
            file.out = paste0("Result_other/sequences/contaminate_features.fasta"))

# Remove contaminants from data. Retain both datasets for comparison.
otu_rel_decontaminated.m <- otu_rel.m[!rownames(otu_rel.m) %in% unique(c(contaminating_features)),]
otu_decontaminated.m <- otu.m[!rownames(otu.m) %in% unique(c(contaminating_features)),]

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# Create contaminant count and abundance tables and write to file
otu_contaminates.m <- otu.m[rownames(otu.m) %in% contaminating_features,]
otu_rel_contaminates.m <- otu_rel.m[rownames(otu_rel.m) %in% contaminating_features,]

write.table(otu_contaminates.m, file = "Result_tables/count_tables/Contaminates_OTU_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(otu_rel_contaminates.m, file = "Result_tables/relative_abundance_tables/Contaminates_OTU_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)

# Generate and write the final counts and abundances for each taxonomy level to file
# Do this for all samples, regardless of cohort
otu_contaminates.df <- m2df(otu_contaminates.m, "OTU.ID")
otu_rel_contaminates.df <- m2df(otu_rel_contaminates.m, "OTU.ID")

# Merge the otu counts back with the taxonomy data
otu_metadata_merged.df <- merge(otu_contaminates.df, otu_taxonomy_map.df, by.x = "OTU.ID", by.y = "OTU.ID")
otu_rel_metadata_merged.df <- merge(otu_rel_contaminates.df, otu_taxonomy_map.df, by.x = "OTU.ID", by.y = "OTU.ID")


# Just use for counts
tax_data_for_counts.l <- generate_tax_level_data(otu_metadata_merged.df, 
                                                 tax_string_levels = c("taxonomy_species",
                                                                       "taxonomy_genus", 
                                                                       "taxonomy_family", 
                                                                       "taxonomy_class",
                                                                       "taxonomy_order", 
                                                                       "taxonomy_phylum"),
                                                 sample_ids = sample_ids)

# Use generate_tax_level_data function to combine at each taxonomy level, will be under 'counts'
tax_data_for_abundances.l <- 
  generate_tax_level_data(otu_rel_metadata_merged.df, 
                          tax_string_levels = c("taxonomy_species",
                                                "taxonomy_genus",
                                                "taxonomy_family",
                                                "taxonomy_class",
                                                "taxonomy_order",
                                                "taxonomy_phylum"),
                          sample_ids = sample_ids)

species.df <- m2df(tax_data_for_counts.l$taxonomy_species$counts,"taxonomy_species")
species_rel.df <-  m2df(tax_data_for_abundances.l$taxonomy_species$abundances,"taxonomy_species")
genus.df <- m2df(tax_data_for_counts.l$taxonomy_genus$counts,"taxonomy_genus")
genus_rel.df <-  m2df(tax_data_for_abundances.l$taxonomy_genus$abundances,"taxonomy_genus")
family.df <- m2df(tax_data_for_counts.l$taxonomy_family$counts,"taxonomy_family")
family_rel.df <-  m2df(tax_data_for_abundances.l$taxonomy_family$abundances,"taxonomy_family")
order.df <- m2df(tax_data_for_counts.l$taxonomy_order$counts,"taxonomy_order")
order_rel.df <-  m2df(tax_data_for_abundances.l$taxonomy_order$abundances,"taxonomy_order")
class.df <- m2df(tax_data_for_counts.l$taxonomy_class$counts,"taxonomy_class")
class_rel.df <-  m2df(tax_data_for_abundances.l$taxonomy_class$abundances,"taxonomy_class")
phylum.df <- m2df(tax_data_for_counts.l$taxonomy_phylum$counts,"taxonomy_phylum")
phylum_rel.df <-  m2df(tax_data_for_abundances.l$taxonomy_phylum$abundances,"taxonomy_phylum")


write.table(species.df, file = "Result_tables/count_tables/Contaminants_Specie_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(species_rel.df, file = "Result_tables/relative_abundance_tables/Contaminants_Specie_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus.df, file = "Result_tables/count_tables/Contaminants_Genus_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus_rel.df, file = "Result_tables/relative_abundance_tables/Contaminants_Genus_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family.df, file = "Result_tables/count_tables/Contaminants_Family_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family_rel.df, file = "Result_tables/relative_abundance_tables/Contaminants_Family_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order.df, file = "Result_tables/count_tables/Contaminants_Order_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order_rel.df, file = "Result_tables/relative_abundance_tables/Contaminants_Order_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class.df, file = "Result_tables/count_tables/Contaminants_Class_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class_rel.df, file = "Result_tables/relative_abundance_tables/Contaminants_Class_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum.df, file = "Result_tables/count_tables/Contaminants_Phylum_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum_rel.df, file = "Result_tables/relative_abundance_tables/Contaminants_Phylum_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# ----------------------------------------
# ---- General ---- 
# How many / what samples are the contaminants in?
# How abundant are the contaminants in each sample?
# What samples types are the contaminants in?
sample_contamination.df <- m2df(melt(round(colSums(otu_rel.m[contaminating_features,])*100,5),value.name = "Contaminant_abundance"),"Sample")
sample_contamination.df <- left_join(sample_contamination.df,m2df(melt(round(colSums(otu_unfiltered_rel.m[contaminating_features,])*100,5),value.name = "Contaminant_abundance_unfiltered"),"Sample"))
rownames(sample_contamination.df) <- sample_contamination.df$Sample

# Get the number of contaminant features in each sample
temp <- otu_rel.m[contaminating_features,sample_contamination.df$Sample]
temp[temp > 0] <- 1
temp <- melt(colSums(temp))
sample_contamination.df$Number_of_contaminant_features <- temp[rownames(sample_contamination.df),]

temp <- otu_rel.m[,sample_contamination.df$Sample]
temp[temp > 0] <- 1
temp <- melt(colSums(temp))
sample_contamination.df$Number_of_features <- temp[rownames(sample_contamination.df),]

# Get the number of features in each sample (unfiltered)
temp <- otu_unfiltered_rel.m[,sample_contamination.df$Sample]
temp[temp > 0] <- 1
temp <- melt(colSums(temp))
sample_contamination.df$Number_of_features_unfiltered <- temp[rownames(sample_contamination.df),]

# Fix entries where if there are no features total, contaminant abundance must be 0
sample_contamination.df[which(sample_contamination.df$Number_of_features == 0),]$Contaminant_abundance <- 0

# Order by contamination
sample_contamination.df <- sample_contamination.df[order(sample_contamination.df$Contaminant_abundance,decreasing = T),]

# Add some metadata
sample_contamination.df$Community <- metadata.df[rownames(sample_contamination.df),]$Community
sample_contamination.df$Gold_Star <- metadata.df[rownames(sample_contamination.df),]$Gold_Star

# Re-normalise the data
sample_contamanination.df <- m2df(melt(round(100-colSums(otu_rel.m)*100,5),value.name = "Percent_contaminant"),"Sample")
sample_contamanination.df <- sample_contamanination.df[order(sample_contamanination.df$Percent_contaminant),]
sample_contamanination.df
write.csv(sample_contamanination.df, "Result_tables/other/sample_contaminant_percentages.csv", quote =F, row.names = F)


min(sample_contamination.df$Contaminant_abundance)
max(sample_contamination.df$Contaminant_abundance)
mean(sample_contamination.df$Contaminant_abundance)
median(sample_contamination.df$Contaminant_abundance)
min(sample_contamination.df$Contaminant_abundance_unfiltered)
max(sample_contamination.df$Contaminant_abundance_unfiltered)
mean(sample_contamination.df$Contaminant_abundance_unfiltered)
median(sample_contamination.df$Contaminant_abundance_unfiltered)

write.csv(sample_contamination.df, "Result_tables/other/sample_contaminant_summary.csv", quote =F, row.names = F)


# Re-normalise the relative abundance data
otu_rel_decontaminated.m <- t(t(otu_rel_decontaminated.m) / colSums(otu_rel_decontaminated.m))
otu_rel_decontaminated.m[is.nan(otu_rel_decontaminated.m)] <- 0

# (Optional) Use de-contaminated data going forward
otu_rel.m <- otu_rel_decontaminated.m
otu.m <- otu_decontaminated.m

temp <- read_counts_and_unique_features(otu.m, sort(c(sample_ids_original, "PNeg1P1")))
df2matrix(temp$sample_read_counts)
df2matrix(temp$sample_feature_counts)

# (Optional) Remove negative samples as they are no longer needed
otu_rel.m <- otu_rel.m[,!colnames(otu_rel.m) %in% negative_sample_ids]
otu.m <- otu.m[,!colnames(otu.m) %in% negative_sample_ids]
negative_metadata.df <- metadata.df[rownames(metadata.df) %in% negative_sample_ids,]
metadata.df <- metadata.df[!rownames(metadata.df) %in% negative_sample_ids,]

# Reassign the sample ids
sample_ids <- colnames(otu.m)
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# Filter those OTUs that are low abundance in all samples

# Check how many OTUs would be removed if we filtered any whose abundance is less than 0.05% (0.0005) in all NON-NEGATIVE samples
filter_percentage = 0.0005
otu_rel_low_abundance_otus.m <- otu_rel.m[apply(otu_rel.m[,not_negative_sample_ids],1,function(z) all(z < filter_percentage)),]

otu_rel_low_abundance_otus.df <- as.data.frame(otu_rel_low_abundance_otus.m)
percent_removed_before_taxa_filtered <- round(dim(otu_rel_low_abundance_otus.m)[1]/length(rownames(otu_unfiltered.m)) * 100,2)
percent_removed_after_taxa_filtered <- round(dim(otu_rel_low_abundance_otus.m)[1]/length(rownames(otu.m)) * 100,2)

print(paste("There are a total of", dim(otu_rel.m)[1], "OTUs before filtering"))
print(paste("A total of", dim(otu_rel_low_abundance_otus.m)[1], 
            "OTUs will be filtered at a threshold of", 
            filter_percentage * 100, "percent"))

# If you are happy with the filtering threshold, apply it.
# Do this for each dataset (with and without contaminants)
otu_rel.m <- otu_rel.m[apply(otu_rel.m[,sample_ids],1,function(z) any(z>=filter_percentage)),]

# Re-normalise the matrix after filtering
otu_rel.m <- t(t(otu_rel.m) / colSums(otu_rel.m))

# Change nans to 0. Occurs when a sample has no hits at this point.
otu_rel.m[is.nan(otu_rel.m)] <- 0

# Also remove low abundance OTUs from the original OTU count matrix
otu.m  <- otu.m[rownames(otu_rel.m),]
temp <- read_counts_and_unique_features(otu.m, sort(c(sample_ids_original, "PNeg1P1")))
df2matrix(temp$sample_read_counts)
df2matrix(temp$sample_feature_counts)

# Discard samples with less than # reads.
count_threshold <- 2000
summary(colSums(otu.m) < count_threshold)
table(metadata_unfiltered.df[colnames(otu.m[,colSums(otu.m) < count_threshold]),]$Otitis_Status)

otu_prior_to_removing_low_read_count_samples.m <- otu.m
otu_prior_to_removing_low_read_count_samples_rel.m <- otu_rel.m

dim(otu.m)
otu.m <- otu.m[,colSums(otu.m) >= 2000]
dim(otu.m)

# The might be many rows whose maximum is 0 at this point. Remove them.
dim(otu.m)
otu.m <- otu.m[apply(otu.m, 1, max) != 0,]
dim(otu.m)

# Reassign the sample ids
sample_ids <- colnames(otu.m)

temp <- read_counts_and_unique_features(otu.m, sort(c(sample_ids_original, "PNeg1P1")))
df2matrix(temp$sample_read_counts)
df2matrix(temp$sample_feature_counts)

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# Get the most abundant unassigned features

# Project table with just unassigned features
unassigned_project_otu_table_unfiltered.df <- project_otu_table_unfiltered.df[project_otu_table_unfiltered.df$Domain == "Unassigned",]
rownames(unassigned_project_otu_table_unfiltered.df) <- unassigned_project_otu_table_unfiltered.df$OTU.ID

# Convert to dataframe
unassigned_otu_unfiltered_rel.df <- m2df(otu_unfiltered_rel.m[as.character(unassigned_project_otu_table_unfiltered.df$OTU.ID),,drop =F],
                                         column_name = "OTU.ID")

# Melt
unassigned_otu_unfiltered_rel.df <- melt(unassigned_otu_unfiltered_rel.df, variable.name = "Sample", value.name = "Relative_abundance")

if (max(unassigned_otu_unfiltered_rel.df$Relative_abundance) != 0){
  
  # Get the top most abundant unassigned features per sample
  most_abundant_unassigned.df <- unassigned_otu_unfiltered_rel.df %>% 
    group_by(Sample) %>%
    filter(Relative_abundance > 0) %>%
    top_n(n = 10, wt = Relative_abundance) %>% 
    mutate(Relative_abundance = round(Relative_abundance*100,3)) %>%
    as.data.frame() 
  
  # Get corresponding representative sequence
  most_abundant_unassigned.df$RepSeq <- unlist(lapply(most_abundant_unassigned.df$OTU.ID, function(x) as.character(unassigned_project_otu_table_unfiltered.df[unassigned_project_otu_table_unfiltered.df$OTU.ID == x,]$RepSeq)))
  write.csv(x = most_abundant_unassigned.df, file = paste0("Result_tables/other/most_abundant_unassigned.csv"), row.names = F)
  
  # Get unique set of features
  unique_most_abundant_unassigned.df <- unique(most_abundant_unassigned.df[c("OTU.ID", "RepSeq")])
  
  # Write fasta file
  write.fasta(sequences = as.list(unique_most_abundant_unassigned.df$RepSeq),open = "w", 
              names = as.character(unique_most_abundant_unassigned.df$OTU.ID),
              file.out = paste0("Result_other/sequences/most_abundant_unassigned_features.fasta"))
  
}


# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
#         Rarefying (capping maximum reads)

# Note - For samples with a read count lower then the sample=# parameter, 
# the rrarefy function in vegan will return the this sample with its existing read count.
# Sample with counts higher than the sample parameter will be rarefied as normal and their counts will be capped at the sample parameter value

# Normally samples with OTU counts lower than our desired threshold need to be removed if we want to rarefy.
# This should have been performed at this point in normal circumstances. However, we are simply using rrarefy to cap the read depth so the fold difference
# between the highest and lowest is not extreme.

# Looking at rarefaction curve and the counts for each sample can be a good
# way to determine a good minimum count threshold for samples

# Counts for each Sample
column_sums <- colSums(otu.m)
column_sums.df <- melt(column_sums[order(column_sums)])
column_sums.df$sample <- rownames(column_sums.df)
rownames(column_sums.df) <- NULL
column_sums.df <- column_sums.df[c("sample", "value")]
column_sums.df$sample <- factor(column_sums.df$sample, levels = column_sums.df$sample)

myplot <- ggplot(column_sums.df, aes(x = sample, y = value)) + 
  geom_histogram(stat = "identity") +
  geom_hline(yintercept = 30000, color = 'red')+
  geom_hline(yintercept = 20000, color = 'red')+
  geom_hline(yintercept = mean(column_sums.df$value), color = 'blue')+
  geom_hline(yintercept = median(column_sums.df$value), color = 'purple')+
  scale_y_continuous(breaks = seq(0,max(column_sums.df$value), 5000)) +
  xlab("Sample") +
  ylab("Read count") +
  common_theme +
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
ggsave(plot = myplot, filename = "./Result_figures/exploratory_analysis/sample_read_depth_distribution.pdf", width=30, height=6)

myplot <- ggplot(column_sums.df, aes(x = value)) + 
  xlab("Read count") +
  ylab("Number of samples") +
  scale_x_continuous(breaks = seq(500,max(column_sums.df$value), 2000)) +
  scale_y_continuous(breaks = seq(0,100, 10)) +
  geom_histogram(stat = "bin", bins = 80, colour = "black",fill = "grey") +
  common_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))
myplot
ggsave(plot = myplot, filename = "./Result_figures/exploratory_analysis/sample_read_depth_distribution_2.pdf", width=20, height=6)

# Rarefaction curve
# rarecurve(t(otu.m[,colSums(otu.m) > 1]),step = 500, label = F,xlim = c(0,30000),sample = 10000)

# Generate a rarefied OTU count matrix. The rrarefy function just caps samples to a maximum of threshold. Samples with less will still have less.
set.seed(1234)
otu_rare_count.m <- t(rrarefy(x = t(otu.m), sample=50000))

# How many samples have reads removed
summary(colSums(otu_rare_count.m - otu.m) < 0)
# otu_rare_count.m <- otu_rare_count.m[apply(otu_rare_count.m, 1, max) > 0,]

# Read stats on rarefied data
rare_stats.df <- data.frame(row.names = colnames(otu.m))
rare_stats.df[,'Unfiltered_data_reads'] <- colSums(otu_unfiltered.m[,colnames(otu.m)]) # reads in unfiltered data
rare_stats.df[,'Filtered_data_reads'] <- colSums(otu.m) # Reads in filtered data (contaminants, low abundance etc.)
rare_stats.df[,'Rarefied_data_reads'] <- colSums(otu_rare_count.m[,colnames(otu.m)]) # Reads in rarefied data (from filtered data)
# rare_stats.df[,'Reads_removed_unfiltered_rarefied'] <-  rare_stats.df[,'Unfiltered_reads'] - rare_stats.df[,'Rarefied_reads'] # Reads removed in rarefied vs unfiltered
rare_stats.df[,'Reads_removed_from_filtered_data_after_rarefaction'] <-  rare_stats.df[,'Filtered_data_reads'] - rare_stats.df[,'Rarefied_data_reads'] # Reads removed in rarefied vs filtered
# rare_stats.df[,'Proportion_reads_removed_unfiltered'] <- rare_stats.df[,'Reads_removed_unfiltered_rarefied'] / rare_stats.df[,'Unfiltered_reads']
rare_stats.df[,'Proportion_reads_removed_from_filtered_data_after_rarefaction'] <- rare_stats.df[,'Reads_removed_from_filtered_data_after_rarefaction'] / rare_stats.df[,'Filtered_data_reads']
rare_stats.df[,"Features_unfiltered_data_total"] <- length(which(rowSums(otu.m[,colnames(otu.m)]) > 0))
rare_stats.df[,'Features_unfiltered_data'] <- apply(otu_unfiltered.m[,colnames(otu.m)], MARGIN = 2, FUN = function(x) { length(which(x > 0)) } )
rare_stats.df[,'Features_filtered_data'] <- apply(otu.m, MARGIN = 2, FUN = function(x) { length(which(x > 0)) } )
rare_stats.df[,'Features_filtered_data_after_rarefaction'] <- apply(otu_rare_count.m, MARGIN = 2, FUN = function(x) { length(which(x > 0)) } )

rare_stats.df[,'Features_removed_from_filtered_data_after_rarefaction'] <- rare_stats.df[,'Features_filtered_data'] - rare_stats.df[,'Features_filtered_data_after_rarefaction']
rare_stats.df[,'Proportion_features_removed_from_filtered_data_after_rarefaction'] <-  rare_stats.df[,'Features_removed_from_filtered_data_after_rarefaction'] / rare_stats.df[,'Features_filtered_data']

rare_stats.df <- m2df(rare_stats.df, "Sample")
rare_stats.df$Community <- metadata.df[rare_stats.df$Sample,]$Community
rare_stats.df$Gold_Star <- metadata.df[rare_stats.df$Sample,]$Gold_Star

# Number of features removed from the samples that lost reads
min(rare_stats.df[which(rare_stats.df$Filtered_data_reads > rare_stats.df$Rarefied_data_reads),]$Features_removed_from_filtered_data_after_rarefaction)
max(rare_stats.df[which(rare_stats.df$Filtered_data_reads > rare_stats.df$Rarefied_data_reads),]$Features_removed_from_filtered_data_after_rarefaction)
median(rare_stats.df[which(rare_stats.df$Filtered_data_reads > rare_stats.df$Rarefied_data_reads),]$Features_removed_from_filtered_data_after_rarefaction)

# Proportion of reads removed from the samples that lost reads
min(rare_stats.df[which(rare_stats.df$Filtered_data_reads > rare_stats.df$Rarefied_data_reads),]$Proportion_reads_removed_from_filtered_data_after_rarefaction)
max(rare_stats.df[which(rare_stats.df$Filtered_data_reads > rare_stats.df$Rarefied_data_reads),]$Proportion_reads_removed_from_filtered_data_after_rarefaction)
median(rare_stats.df[which(rare_stats.df$Filtered_data_reads > rare_stats.df$Rarefied_data_reads),]$Proportion_reads_removed_from_filtered_data_after_rarefaction)

write.csv(rare_stats.df,"Result_tables/other/rarefaction_stats.csv",quote = F, row.names = F)

# (optional) only use rrarefied capped counts
otu.m <- otu_rare_count.m
temp <- read_counts_and_unique_features(otu.m, sort(c(sample_ids_original, "PNeg1P1")))
df2matrix(temp$sample_read_counts)
df2matrix(temp$sample_feature_counts)


# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# And re-calculate the abundances after filtering
otu_rel.m <- t(t(otu.m)/ colSums(otu.m))
otu_rel.m[is.nan(otu_rel.m)] <- 0
otu_rel_rare.m <- t(t(otu_rare_count.m) /colSums(otu_rare_count.m))
otu_rel_rare.m[is.nan(otu_rel_rare.m)] <- 0
otu_rel.m <- t(t(otu.m)/ colSums(otu.m))
otu_rel.m[is.nan(otu_rel.m)] <- 0

# Reassign sample IDs
sample_ids <- colnames(otu_rel.m)

## (Re)build dataframes

# Un-rarified
otu_rel.df <- data.frame("OTU.ID" = rownames(otu_rel.m))
otu_rel.df <- cbind(otu_rel.df, otu_rel.m[,colnames(otu_rel.m)])
rownames(otu_rel.df) <- c()

otu.df <- data.frame("OTU.ID" = rownames(otu.m))
otu.df <- cbind(otu.df, otu.m[,colnames(otu.m)])
rownames(otu.df) <- c()

# Write the final otu counts and abundances to file
write.table(otu.df, file = "Result_tables/count_tables/OTU_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(otu_rel.df, file = "Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)

# NOTE - otu.df, otu.m, otu_rel.df and otu_rel.m are the final, filtered OTU count / relative abundance dataframes and matrices. These can be used
# elsewhere for a variety of analyses at the OTU level, or, as is shown below, used to calculate abundances at different taxa levels


# Label those samples from the metadata.df that are not in the OTU table. A sample that has been filtered out
# will generally not be of interest for the study, though we may wish to have the metadata at hand
# Sample removed
samples_retained <- colnames(otu.df)[2:length(colnames(otu.df))]
samples_lost <- metadata.df$Index[!metadata.df$Index %in% samples_retained]
print(paste(length(samples_retained), "samples retained"))
print(paste(length(samples_lost), "samples lost"))

metadata.df$Sample_retained <- "no"
metadata.df[metadata.df$Index %in% samples_retained,]$Sample_retained <- "yes"

write.table(metadata.df[metadata.df$Sample_retained == "yes",], file = "Result_tables/other/processed_metadata.csv", sep = ",", quote = F, row.names = F)
write.table(metadata.df[metadata.df$Sample_retained == "no",], file = "Result_tables/other/metadata_samples_removed.csv", sep = ",", quote = F, row.names = F)


# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
#           Reads stats

# Rarefied counts                      = otu_rare_count.m 
# Rarefied relative abundances         = otu_rare_rel.m
# Relative abundances                  = otu_rel.m
# Counts                               = otu.m
# Relative abundances decontaminated  = otu_rel.m
# Counts decontaminated               = otu.m

# Unfiltered counts               = otu_unfiltered.m
# Unfiltered relative abundances  = otu_unfiltered_rel.m

# Only focus on those samples that have passed QC.
samples_passing_QC <- colnames(otu.m)
samples_in_unfiltered <- colnames(otu_unfiltered.m)
samples_in_unfiltered <- samples_in_unfiltered[!samples_in_unfiltered %in% negative_sample_ids] # Remove negative samples

# stats.df <- data.frame(Sample = samples_passing_QC)
stats.df <- data.frame(Sample = samples_in_unfiltered)
rownames(stats.df) <- stats.df$Sample

stats.df$Sample_retained <- "no"
stats.df[samples_passing_QC,]$Sample_retained <- "yes"
samples_not_retained <- rownames(stats.df[stats.df$Sample_retained == "no",])

stats.df$Community <- metadata_unfiltered.df[samples_in_unfiltered,"Community"]
stats.df$Gold_Star <- metadata_unfiltered.df[samples_in_unfiltered,"Gold_Star"]
stats.df$Season <- metadata_unfiltered.df[samples_in_unfiltered,"Season"]
stats.df$Nose <- metadata_unfiltered.df[samples_in_unfiltered,"Nose"]

# Add the raw read counts
stats.df[,"Raw_R1_read_counts"] <- metadata_unfiltered.df[rownames(stats.df),]$R1_read_count_raw

# Read counts prior to filtering out taxonomy / contaminants / low abundance features
stats.df[,"Original_read_counts"] <- colSums(otu_unfiltered.m[,samples_in_unfiltered]) 
stats.df[samples_passing_QC,"Filtered_read_counts"] <- colSums(otu.m[,samples_passing_QC])  # Read counts after filtering (contaminants, low read depth samples and low abundance features removed)
stats.df[samples_not_retained,"Filtered_read_counts"] <- colSums(otu_prior_to_removing_low_read_count_samples.m[,rownames(stats.df[samples_not_retained,])])
stats.df[samples_passing_QC,"Filtered_rarefied_read_counts"] <- colSums(otu_rare_count.m[,samples_passing_QC]) # "" and rarefied

# Reads removed from filtering
stats.df[,"Reads_removed_from_filtering"] <- stats.df[,"Original_read_counts"] - stats.df[,"Filtered_read_counts"]
stats.df[,"Proportion_reads_removed_from_filtering"] <- stats.df[,"Reads_removed_from_filtering"] / stats.df[,"Original_read_counts"]

# Reads removed from filtering and rarefaction
stats.df[,"Reads_removed_from_filtering_rarefied"] <- stats.df[,"Original_read_counts"] - stats.df[,"Filtered_rarefied_read_counts"]
stats.df[,"Proportion_reads_removed_from_filtering_rarefied"] <- stats.df[,"Reads_removed_from_filtering_rarefied"] / stats.df[,"Original_read_counts"]


# ---------------------------------------------
# Read counts and proportions original, unprocessed data. Summing each domain + unassigned should equal one

# Proportion of reads that are mammal (generally human)
stats.df[,"Mammal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("Mammal", project_otu_table_unfiltered.df$taxonomy_species),samples_in_unfiltered])
stats.df[,"Mammal_proportion_original"] <- stats.df[,"Mammal_read_count_original"] / stats.df[,"Original_read_counts"]

# Proportion of reads that are fungal
stats.df[,"Fungal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl(fungal_pattern, project_otu_table_unfiltered.df$taxonomy_species),samples_in_unfiltered])
stats.df[,"Fungal_proportion_original"] <- stats.df[,"Fungal_read_count_original"] / stats.df[,"Original_read_counts"]

# Proportion of reads that are bacterial
stats.df[,"Bacterial_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("d__Bacteria", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered])
stats.df[,"Bacterial_proportion_original"] <- stats.df[,"Bacterial_read_count_original"] / stats.df[,"Original_read_counts"]

# Proportion of reads that are Archaea
stats.df[,"Archaeal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("d__Archaea", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered])
stats.df[,"Archaeal_proportion_original"] <- stats.df[,"Archaeal_read_count_original"] / stats.df[,"Original_read_counts"]

# Proportion of reads that are Eukaryota
stats.df[,"Eukaryal_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("d__Eukaryota", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered])
stats.df[,"Eukaryal_proportion_original"] <- stats.df[,"Eukaryal_read_count_original"] / stats.df[,"Original_read_counts"]

# Proportion of reads that are Unassigned a taxonomy
stats.df[,"Unassigned_read_count_original"] <- colSums(project_otu_table_unfiltered.df[grepl("Unassigned", project_otu_table_unfiltered.df$Domain),samples_in_unfiltered])
stats.df[,"Unassigned_proportion_original"] <- stats.df[,"Unassigned_read_count_original"] / stats.df[,"Original_read_counts"]

# -------------------------------------
# Read counts and proportions filtered data
# otu_prior_to_removing_low_read_count_samples.m should be the same as otu.m except the low abundance samples have not been removed
temp <- otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% rownames(otu_prior_to_removing_low_read_count_samples.m),]
bacterial_otus_filtered <- temp[temp$Domain == "d__Bacteria",]$OTU.ID
fungal_otus_filtered <- temp[grepl("Fungi", temp$taxonomy_species),]$OTU.ID

# Proportion of reads in the filtered data that are bacterial
# temp <- stats.df
stats.df[,"Bacterial_read_count_after_filtering"] <- colSums(otu_prior_to_removing_low_read_count_samples.m[which(rownames(otu_prior_to_removing_low_read_count_samples.m) %in% bacterial_otus_filtered),])
stats.df[,"Bacterial_proportion_after_filtering"] <- stats.df[,"Bacterial_read_count_after_filtering"] / stats.df[,"Filtered_read_counts"]

# Proportion of reads in the filtered data that are Fungal
stats.df[,"Fungal_read_count_after_filtering"] <- colSums(otu_prior_to_removing_low_read_count_samples.m[which(rownames(otu_prior_to_removing_low_read_count_samples.m) %in% fungal_otus_filtered),])
stats.df[,"Fungal_proportion_after_filtering"] <- stats.df[,"Fungal_read_count_after_filtering"] / stats.df[,"Filtered_read_counts"]

# -------------------------------------
# Read counts and proportions filtered + rarefied data
temp <- otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID %in% rownames(otu_rare_count.m),]
bacterial_otus_filtered_rarefied <- temp[temp$Domain == "d__Bacteria",]$OTU.ID
# fungal_otus_filtered_rarefied <- temp[grepl("Fungi", temp$taxonomy_species),]$OTU.ID
fungal_otus_filtered_rarefied <- temp[grepl(fungal_pattern, temp$taxonomy_species),]$OTU.ID

# Proportion of reads in the filtered data that are bacterial
stats.df[samples_passing_QC,"Bacterial_read_count_after_filtering_rarefied"] <- colSums(otu_rare_count.m[which(rownames(otu_rare_count.m) %in% bacterial_otus_filtered_rarefied),samples_passing_QC])
stats.df[samples_passing_QC,"Bacterial_proportion_after_filtering_rarefied"] <- stats.df[samples_passing_QC,"Bacterial_read_count_after_filtering_rarefied"] / stats.df[samples_passing_QC,"Filtered_rarefied_read_counts"]

# Proportion of reads in the filtered data that are Fungal
stats.df[samples_passing_QC,"Fungal_read_count_after_filtering_rarefied"] <- colSums(otu_rare_count.m[which(rownames(otu_rare_count.m) %in% fungal_otus_filtered_rarefied),samples_passing_QC])
stats.df[samples_passing_QC,"Fungal_proportion_after_filtering_rarefied"] <- stats.df[samples_passing_QC,"Fungal_read_count_after_filtering_rarefied"] / stats.df[samples_passing_QC,"Filtered_rarefied_read_counts"]

## ------------------------------------------------------------------------------------
## This will calculate the total number of features across all samples and the number of non-zero features in each sample
# Can either calculate the feature numbers on just samples passing QC
# stats.df[,"Features_total"] <- length(which(rowSums(otu_unfiltered.m[,samples_passing_QC]) > 0 ))
# stats.df[,"Features_original"] <- apply(otu_unfiltered.m[,samples_passing_QC], 2, function(x) { length(which(x > 0)) } )

# Or calculate on all samples

stats.df[,"Features_total"] <- length(which(rowSums(otu_unfiltered.m[,samples_in_unfiltered]) > 0 ))
stats.df[,"Features_original"] <- apply(otu_unfiltered.m[,samples_in_unfiltered], 2, function(x) { length(which(x > 0)) } )
## ------------------------------------------------------------------------------------
stats.df[samples_passing_QC,"Features_filtered"] <- apply(otu.m[,samples_passing_QC], 2, function(x) { length(which(x > 0)) } )
stats.df[samples_not_retained,"Features_filtered"] <- apply(otu_prior_to_removing_low_read_count_samples.m[,samples_not_retained], 2, function(x) { length(which(x > 0)) } )
stats.df[samples_passing_QC,"Features_filtered_rarefied"] <- apply(otu_rare_count.m[,samples_passing_QC], 2, function(x) { length(which(x > 0)) } )

stats.df[,"Features_removed_from_filtering"] <- stats.df[,"Features_original"] - stats.df[,"Features_filtered"] 
stats.df[,"Features_removed_from_filtering_rarefied"] <- stats.df[,"Features_original"] - stats.df[,"Features_filtered_rarefied"] 

stats.df[,"Proportion_features_removed_from_filtering"] <- stats.df[,"Features_removed_from_filtering"] / stats.df[,"Features_original"]
stats.df[,"Proportion_features_removed_from_filtering_rarefied"] <- stats.df[,"Features_removed_from_filtering_rarefied"] / stats.df[,"Features_original"]


write.csv(stats.df, "Result_tables/other/QC_summary.csv", row.names = F, quote = F)



# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

# Above we processed the frequencies for each OTU to calculate the relative abundances.
# However, we often want the abundances at not just the individual OTU level, but also different taxonomy levels.
# For example, we may want to know the abundance of a particular Family.
# Now we will generate the abundance tables at each taxonomy level from Phylum, Class, Order, Family and Genus

# Merge the otu counts back with the taxonomy data
otu_metadata_merged.df <- merge(otu.df, otu_taxonomy_map.df, by.x = "OTU.ID", by.y = "OTU.ID")

tax_data.l <- generate_tax_level_data(otu_metadata_merged.df, 
                                      tax_string_levels = c("taxonomy_species", "taxonomy_genus", "taxonomy_family", "taxonomy_class", "taxonomy_order", "taxonomy_phylum"),
                                      sample_ids = sample_ids)

species.df <- m2df(tax_data.l$taxonomy_species$counts,"taxonomy_species")
species_rel.df <-  m2df(tax_data.l$taxonomy_species$abundances,"taxonomy_species")
genus.df <- m2df(tax_data.l$taxonomy_genus$counts,"taxonomy_genus")
genus_rel.df <-  m2df(tax_data.l$taxonomy_genus$abundances,"taxonomy_genus")
family.df <- m2df(tax_data.l$taxonomy_family$counts,"taxonomy_family")
family_rel.df <-  m2df(tax_data.l$taxonomy_family$abundances,"taxonomy_family")
order.df <- m2df(tax_data.l$taxonomy_order$counts,"taxonomy_order")
order_rel.df <-  m2df(tax_data.l$taxonomy_order$abundances,"taxonomy_order")
class.df <- m2df(tax_data.l$taxonomy_class$counts,"taxonomy_class")
class_rel.df <-  m2df(tax_data.l$taxonomy_class$abundances,"taxonomy_class")
phylum.df <- m2df(tax_data.l$taxonomy_phylum$counts,"taxonomy_phylum")
phylum_rel.df <-  m2df(tax_data.l$taxonomy_phylum$abundances,"taxonomy_phylum")



### Write the final counts and abundances for each taxonomy level to file
# Not rarefied
write.table(species.df, file = "Result_tables/count_tables/Specie_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(species_rel.df, file = "Result_tables/relative_abundance_tables/Specie_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus.df, file = "Result_tables/count_tables/Genus_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus_rel.df, file = "Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family.df, file = "Result_tables/count_tables/Family_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family_rel.df, file = "Result_tables/relative_abundance_tables/Family_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order.df, file = "Result_tables/count_tables/Order_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order_rel.df, file = "Result_tables/relative_abundance_tables/Order_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class.df, file = "Result_tables/count_tables/Class_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class_rel.df, file = "Result_tables/relative_abundance_tables/Class_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum.df, file = "Result_tables/count_tables/Phylum_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum_rel.df, file = "Result_tables/relative_abundance_tables/Phylum_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)

# ------------------------------------------------------------------------------------------------------------------------------
# Finally create and save a dataframe, separately for each Phylum, Class, Order, Family,Genus,Species and OTU ,
# containing the abundances/counts/log(counts, 10) in each sample, metadata and taxonomy information.


# Generate the combined datasets
otu_combined <- create_combined_dataframe_no_rare(counts.df = otu.df, 
                                                  abundances.df = otu_rel.df,
                                                  mylevel = "OTU.ID",
                                                  mymetadata = metadata.df,
                                                  otu_map.df = reduced_tax_map)

species_combined <- create_combined_dataframe_no_rare(counts.df = species.df, 
                                                      abundances.df = species_rel.df,
                                                      mylevel = "Species",
                                                      mymetadata = metadata.df,
                                                      otu_map.df = reduced_tax_map)

genus_combined <- create_combined_dataframe_no_rare(counts.df = genus.df, 
                                                    abundances.df = genus_rel.df,
                                                    mylevel = "Genus",
                                                    mymetadata = metadata.df,
                                                    otu_map.df = reduced_tax_map)

family_combined <- create_combined_dataframe_no_rare(counts.df = family.df, 
                                                     abundances.df = family_rel.df,
                                                     mylevel = "Family",
                                                     mymetadata = metadata.df,
                                                     otu_map.df = reduced_tax_map)

order_combined <- create_combined_dataframe_no_rare(counts.df = order.df, 
                                                    abundances.df = order_rel.df,
                                                    mylevel = "Order",
                                                    mymetadata = metadata.df,
                                                    otu_map.df = reduced_tax_map)

class_combined <- create_combined_dataframe_no_rare(counts.df = class.df, 
                                                    abundances.df = class_rel.df,
                                                    mylevel = "Class",
                                                    mymetadata = metadata.df,
                                                    otu_map.df = reduced_tax_map)

phylum_combined <- create_combined_dataframe_no_rare(counts.df = phylum.df, 
                                                     abundances.df = phylum_rel.df,
                                                     mylevel = "Phylum",
                                                     mymetadata = metadata.df,
                                                     otu_map.df = reduced_tax_map)

write.table(otu_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(species_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Specie_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Family_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Order_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Class_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Phylum_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)

