# Calculate the correlations in abundances between taxa

library(devtools)
library(Matrix)
library(igraph)
library(psych)
library(ggraph)
library(tidygraph)

otu_relabeller_function <- function(my_labels){
  taxonomy_strings <- unlist(lapply(my_labels, function(x) {
    as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_genus)
  }
  ))
  unlist(lapply(taxonomy_strings, function(x) {
    phylostring <- unlist(strsplit(x, split = ";"))
    paste(phylostring[3], phylostring[6], sep = ";")
  }))
}

genus_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  # paste(phylostring[2],phylostring[3], phylostring[6], sep = ";")
                  # paste(phylostring[3], phylostring[6], sep = ";")
                  paste(phylostring[3], phylostring[5], phylostring[6], sep = ";")
                }))
}

combined_otu_labeller <- function(x){
  # print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
  first_resolved_taxonomy(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
}

genus_relabeller_network <- function(x){
  # print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
  gsub(".*__(.*)","\\1", first_resolved_taxonomy(x))
}

source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T, row.names = "Sequence_file_ID_clean")
metadata.df <- metadata.df[!is.na(metadata.df$Otitis_Status),]

# Load feature taxonomy map
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)
rownames(otu_taxonomy_map.df) <- otu_taxonomy_map.df$OTU.ID

# Define descrete variables
discrete_variables <- c("Nose","Tympanic_membrane", "Otitis_Status",
                        "Season","Community","Gold_Star",
                        "H.influenzae_culture","M.catarrhalis_culture","S.pneumoniae_culture",
                        "Otitis_Status__Gold_Star", "Tympanic_membrane__Gold_Star",
                        "Community__Season","Community__Gold_Star","Community__Otitis_Status",
                        "H.Influenzae_qPCR", "M.catarrhalis_qPCR", "S.pneumoniae_qPCR",
                        "Corynebacterium_pseudodiphtheriticum","Dolosigranulum_pigrum","N_HRV")
discrete_variables_to_add_with_counts <- c("Community","Gold_Star","Season","Nose")


# Load count matrices
otu.df <- read.csv("Result_tables/count_tables/OTU_counts.csv", header =T)
genus.df <- read.csv("Result_tables/count_tables/Genus_counts.csv", header =T)
otu.m <- df2matrix(otu.df)
genus.m <- df2matrix(genus.df)



# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                           Generate fastspar inputs, only required for group subsets

prepare_input_variables <- c("Community", "Nose", "Otitis_Status", "Community__Gold_Star", "Community__Otitis_Status")

for (variable in prepare_input_variables){
  for (group in as.character(unique(metadata.df[,variable]))){
    sample_list <- as.character(subset(metadata.df, get(variable) == group)$Index)
    if (length(sample_list) < 4){
      print(paste0("Variable ", variable, ", group ", group, " has less than 2 samples, skipping"))
      next()
    }
    otu_subset_data.df <- otu.df[,c("OTU.ID", sample_list)]
    genus_subset_data.df <- genus.df[,c("taxonomy_genus", sample_list)]
    otu_subset_data.df <- otu_subset_data.df[otu_subset_data.df$OTU.ID %in% rownames(df2matrix(otu_subset_data.df)[which(apply(df2matrix(otu_subset_data.df), 1, sum) >= 50),]),]
    genus_subset_data.df <- genus_subset_data.df[genus_subset_data.df$taxonomy_genus %in% rownames(df2matrix(genus_subset_data.df)[which(apply(df2matrix(genus_subset_data.df), 1, sum) >= 50),]),]
    names(otu_subset_data.df)[1] <- "#OTU ID"
    names(genus_subset_data.df)[1] <- "#OTU ID"

    write.table(x = otu_subset_data.df, file = paste0("Result_tables/fastspar_inputs/", variable, "/", variable, "___",gsub("/|\\s", "_",group), "___otu_counts_fastspar.tsv"), sep = "\t", quote = F, row.names = F)
    write.table(x = genus_subset_data.df, file = paste0("Result_tables/fastspar_inputs/", variable, "/", variable, "___",gsub("/|\\s", "_",group), "___genus_counts_fastspar.tsv"), sep = "\t", quote = F, row.names = F)
  }
}

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                                             NETWORKS
# 
# NOTE - FASTSPAR SHOULD HAVE BEEN RUN SEPARATELY USING INPUTS GENERATED ABOVE
otu_cor_files <- list.files("Additional_results/fastspar/")[grepl("___otu___correlation.tsv", list.files("Additional_results/fastspar/"))]

source("code/helper_functions.R")
for (cor_file in otu_cor_files){

  otu_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  otu_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(cor_file)
  print(variable)
  print(dim(otu_fastspar_cor.m))
  otu_correlation_network.l <- generate_correlation_network(cor_matrix = otu_fastspar_cor.m,
                                                            p_matrix = otu_fastspar_pval.m,
                                                            relabeller_function = combined_otu_labeller,
                                                            p_value_threshold = 0.01,
                                                            cor_threshold = 0.6,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            node_label_segment_colour = "purple",
                                                            label_colour = "black",
                                                            label_size = 3,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            edge_alpha = 1,
                                                            # network_layout = "stress",
                                                            network_layout = "fr",
                                                            # network_layout = "kk",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            plot_title = paste0(variable, ": ", group, "; ASV correlation"),
                                                            filename= paste0("Result_figures/correlation_analysis/networks/otu/",variable,"___",group,"___feature_correlation_network.pdf"),
                                                            myseed = 1,
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = "pdf")
}

genus_cor_files <- list.files("Additional_results/fastspar/")[grepl("___genus___correlation.tsv", list.files("Additional_results/fastspar/"))]
genus_cor_files <- grep("Nose|Otitis", genus_cor_files,value =T)
file_type <- "svg"
for (cor_file in genus_cor_files){
  genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                               sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                                sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(variable)
  genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                              p_matrix = genus_fastspar_pval.m,
                                                              # relabeller_function = genus_relabeller_function,
                                                              relabeller_function = first_resolved_taxonomy,
                                                              # relabeller_function = genus_relabeller_network,
                                                              p_value_threshold = 0.05,
                                                              cor_threshold = 0.5,
                                                              node_size = 4,
                                                              node_colour = "grey20",
                                                              node_fill = "grey20",
                                                              node_label_segment_colour = "purple",
                                                              label_colour = "black",
                                                              label_size = 4,
                                                              plot_height = 10,
                                                              plot_width = 10,
                                                              edge_width_min = .5,
                                                              edge_width_max = 2.5,
                                                              edge_alpha = 1,
                                                              # network_layout = "stress",
                                                              network_layout = "fr",
                                                              # network_layout = "kk",
                                                              # exclude_to_from_df = edges_to_remove.df,
                                                              # plot_title = paste0(variable, ": ", group, "; Genus correlation"),
                                                              filename= paste0("Result_figures/correlation_analysis/networks/genus/",variable,"___",group,"___genus_correlation_network.",file_type),
                                                              myseed = 1,
                                                              edgetype = "link",
                                                              show_p_label = F,
                                                              file_type = file_type)
}

genus_fastspar_cor.m <- as.matrix(read.table("Additional_results/fastspar/Genus_correlation.tsv",
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table("Additional_results/fastspar/Genus_pvalues.tsv",
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
file_type <- "pdf"
genus_correlation_network.l <- generate_correlation_network(cor_matrix = genus_fastspar_cor.m,
                                                            p_matrix = genus_fastspar_pval.m,
                                                            relabeller_function = first_resolved_taxonomy,
                                                            # relabeller_function = genus_relabeller_network,
                                                            p_value_threshold = 0.05,
                                                            cor_threshold = 0.3,
                                                            node_size = 4,
                                                            node_colour = "grey20",
                                                            node_fill = "grey20",
                                                            node_label_segment_colour = "purple",
                                                            label_colour = "black",
                                                            label_size = 4,
                                                            plot_height = 10,
                                                            plot_width = 10,
                                                            edge_width_min = .5,
                                                            edge_width_max = 2.5,
                                                            edge_alpha = 1,
                                                            # network_layout = "stress",
                                                            network_layout = "fr",
                                                            # network_layout = "kk",
                                                            # exclude_to_from_df = edges_to_remove.df,
                                                            # plot_title = paste0(variable, ": ", group, "; Genus correlation"),
                                                            filename= paste0("Result_figures/correlation_analysis/networks/genus/Genus_correlation_network.",file_type),
                                                            myseed = 1,
                                                            edgetype = "link",
                                                            show_p_label = F,
                                                            file_type = file_type
                                                            )

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                           vs g__Dolosigranulum
#
# Break down Remote + Otitis status, 
# genus level, 
# abundances for "g__Moraxella", "g__Haemophilus", "g__Corynebacterium","g__Streptococcus" vs "g__Dolosigranulum"
# Include in insets plot_feature_correlations_external() results for g__Dolosigranulum

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

source("code/helper_functions.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

taxa_of_interest <- c("g__Dolosigranulum","g__Moraxella", "g__Haemophilus", "g__Corynebacterium","g__Streptococcus")

# ------------------------------------------------------------------------
# Loop over each result and save to file
genus_cor_files <- list.files("Additional_results/fastspar/")[grepl("___genus___correlation.tsv", list.files("Additional_results/fastspar/"))]
# genus_cor_files <- grep("Community__Otitis_Status.*genus", genus_cor_files,value =T)
genus_cor_files <- grep("^Otitis_Status.*genus", genus_cor_files,value =T)
# genus_cor_files <- grep("^Community___.*genus", genus_cor_files,value =T)

file_type = "svg"
for (cor_file in genus_cor_files){
  genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",cor_file),
                                               sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
  genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/",gsub("___correlation.tsv", "___pvalues.tsv",cor_file)),
                                                sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
  file_name_split <- strsplit(cor_file, split = "___")[[1]]
  variable <- file_name_split[1]
  group <- file_name_split[2]
  print(variable)
  print(group)
  
  rownames(genus_fastspar_cor.m) <- unlist(lapply(rownames(genus_fastspar_cor.m), first_resolved_taxonomy))
  colnames(genus_fastspar_cor.m) <- unlist(lapply(colnames(genus_fastspar_cor.m), first_resolved_taxonomy))
  rownames(genus_fastspar_pval.m) <- unlist(lapply(rownames(genus_fastspar_pval.m), first_resolved_taxonomy))
  colnames(genus_fastspar_pval.m) <- unlist(lapply(colnames(genus_fastspar_pval.m), first_resolved_taxonomy))
  
  print(dim(genus_fastspar_pval.m))
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                     feature = "g__Dolosigranulum",
                                     p_value_matrix = genus_fastspar_pval.m,
                                     top_n = 25,
                                     plot_width = 7, plot_height = 7,
                                     include_self = F,
                                     include_title = F,
                                     # filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                     filename= paste0("Result_figures/correlation_analysis/by_feature/",variable,"___",group,"___genus_dolosigranulum_correlations_top_25.",file_type),
                                     format = file_type)
  
  genus_fastspar_cor.m <- genus_fastspar_cor.m[taxa_of_interest,taxa_of_interest]
  print(genus_fastspar_cor.m)
  genus_fastspar_pval.m <- genus_fastspar_pval.m[taxa_of_interest,taxa_of_interest]
  
  plot_feature_correlations_external(cor_matrix = genus_fastspar_cor.m,
                                     feature = "g__Dolosigranulum",
                                     p_value_matrix = genus_fastspar_pval.m,
                                     top_n = 10,
                                     plot_width = 4, plot_height = 3,
                                     include_self = F,
                                     include_title = F,
                                     # filename = "Result_figures/correlation_analysis/by_feature/test.pdf",
                                     filename= paste0("Result_figures/correlation_analysis/by_feature/",variable,"___",group,"___genus_dolosigranulum_correlations.",file_type),
                                     format = file_type)
}


genus_fastspar_cor.m <- as.matrix(read.table(paste0("Additional_results/fastspar/Genus_correlation.tsv"),
                                             sep ="\t",header = T,row.names = 1,comment.char = "", check.names = F))
genus_fastspar_pval.m <- as.matrix(read.table(paste0("Additional_results/fastspar/Genus_pvalues.tsv"),
                                              sep ="\t",header = T,row.names = 1,comment.char = "",check.names = F))
rownames(genus_fastspar_cor.m) <- unlist(lapply(rownames(genus_fastspar_cor.m), first_resolved_taxonomy))
colnames(genus_fastspar_cor.m) <- unlist(lapply(colnames(genus_fastspar_cor.m), first_resolved_taxonomy))
rownames(genus_fastspar_pval.m) <- unlist(lapply(rownames(genus_fastspar_pval.m), first_resolved_taxonomy))
colnames(genus_fastspar_pval.m) <- unlist(lapply(colnames(genus_fastspar_pval.m), first_resolved_taxonomy))
