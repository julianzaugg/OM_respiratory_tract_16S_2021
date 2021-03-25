library(decontam)
library(microDecon)

library(DESeq2)

library(ComplexHeatmap)

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(gplots)
library(grid)
library(phyloseq)
library(seqinr) # For writing fasta files
library(corrplot)
library(svglite)
library(FSA)
library(ggsignif)




####################################
# Define various colour palettes
my_colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_palette_25_distinct <- c("#005e27","#34adff","#352138","#f10039","#92bb00","#99e0ff","#373300","#680015","#fca000","#d0dbb1","#507400","#9900b7","#3a32cd","#ff6f21","#ff9794","#e700a6","#54032d","#ffc74c","#dfda84","#8d2700","#271472","#eecbf4","#013e77","#d59bff","#ff81ce")
my_colour_palette_25_distinct_b <- c("#82ff7a","#f1beff","#26324c","#c82be3","#006f91","#7a003e","#efe5ff","#ff5bbc","#9bff42","#004e1c","#ce8c00","#3a246d","#c6b700","#7c1100","#2affbc","#92008e","#9a9900","#c60024","#ff753e","#009b8c","#ffd68a","#009008","#489aff","#6e19cc","#ecffc8")
my_colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
my_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
# my_colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#a85bd2","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
my_colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")

####################################


read_counts_and_unique_features <- function(count_matrix, sample_list){
  # Number of features and read counts
  temp <- count_matrix
  temp[temp > 0] <- 1
  read_sums <- colSums(count_matrix)
  feature_sums <- colSums(temp)
  sample_read_counts <- data.frame(sample_list, unlist(lapply(sample_list, 
                                                              function(x) ifelse(x %in% names(read_sums), read_sums[x],0))))
  sample_feature_counts <- data.frame(sample_list, unlist(lapply(sample_list, 
                                                                 function(x) ifelse(x %in% names(feature_sums), feature_sums[x],0))))
  out <- list("sample_read_counts" = sample_read_counts,
              "sample_feature_counts" = sample_feature_counts)
  out
}

generate_tax_level_data <- function(feature_count_taxonomy_data.df, sample_ids, tax_string_levels){
  
  df2matrix <- function(mydataframe){
    mymatrix <- mydataframe
    rownames(mymatrix) <- mydataframe[,1]
    mymatrix[,1] <- NULL
    mymatrix <- as.matrix(mymatrix)
    mymatrix
  }
  
  m2df <- function(mymatrix, column_name = "Row_variable"){
    mydf <- as.data.frame(mymatrix)
    cur_names <- names(mydf)
    mydf[, column_name] <- rownames(mydf)
    rownames(mydf) <- NULL
    mydf <- mydf[,c(column_name,cur_names)]
    return(mydf)
  }
  
  output <- list()
  for (tax_string_level in tax_string_levels){
    counts.df <- 
      feature_count_taxonomy_data.df[c(tax_string_level, sample_ids)] %>% 
      group_by_(tax_string_level) %>% # Group the dataframe by the taxonomy string 
      # Summarise each group by applying the the 'sum' function to the counts of each member of the group, 
      # i.e. duplicate entries for each taxa level are collapsed into a single entry and their counts summed together
      dplyr::summarise_all(funs(sum)) %>%  
      as.data.frame() # convert back to dataframe
    
    # Now create the relative abundance matrix at the current taxa level
    abundances.m <- counts.df
    
    rownames(abundances.m) <- abundances.m[[tax_string_level]]
    abundances.m[tax_string_level] <- NULL
    abundances.m <- as.matrix(abundances.m)
    abundances.m <- t(t(abundances.m) / colSums(abundances.m))
    abundances.m[is.nan(abundances.m)] <- 0
    
    output[[tax_string_level]] <- list("counts" = df2matrix(counts.df), "abundances" = abundances.m)
  }
  output
}


# Function
log_matrix <- function(mymat){
  out <- log(mymat, 10)
  out[is.infinite(out)] <- 0
  return(out)
}

# Function that calculates the geometric mean with some error-protection bits. 
# DESeq2 does not appear to work (will throw an error) if every OTU (or genus or genome etc.) 
# contains at least one count of zero in every row of the count data.
# Specifically, the function "dds<-DESeq(dds, betaPrior = FALSE)" will fail
# One way to address this is to use the function below as input to DESeq2 to transform the data.
# Calculate the geometric means prior to estimating the size factors
gm_mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# Center log ratio transform
clr <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

# For each rowname (OTU has), get the corresponding taxonomy_species
# Assumes "OTU.ID" and "taxonomy_species" columns in the provided map dataframe
assign_taxonomy_to_otu <- function(otutable, taxon_map){
  taxonomies <- c()
  for (otuid in rownames(otutable)){
    taxonomies <- c(taxonomies, as.character(taxon_map[taxon_map$OTU.ID == otuid,]$taxonomy_species))
  }
  return(taxonomies)
}

# Function that takes the metadata and a list of variables (column names) and returns those samples (rownames) with NA entries
get_samples_missing_data <- function(my_metadata, variables){
  samples_missing_data <- c()
  for (name in variables) {
    samples_missing_data <- c(samples_missing_data, rownames(my_metadata[is.na(my_metadata[[name]]),]))
  }
  return(unique(samples_missing_data))
}

bin_my_variable <- function(mydataframe, variable, breaks){
  temp <- cut(mydataframe[,variable], breaks = breaks,dig.lab = 4)
  temp <- gsub("\\(", "", temp)
  temp <- gsub("\\]", "", temp)
  temp <- gsub(",", "-", temp)
  return(temp)
}

first_resolved_taxonomy <- function(x) {
  # Get the lowest taxonomy level (or 'first') that has been resolved from a taxonomy string
  ranks.v <- unlist(strsplit(x, split = ';'))
  # print(x)
  # print(ranks.v)
  for (i in length(ranks.v):1) {
    split.v <- unlist(strsplit(ranks.v[i], split = '__'))
    if (!is.na(split.v[2]) & split.v[2] != "") {
      if (!grepl(split.v[2], pattern = "Unassigned|uncultured")) {
        return(ranks.v[i])
      }
    }
  }
  return(ranks.v[1])
}


filter_heatmap_matrix <- function(myheatmap, row_max = 0, prevalence = 0){
  internal_heatmap <- myheatmap
  internal_heatmap <- internal_heatmap[which(apply(internal_heatmap, 1, max) >= row_max), ]
  # keep only OTUs/taxa that are in more than this fraction of samples
  filter_fraction <- prevalence
  entry_prevalences <- apply(internal_heatmap, 1, function(x) {length(which(x > 0))})/dim(internal_heatmap)[2]
  entries_from_prevalences <- names(entry_prevalences)[entry_prevalences > filter_fraction]
  entries_from_prevalences <- entries_from_prevalences[!is.na(entries_from_prevalences)]
  return(internal_heatmap[entries_from_prevalences,])
}


df2matrix <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}

# Taxonomy-sample matrix to dataframe convertor
# Just converts a matrix where the row name is a taxonomy label (really can be anything)
m2df <- function(mymatrix, column_name = "Row_variable"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, column_name] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(column_name,cur_names)]
  return(mydf)
}

# Make row names of a dataframe the value of a defined column (default 1st) and then remove the column
clean_dataframe <- function(mydf, rowname_col = 1){
  my_clean.df <- mydf
  rownames(my_clean.df) <- my_clean.df[,rowname_col]
  my_clean.df[,1] <- NULL
  return(my_clean.df)
}


abundance_summary <- function(abundance.df, top_n = NULL){
  # Generate summary from abundance dataframe. Assumes first column is the unique taxa/feature ID
  abundance_melt.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
  abundance_summary.df <- 
    abundance_melt.df %>% 
    dplyr::mutate(Samples_total = n_distinct(Sample)) %>% # number of unique samples/index
    dplyr::group_by_at(1,.drop = F) %>% 
    dplyr::mutate(Samples_present = sum(Relative_abundance > 0, na.rm = TRUE)) %>%
    dplyr::summarise(
      Samples_present = max(Samples_present),
      Samples_total = max(Samples_total),
      # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
      Min_relative_abundance = round(min(Relative_abundance),5),
      Max_relative_abundance = round(max(Relative_abundance),5),
      Mean_relative_abundance = round(mean(Relative_abundance), 5),
      Median_relative_abundance = round(median(Relative_abundance), 5),
    ) %>%
    arrange(dplyr::desc(Mean_relative_abundance))
  if (!is.null(top_n)){
    abundance_summary.df <- 
      abundance_summary.df %>% 
      arrange(dplyr::desc(Mean_relative_abundance)) %>%
      head(top_n)
  }
  as.data.frame(abundance_summary.df)
}

collapse_abundance_dataframe <- function(abundance.df, filtering_taxa.v = NULL){
  # Collapse abundance dataframe. 
  # If filtering_taxa.v provided, other taxa are assigned 'Other' and values summed together
  abundance.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
  if (!is.null(filtering_taxa.v)){
    abundance.df[!abundance.df[,1] %in% filtering_taxa.v,][,1] <- "Other"
  }
  abundance.df <- 
    abundance.df %>% 
    group_by_at(c(2,1)) %>% 
    summarise_all(funs(sum)) %>%
    as.data.frame()
  abundance.df
}




# Calculate the (min, max, mean, median, stdev, #samples) abundances of each taxa at each taxa level
generate_taxa_summary <- function(mydata, taxa_column, group_by_columns = NULL){
  # if (is.null(group_by_columns)){
  #   select_columns <- c(taxa_column, "Read_count", "Relative_abundance")
  # } else{
  select_columns <- c(taxa_column, group_by_columns, "Sample", "Read_count", "Relative_abundance")
  total_samples <- length(unique(mydata$Sample))
  
  taxa_group_summary <- 
    mydata %>% 
    # dplyr::select_(.dots = select_columns) %>%
    dplyr::select(select_columns) %>%
    dplyr::group_by_(.dots = c(taxa_column, group_by_columns)) %>%
    dplyr::mutate(N_samples = n_distinct(Sample)) %>% # number of unique samples/index
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::mutate(N_total_samples_in_group = n_distinct(Sample))  %>%
    dplyr::group_by_(.dots = c(group_by_columns, taxa_column)) %>%
    dplyr::select(-Sample) %>%
    dplyr::summarise(N_samples = max(N_samples),
                     N_total_samples_in_group = max(N_total_samples_in_group),
                     Percent_group_samples = round((max(N_samples) / max(N_total_samples_in_group))*100, 2),
                     Percent_total_samples = round((max(N_samples) / total_samples)*100, 2),
                     
                     Mean_read_count = round(mean(Read_count), 2),
                     Median_read_count = median(Read_count),
                     Min_read_count = min(Read_count),
                     Max_read_count = max(Read_count),
                     Summed_read_count = sum(Read_count),
                     
                     Mean_relative_abundance_non_zero = round(mean(Relative_abundance), 5),
                     Mean_relative_abundance = round(sum(Relative_abundance)/N_total_samples_in_group, 5),
                     Median_relative_abundance = round(median(Relative_abundance), 5),
                     Min_relative_abundance = round(min(Relative_abundance),5),
                     Max_relative_abundance = round(max(Relative_abundance),5),
                     Summed_relative_abundance = round(sum(Relative_abundance),5)
    ) %>%
    as.data.frame()
  return(taxa_group_summary)
}


filter_summary_to_top_n <- function(taxa_summary, grouping_variables, abundance_column, my_top_n = 10){
  # Get the top N taxa as described in a provided taxa summary table.
  out <- 
    taxa_summary %>%
    dplyr::group_by_(.dots = c(grouping_variables)) %>%
    dplyr::arrange(dplyr::desc(get(abundance_column))) %>%
    dplyr::top_n(my_top_n, get(abundance_column)) %>% 
    # dplyr::arrange_(.dots = c(grouping_variables),abundance_column) %>%
    dplyr::arrange_(.dots = c(grouping_variables)) %>%
    as.data.frame()
  return(out)
}

# ---------------------------------------------------------------------------------------------------------
generate_pca_plot <- function(pca_object, # PCA object generated by the rda() function from vegan package
                              my_metadata.df, # Dataframe containing metadata
                              variable_to_plot, # Variable (column) that we are plotting
                              variable_colours_available = F, 
                              my_colour_palette = NULL,
                              my_levels = NULL,
                              
                              axis_limits = NULL, # list of axes limits c(xmin, xmax,ymin,ymax)
                              hide_grid = F,
                              show_x_label = T, # Show x-axis label
                              show_y_label = T, # Show y-axis label
                              plot_x_ticks = T, # Show x-axis tick marks
                              plot_y_ticks = T, # Show y-axis tick marks
                              plot_x_tick_labels = T, # Show x-axis tick labels
                              plot_y_tick_labels = T, # Show x-axis tick labels
                              plot_title = NULL, # Show plot title
                              title_cex = 1, # Size of title
                              
                              component_choices = c(1,2),
                              
                              plot_arrows = F, # Show contributing species/OTUs vectors
                              num_top_species = 5, # How many species/OTUs to show with the largest contributions to the variance
                              arrow_colour = "black", # Color of arrows
                              arrow_alpha = 1, # Transparency of arrows, 0 = invisible, 1= fully, visible, 0.5 = 50% transparent
                              arrow_thickness = .2, # Thickness of arrow
                              label_arrows=T, # Label the arrows
                              arrow_label_size = .5, # Size of arrow label
                              arrow_scalar = 1, # Scale of the arrows
                              arrow_label_colour = "black", # Colour of arrow label
                              arrow_label_font_type = 1, # Font type of arrow label : 1 Normal, 2 bold, 3 italic, 4 bold italic
                              arrow_label_offset = 0, # How much to offset the arrow label
                              arrow_label_alpha = 1, # The transparency of the arrow labels
                              
                              plot_hulls = F, # Plot hulls
                              hull_alpha = 1, # The transparency of the hulls
                              
                              plot_spiders = F, # Plot spiders
                              label_spider = F, # Label the spider centroid
                              spider_label_size = 0.5, # Size of spider label
                              spider_alpha = 1, # The transparency of the spiders
                              
                              plot_ellipses = F, # Plot ellipse
                              label_ellipse = F, # Label ellipse
                              ellipse_label_size = 0.5, # Size of ellipse label
                              ellipse_border_width = 1, # Width of ellipse border
                              ellipse_alpha = 1, # The transparency of the ellipses
                              
                              
                              label_sites = F, # Label the individual sites/samples
                              label_species = F, # Label the individual species/OTUs
                              
                              point_alpha = 1, # Transparency of points, 0 = invisible, 1= fully, visible, 0.5 = 50% transparent
                              point_size = 0.8, # Size of points
                              point_line_thickness = 1, # Thickness of point outline
                              use_shapes = T,
                              variable_shapes_available = F,
                              
                              include_legend = T, # Show the legend
                              legend_x = NULL, # x axis position of legend
                              legend_y = NULL, # y axis position of legend
                              legend_x_offset = 0, # How much to offset the legend on the x axis (assumes legend_x = NULL)
                              legend_y_offset = 0, # How much to offset the legend on the y axis (assumes legend_y = NULL)
                              legend_cex = 0.6, # Size of text in legend
                              legend_columns = 2, # Number of columns used for legend key
                              legend_key_text_distance = 0.5, # Distance between legend keys and text
                              legend_column_spacing = 0.5, # Distance between legend columns
                              legend_title = NULL, # Title of the legend. If NULL, just uses the variable
                              legend_fill_colour = NULL,
                              specie_labeller_function = NULL, # Function to re-label the species/OTUs,
                              
                              is_constrained = F, # This is a constrained ordination object
                              
                              file_type = "pdf",
                              filename = NULL,
                              plot_width = 10,
                              plot_height= 10
){
  # Function to generate a PCA plot
  
  pca.scores <- try(scores(pca_object, choices=c(1,2,3)))
  if(inherits(pca.scores, "try-error")) {
    return()
  }
  # Get component x,y coordinates
  pca_site_scores <- scores(pca_object, display = "sites", choices = component_choices)
  pca_specie_scores <- scores(pca_object, display = "species", choices = component_choices)
  
  # Check all entries in the PCA object are also in the metadata
  if (!all(rownames(pca_site_scores) %in% rownames(my_metadata.df))){
    error_message <- paste0("There are rows in the PCA object that are not defined in the metadata.
                             Rownames in the metadata need to match those used in the PCA object")
    stop(error_message)
  }
  # Remove NA entries from the metadata and from the PCA
  my_metadata.df <- my_metadata.df[!is.na(my_metadata.df[[variable_to_plot]]),]
  pca_site_scores <- pca_site_scores[rownames(pca_site_scores) %in% rownames(my_metadata.df),]
  
  if (is_constrained){
    pca_percentages <- (pca_object$CCA$eig/sum(pca_object$CCA$eig)) * 100  
  } else{
    pca_percentages <- (pca_object$CA$eig/sum(pca_object$CA$eig)) * 100  
  }
  
  # Assign a percentage of zero to NA values
  pca_percentages[is.na(pca_percentages)] <- 0
  if (length(pca_percentages) == 1){
    pca_percentages[2] <- 0
  }
  
  if (!is.null(axis_limits)){
    x_min <- axis_limits[1]
    x_max <- axis_limits[2]
    y_min <- axis_limits[3]
    y_max <- axis_limits[4]
  }else {
    x_min <- round(lapply(min(pca_site_scores[,1]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    x_max <- round(lapply(max(pca_site_scores[,1]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    y_min <- round(lapply(min(pca_site_scores[,2]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    y_max <- round(lapply(max(pca_site_scores[,2]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    
  }
  my_xlab <- ""
  my_ylab <- ""
  if (show_x_label){
    my_xlab = paste("PC",component_choices[1],"(", round(pca_percentages[component_choices[1]],1), "%)", sep = "")  
  }
  if (show_y_label){
    my_ylab = paste("PC", component_choices[2],"(", round(pca_percentages[component_choices[2]],1), "%)", sep = "")
  }
  
  my_metadata.df <- my_metadata.df[order(rownames(my_metadata.df)),]
  my_metadata.df <- my_metadata.df[order(my_metadata.df[[variable_to_plot]]),]
  
  # ------------------------------------------------------------------------------------
  # Ensure outcome variable is factored
  # Refactor the variable column so that the levels are consistent
  if (!is.null(my_levels)){
    my_metadata.df[,variable_to_plot] <- factor(my_metadata.df[,variable_to_plot], levels = my_levels)
  } else{
    # Uncomment to factorise and order levels alphabetically
    my_metadata.df[,variable_to_plot] <- factor(my_metadata.df[,variable_to_plot], 
                                                levels = sort(unique(as.character(my_metadata.df[,variable_to_plot]))))
    # Uncomment to factorise and inherit the current ordering
    # my_metadata.df[,variable_to_plot] <- factor(my_metadata.df[,variable_to_plot])
  }
  # ------------------------------------------------------------------------------------
  
  
  if (!is.null(filename)){
    if (file_type == "pdf"){
      pdf(file = filename,height = plot_height, width = plot_width)
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      svg(filename = filename,height = plot_height, width = plot_width)
      # svglite(file = filename,height = plot_height, width = plot_width)
    }
  }
  
  
  # dataframe.df <- dataframe.df[order(rownames(dataframe.df)),,drop = F]
  # dataframe.df <- dataframe.df[order(dataframe.df[[variable_to_plot]]),,drop = F]
  # dataframe.df[,variable_to_plot] <- factor(dataframe.df[,variable_to_plot], 
  #                                           levels = sort(unique(as.character(dataframe.df[,variable_to_plot]))))
  plot(0,
       type='n',
       # x = 0, y=0,
       xlim = c(x_min,x_max),
       ylim = c(y_min,y_max),
       xlab = my_xlab,
       ylab = my_ylab,
       xaxt = "n",
       yaxt = "n",
       # frame.plot = F,
       frame.plot = T
  )
  # Make grid
  if (hide_grid == F){
    grid(NULL,NULL, lty = 2, lwd = 1, col = "grey80")
  }
  
  # Add axes 
  axis(side = 1, labels = ifelse(plot_x_tick_labels, T, F), tck = -0.01,tick = ifelse(plot_x_ticks,T,F))
  axis(side = 2, labels = ifelse(plot_y_tick_labels, T, F), tck = -0.01, tick = ifelse(plot_x_ticks,T,F))
  
  # Assign (unique) colours and shapes for each grouping variable
  variable_values <- levels(my_metadata.df[[variable_to_plot]])
  
  # If variable colour column "variable_colour" in metadata, use colours from there
  if (variable_colours_available == T){
    colour_col_name <- paste0(variable_to_plot, "_colour")
    variable_colours <- setNames(as.character(unique(my_metadata.df[[colour_col_name]])), 
                                 as.character(unique(my_metadata.df[[variable_to_plot]])))
  } else{ # Otherwise use default palette
    colour_palette_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
    variable_colours <- setNames(rep(colour_palette_distinct,length(variable_values))[1:length(variable_values)], variable_values)
  }
  
  # Whether to use shapes
  if (use_shapes == T){
    # If variable shape column "variable_shape" in metadata, use shapes from there
    if (variable_shapes_available == T){
      shape_col_name <- paste0(variable_to_plot, "_shape")
      variable_shapes <- setNames(unique(my_metadata.df[[shape_col_name]]), 
                                  as.character(unique(my_metadata.df[[variable_to_plot]])))
      
    } else{ # Use default pre-defined shapes
      # variable_shapes <- setNames(c(25,24,23,22,21,8,6,5,4,3,2,1)[1:length(variable_values)],variable_values)
      variable_shapes <- setNames(rep(c(25,24,23,22,21),
                                      length(variable_values))[1:length(variable_values)],variable_values)      
    }
  } else{
    variable_shapes <- setNames(rep(c(21),length(variable_values))[1:length(variable_values)],variable_values)  
  }
  
  # variable_colours <- colours.l[[variable_to_plot]] # Assumes colours.l a global variable
  # variable_shapes <- pchs.l[[variable_to_plot]]  # Assumes pchs.l a global variable, also assumes shapes are only 25,24,23,22,21 (I think)
  # variable_shapes <- setNames(rep(c(25,24,23,22,21),length(variable_values))[1:length(variable_values)],variable_values)      
  
  # print(variable_values)
  # print(variable_colours)
  # print(variable_shapes)
  annotation_dataframe <- data.frame(variable_colours, variable_shapes, stringsAsFactors = F)
  annotation_dataframe$variable_outline_colours <- as.character(annotation_dataframe$variable_colours)
  if (point_line_thickness != 0){
    annotation_dataframe[annotation_dataframe$variable_shapes > 15,"variable_outline_colours"] <- "black"  
  }
  
  # Filter annotation dataframe to levels of interest
  if (!is.null(my_levels)){
    annotation_dataframe <- annotation_dataframe[my_levels,]
  }
  
  pca_site_scores <- pca_site_scores[rownames(my_metadata.df),]
  
  # Get the colours for all samples
  all_sample_colours <- as.character(
    lapply(
      as.character(my_metadata.df[rownames(pca_site_scores),variable_to_plot]), 
      function(x) as.character(annotation_dataframe[x,"variable_colours"])
    )
  )
  
  # Get the shapes for all samples
  all_sample_shapes <- as.numeric(
    lapply(
      as.character(my_metadata.df[rownames(pca_site_scores),variable_to_plot]), 
      function(x) annotation_dataframe[x,"variable_shapes"][[1]]
    )
  )
  
  # Set the outline colours for all samples based on the sample colours and refering to the annotation dataframe created above
  all_sample_outline_colours <- as.character(unlist(lapply(all_sample_colours, 
                                                           function(x) annotation_dataframe[annotation_dataframe$variable_colours == x, "variable_outline_colours"])))
  
  points(pca_site_scores, 
         cex = point_size,
         lwd = point_line_thickness,
         pch = all_sample_shapes,
         col = alpha(all_sample_outline_colours,point_alpha),
         # col = alpha("black",point_alpha),
         bg = alpha(all_sample_colours, point_alpha)
  )
  
  # Plot arrows for species / variables
  plot_arrows_func <- function(){
    
    left_pc1.v <- rownames(pca_specie_scores[order(pca_specie_scores[,1]),][1:num_top_species,])
    right_pc1.v <- rownames(pca_specie_scores[order(pca_specie_scores[,1]),][(length(pca_specie_scores[,1]) - num_top_species):length(pca_specie_scores[,1]),])
    
    left_pc2.v <- rownames(pca_specie_scores[order(pca_specie_scores[,2]),][1:num_top_species,])
    right_pc2.v <- rownames(pca_specie_scores[order(pca_specie_scores[,2]),][(length(pca_specie_scores[,2]) - num_top_species):length(pca_specie_scores[,2]),])
    
    top_vars.v <- unique(c(left_pc1.v, right_pc1.v, left_pc2.v, right_pc2.v))
    arrows(0,0, 
           arrow_scalar * pca_specie_scores[top_vars.v,1], 
           arrow_scalar * pca_specie_scores[top_vars.v,2], 
           length =0.05, 
           col = alpha(arrow_colour, arrow_alpha),
           lwd = arrow_thickness,
           lty = 1)
    
    if (label_arrows){
      if (!is.null(specie_labeller_function)){
        for (tv in top_vars.v){
          text(x = pca_specie_scores[tv,1]* arrow_scalar,
               y = pca_specie_scores[tv,2]* arrow_scalar,
               labels = specie_labeller_function(tv),
               cex = arrow_label_size,
               # Values of 1, 2, 3 and 4, respectively indicate positions below, 
               # to the left of, above and to the right of the specified (x,y) coordinates.
               pos = sample(c(1,2,3,4),1),
               # pos = 2,
               offset = arrow_label_offset,
               col = alpha(arrow_label_colour,alpha = arrow_label_alpha),
               font = arrow_label_font_type)
        }
      } else{
        for (tv in top_vars.v){
          text(x = pca_specie_scores[tv,1] * arrow_scalar,
               y = pca_specie_scores[tv,2] * arrow_scalar,
               labels = tv,
               cex = arrow_label_size,
               # Values of 1, 2, 3 and 4, respectively indicate positions below, 
               # to the left of, above and to the right of the specified (x,y) coordinates.
               pos = sample(c(1,2,3,4),1),
               # pos = 2,
               offset = arrow_label_offset,
               col = alpha(arrow_label_colour,alpha = arrow_label_alpha),
               font = arrow_label_font_type)
        }
      }
    }
  }
  
  if (plot_arrows == T){
    plot_arrows_func()
  }
  
  # Plot ellipses that are filled
  plot_ellipses_func <- function () {
    for (member in variable_values) {
      if (nrow(my_metadata.df[my_metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse
        ordiellipse(pca_site_scores,
                    groups = my_metadata.df[[variable_to_plot]],
                    kind = "ehull",
                    lwd = ellipse_border_width,
                    # border = variable_colours[member][[1]],
                    border = alpha(variable_colours[member][[1]],ellipse_alpha),
                    # col = variable_colours[member][[1]],
                    col = alpha(variable_colours[member][[1]],ellipse_alpha),
                    show.groups = member,
                    alpha = .05,
                    draw = "polygon",
                    label = F,
                    cex = .5)
      }
    }
  }
  
  # Plot hulls that are filled
  plot_hulls_func <- function () {
    for (member in variable_values){
      if (nrow(my_metadata.df[my_metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse}
        ordihull(pca_site_scores,
                 groups = my_metadata.df[[variable_to_plot]],
                 lwd = ellipse_border_width,
                 # border = variable_colours[member][[1]],
                 border = alpha(variable_colours[member][[1]],hull_alpha),
                 # col = variable_colours[member][[1]],
                 col = alpha(variable_colours[member][[1]],hull_alpha),
                 show.groups = member,
                 alpha = .05,
                 draw = "polygon",
                 label = F,
                 cex = .5)
      }
    }
  }
  if (hasArg(plot_hulls)){
    if (plot_hulls == T){
      plot_hulls_func()    
    }
  }
  
  #Plot spiders
  plot_spiders_func <- function (label_spider = F) {
    for (member in variable_values){
      if (nrow(my_metadata.df[my_metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse
        ordispider(pca_site_scores,
                   groups = my_metadata.df[[variable_to_plot]],
                   # col = variable_colours[member][[1]],
                   col = alpha(variable_colours[member][[1]],spider_alpha),
                   show.groups = member,
                   #alpha = .05,
                   label = label_spider,
                   cex = spider_label_size)
      }
    }
  }
  if (hasArg(plot_spiders)){
    if (plot_spiders == T){
      plot_spiders_func(label_spider = label_spider)    
    }
  }
  
  plot_ellipses_labels_func <- function(label_ellipse = F){
    # Repeat to have labels clearly on top of all ellipses
    for (member in variable_values){
      if (nrow(my_metadata.df[my_metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse
        ordiellipse(pca_site_scores,
                    groups = my_metadata.df[[variable_to_plot]],
                    kind = "ehull",
                    # border = variable_colours[member][[1]],
                    border = NA,
                    # col = variable_colours[member][[1]],
                    col = NA,
                    show.groups = member,
                    alpha = 0,
                    draw = "polygon",
                    label = label_ellipse,
                    cex = ellipse_label_size)
      }
    }
  }
  
  if (hasArg(plot_ellipses) | hasArg(label_ellipse)){
    if (plot_ellipses == T){
      plot_ellipses_func()    
      plot_ellipses_labels_func(label_ellipse = label_ellipse)
    } else if (label_ellipse == T){
      plot_ellipses_labels_func(label_ellipse = label_ellipse)
    }
  } 
  
  if (label_sites == T){
    text(x = pca_site_scores[,1],
         y = pca_site_scores[,2],
         labels = rownames(pca_site_scores),
         cex = .5,
         pos = 2)
  }
  if (label_species == T){
    text(x = pca_specie_scores[,1],
         y = pca_specie_scores[,2],
         labels = rownames(pca_specie_scores),
         cex = .5,
         pos = 2)
  }
  if (is.null(legend_x)){
    legend_x <- x_min + legend_x_offset
  }
  if (is.null(legend_y)){
    legend_y <- y_max + legend_y_offset
  }
  if (is.null(legend_title)){
    legend_title <- variable_to_plot
  }
  
  if (!is.null(plot_title)){
    title(main = plot_title, cex.main = title_cex)
  } 
  
  legend_bty = "n"
  if (!is.null(legend_fill_colour)){
    legend_bty <- "o"
    legend_fill <- legend_fill_colour
    legend_colour <- legend_fill_colour
  }
  
  if (include_legend){
    legend(
      # title = bold(variable_to_plot),
      title = as.expression(bquote(bold(.(legend_title)))),
      # title.adj = 0.5,
      title.col="black",
      # x = x_min-4,
      # y = y_max-6,
      x = legend_x,
      y = legend_y,
      # legend= variable_values,
      # pch= unique(all_sample_shapes),
      # col= legend_point_outline_colours,
      # pt.bg = unique(all_sample_colours),
      legend= rownames(annotation_dataframe),
      pch= annotation_dataframe$variable_shapes,
      col= as.character(annotation_dataframe$variable_outline_colours),
      pt.bg = as.character(annotation_dataframe$variable_colours),
      #bg = "white",
      bty = legend_bty,
      bg = legend_fill,
      box.col = legend_colour,
      ncol = legend_columns,
      cex = legend_cex,
      # pt.cex = 0.6,
      pt.lwd = point_line_thickness,
      y.intersp =1,
      x.intersp = legend_key_text_distance, # Distance between legend keys and text
      xjust = 0,
      # yjust = 1,
      title.adj = 0.5,
      text.width = legend_column_spacing # Distance between legend columns
    )
  }
  if (!is.null(filename)){
    dev.off()
  }
}



calculate_PC_taxa_contributions <- function(pca_object){
  # Calculate the percentage contribution from each taxa for PC1-3. Requires unscaled values that are squared
  pc1_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,1]^2, 3),value.name = "PC1_contribution_percentage")
  pc2_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC2_contribution_percentage")
  pc3_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC3_contribution_percentage")
  
  data.frame(pc1_contribution, pc2_contribution, pc3_contribution)
}

calculate_PC_abundance_correlations <- function(pca_object, mydata.df, taxa_column, variables = NULL){
  # pca_object <- genus_pca
  # mydata.df <- genus_data.df
  # taxa_column <-  "taxonomy_genus"
  # variables <- discrete_variables
  
  # Assume "Sample" and taxa_column are in mydata.df
  # mydata.df should be the combined table from the specific taxa level
  taxa_contributions <- calculate_PC_taxa_contributions(pca_object)
  
  # PC scores for each sample (site)
  pca_site_scores <- m2df(scores(pca_object, display = "sites"),"Sample") # These will be scaled scores!
  pca_species_scores <- m2df(scores(pca_object, display = "species"),taxa_column) # These will be scaled scores!
  
  # Filter to taxonomy in pca object (because we have filtered during the PCA calculations)
  abundance_pc_scores.df <- subset(mydata.df, get(taxa_column) %in% rownames(pca_object$CA$v))
  abundance_pc_scores.df <- left_join(abundance_pc_scores.df, pca_site_scores, by = "Sample")
  abundance_pc_scores.df <- abundance_pc_scores.df[,c("Sample","Relative_abundance", variables,taxa_column, "PC1", "PC2")]
  
  # Filter out entries where there are no PC scores
  abundance_pc_scores.df <- abundance_pc_scores.df[!is.na(abundance_pc_scores.df$PC1),]
  abundance_pc_scores.df$Relative_abundance <- abundance_pc_scores.df$Relative_abundance*100
  
  # Calculate the correlation between the abundances for each taxa and the PC1 and PC2 scores
  abundance_pc_correlations.df <- 
    abundance_pc_scores.df %>% 
    dplyr::group_by_(taxa_column) %>% dplyr::summarise(Pearson_PC1 = cor(PC1, Relative_abundance, method = "pearson"),
                                                       Pearson_PC2 = cor(PC2, Relative_abundance, method = "pearson"),
                                                       Spearman_PC1 = cor(PC1, Relative_abundance, method = "spearman"),
                                                       Spearman_PC2 = cor(PC2, Relative_abundance, method = "spearman"), 
                                                       N_Samples = n_distinct(Sample)) %>% 
    # filter(N_Samples >= 5) %>% 
    arrange(desc(abs(Pearson_PC1))) %>%
    as.data.frame()
  
  # Combine the correlation table with the PCA specie scores
  abundance_pc_correlations.df <- left_join(abundance_pc_correlations.df, pca_species_scores, by = taxa_column)
  
  
  # Add percentage contributions
  abundance_pc_correlations.df$PC1_contribution_percentage <- as.numeric(lapply(abundance_pc_correlations.df[,taxa_column],  
                                                                                function(x) taxa_contributions[x,]$PC1_contribution_percentage))
  abundance_pc_correlations.df$PC2_contribution_percentage <- as.numeric(lapply(abundance_pc_correlations.df[,taxa_column],  
                                                                                function(x) taxa_contributions[x,]$PC2_contribution_percentage))
  abundance_pc_correlations.df$PC3_contribution_percentage <- as.numeric(lapply(abundance_pc_correlations.df[,taxa_column], 
                                                                                function(x) taxa_contributions[x,]$PC3_contribution_percentage))
  abundance_pc_correlations.df
}
# ---------------------------------------------------------------------------------------------------------



# Function to create heatmap
make_heatmap <- function(myheatmap_matrix,
                         mymetadata = NULL,
                         filename= NULL,
                         #...,
                         # Dataframe with two columns. First must match row entry, second the new label
                         my_row_labels = NULL, 
                         my_col_labels = NULL, # Same as my_row_labels, though with columns
                         height = 10, # Not currently used
                         width = 10, # Not currently used
                         heatmap_height = 10, # Not currently used
                         heatmap_width = 10, # Not currently used
                         plot_height =10,
                         plot_width =10,
                         column_title_size = 10,
                         row_title_size = 10,
                         annotation_bar_name_size = 10,
                         variables = NULL, # Annotations
                         cluster_columns = T,
                         cluster_rows = T,
                         my_breaks = NULL,
                         my_palette = NULL,
                         palette_choice = NULL,
                         legend_title = NULL,
                         column_title = "Sample",
                         row_title = "Taxa",
                         legend_labels = NULL,
                         my_annotation_palette = NULL,
                         discrete_legend = FALSE, # Whether or not to display continuous legend as discrete
                         simple_anno_size = unit(.5, "cm"), # size of annotations
                         show_column_dend = F,
                         show_row_dend = F,
                         do_not_order = F,
                         show_cell_values = F,
                         # If show_cell_values =T, cells less than this will have a black font colour
                         # and above white
                         cell_fun_value_col_threshold = 15,
                         my_cell_fun = NULL,
                         show_legend = T,
                         show_top_annotation = T,
                         row_name_size = 6,
                         col_name_size = 6,
                         col_annotation_label_size = 6,
                         col_annotation_title_size = 6,
                         col_annotation_legend_grid_height = .2,
                         col_annotation_legend_grid_width = .2,
                         scale_legend_label_size = 6,
                         scale_legend_title_size = 6,
                         grid_thickness = 1,
                         grid_colour = "white",
                         ...){
  # print(list(...))
  argList<-list(...) # argument list for checking unspecified optional parameters
  # print(argList$cell_fun)
  # return(1)
  
  # Assign internal objects
  internal_heatmap_matrix.m <- myheatmap_matrix
  internal_metadata.df <- mymetadata
  
  
  if (!is.null(internal_metadata.df)){
    # Order/filter the heatmap matrix to order/entries of metadata
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,rownames(internal_metadata.df),drop = F]
    if (!is.null(variables)){
      # Order the heatmap matrix by the variables
      internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,do.call(order, internal_metadata.df[,variables,drop=F]),drop =F]
      # Order the metadata by the variables
      internal_metadata.df <- internal_metadata.df[do.call(order, internal_metadata.df[,variables,drop=F]),,drop=F]
      # Create metadata just containing the variables
      metadata_just_variables <- internal_metadata.df[,variables, drop = F]    
    }
    # Check that rownames match colnames
    if (!all(rownames(internal_metadata.df) == colnames(internal_heatmap_matrix.m))){
      stop("Row names in metadata do not match column names in matrix")
    }
    
    # Create annotations
    colour_lists <- list()
    for (myvar in variables){
      var_colour_name <- paste0(myvar, "_colour")
      # Assumes there is a colour column for each variable in the metadata
      # If there is no colour column, create one and assign from palette
      # internal_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
      # internal_colour_palette_10_distinct <- my_colour_palette_20_distinct
      if (is.null(my_annotation_palette)){
        internal_colour_palette <- my_colour_palette_206_distinct
      } else{
        internal_colour_palette <- my_annotation_palette
      }
      if (!var_colour_name %in% names(internal_metadata.df)){
        myvar_values <- factor(as.character(sort(unique(internal_metadata.df[,myvar]))))
        myvar_colours <- setNames(internal_colour_palette[1:length(myvar_values)], myvar_values)
        all_variable_colours <- as.character(lapply(as.character(internal_metadata.df[,myvar]), function(x) myvar_colours[x]))
        internal_metadata.df[,paste0(myvar,"_colour")] <- all_variable_colours
      }
      
      metadata_subset <- unique(internal_metadata.df[,c(myvar, var_colour_name)])
      # Order by the variable column
      metadata_subset <- metadata_subset[order(metadata_subset[,myvar]),]
      # Factorise the variable column
      metadata_subset[,myvar] <- factor(metadata_subset[,myvar])
      metadata_subset <- metadata_subset[!is.na(metadata_subset[,myvar]),]
      named_colour_list <- setNames(as.character(metadata_subset[, var_colour_name]), as.character(metadata_subset[,myvar]))
      colour_lists[[myvar]] <- named_colour_list
    }
    
    # Appearance of the column annotations
    #HeatmapAnnotation
    if (show_top_annotation == T){
      if(is.null(variables)){
        print("No variables specified, cannot add annotations")
        show_top_annotation <- F
      } else{
        ha <- columnAnnotation(df = metadata_just_variables,
                               # which = "column",
                               col = colour_lists,
                               gp = gpar(col = "black",lwd =.2),
                               gap = unit(.1,"cm"),
                               show_annotation_name = T,
                               # annotation_legend_param, # ?color_mapping_legend for options
                               show_legend = show_legend,
                               simple_anno_size = simple_anno_size,
                               annotation_legend_param = list(labels_gp = gpar(fontsize = col_annotation_label_size),
                                                              title_gp = gpar(fontsize = col_annotation_title_size,fontface = "bold"),
                                                              grid_height = unit(col_annotation_legend_grid_height,"cm"),
                                                              grid_width = unit(col_annotation_legend_grid_width,"cm")),
                               annotation_name_gp = gpar(fontsize = annotation_bar_name_size))      
      }
      # HeatmapAnnotation(annotation_)
      # ?color_mapping_legend
      
    }
  }
  
  # TODO - add option for row annotation
  if (is.null(my_palette)){
    if (is.null(palette_choice)) {palette_choice <- "blue"}
    if (!palette_choice %in% c("blue", "purple","red","dark_bluered","bluered")) { palette_choice <- "blue"}
    if (palette_choice == "blue"){
      my_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))
    } 
    else if (palette_choice == "purple"){
      my_palette <- colorRampPalette(c("white", "#f9cdac","#f3aca2", "#ee8b97", "#e96a8d","#db5087","#b8428c", "#973490", "#742796","#5e1f88", "#4d1a70", "#3d1459","#2d0f41"))
    } else if (palette_choice == "red"){
      my_palette <- colorRampPalette(c("white", "#fded86","#fde86e", "#f9d063", "#f5b857","#f0a04b","#eb8a40", "#e77235","#e35b2c", "#c74e29","#9d4429","#753c2c","#4c3430"))
    } else if (palette_choice == "dark_bluered"){
      my_palette <- colorRampPalette(c("#08306B","#FFD92F","#67001F"))
    } else if (palette_choice == "bluered"){
      my_palette <- colorRampPalette(c("#17468a","#ffdd47","#99113a"))
    }
  } else{
    my_palette <- colorRampPalette(my_palette)
  }
  
  if (!is.null(my_breaks)){
    internal_breaks <- my_breaks
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = my_palette(length(internal_breaks)))
  } else{
    internal_breaks <- seq(min(internal_heatmap_matrix.m), max(internal_heatmap_matrix.m), length.out = 6)
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = my_palette(length(internal_breaks)))
  }
  
  my_row_labels.v = rownames(internal_heatmap_matrix.m)
  if (!is.null(my_row_labels)){
    my_row_labels.v <- as.character(lapply(my_row_labels.v, function(x) as.character(my_row_labels[my_row_labels[,1] == x,][,2])))
  }
  my_col_labels.v = colnames(internal_heatmap_matrix.m)
  if (!is.null(my_col_labels)){
    my_col_labels.v <- as.character(lapply(my_col_labels.v, function(x) as.character(my_col_labels[my_col_labels[,1] == x,][,2])))
  }
  
  if (do_not_order != T){
    # Order the heatmap rows by the row labels names
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[order(my_row_labels.v),]
    my_row_labels.v <- my_row_labels.v[order(my_row_labels.v)]    
  }
  
  # if show values and no function provided
  if (show_cell_values == T & is.null(my_cell_fun)){ 
    my_cell_fun <- function(j, i, x, y, width, height, fill) {
      # if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold & internal_heatmap_matrix.m[i, j] != 0){
      if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold){
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "black"))}
      else if(internal_heatmap_matrix.m[i, j] >= cell_fun_value_col_threshold ) {
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "white"))
      }
    }
  }
  
  # Legend appearance
  if (is.null(legend_labels)){
    my_labels <- internal_breaks
  } else{
    my_labels <- legend_labels
  }
  if (discrete_legend == TRUE){
    hm_legend <- Legend(
      labels = rev(my_labels),
      at = internal_breaks,
      labels_gp = gpar(fontsize = scale_legend_label_size),
      legend_gp = gpar(fill = rev(col_fun(internal_breaks))), # For discrete
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = scale_legend_title_size,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  } else{
    hm_legend <- Legend(
      col_fun = col_fun, # For continuous
      labels = my_labels,
      at = internal_breaks,
      labels_gp = gpar(fontsize = scale_legend_label_size),
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = scale_legend_title_size,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  }
  
  hm <- Heatmap(matrix = internal_heatmap_matrix.m,
                
                # top_annotation = ha,
                
                # Colours
                col = col_fun,
                na_col = "grey",
                
                # Sizing
                show_heatmap_legend = F,
                heatmap_legend_param = list(hm_legend),
                row_names_max_width = unit(35,"cm"),
                row_labels = my_row_labels.v,
                column_labels = my_col_labels.v,
                # row_names_side = "left",
                # height = unit(height,"cm"),
                # width = unit(width,"cm"),
                # heatmap_height = unit(heatmap_height,"cm"),
                # heatmap_width = unit(heatmap_width,"cm"),
                # heatmap_width = unit(15,"cm"),
                
                # Titles
                column_title = column_title,
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize = column_title_size),
                row_title = row_title,
                row_title_side = "left",
                row_title_gp = gpar(fontsize = row_title_size),
                
                # Clustering
                cluster_columns = cluster_columns,
                cluster_rows = cluster_rows,
                clustering_method_columns = "average",
                clustering_method_rows = "average",
                show_column_dend = show_column_dend, 
                show_row_dend = show_row_dend,
                # column_dend_height = unit(2, "cm"),
                # row_dend_width = unit(3, "cm"),
                
                # Borders
                border = F,
                rect_gp = gpar(col = grid_colour, lwd = grid_thickness),
                
                # Text appearance
                row_names_gp = gpar(fontsize = row_name_size),
                column_names_gp = gpar(fontsize = col_name_size),
                cell_fun = my_cell_fun,
                ...
  )
  
  if (show_top_annotation == T & exists("ha")){
    hm <- ha %v% hm
  }
  if (!is.null(filename)){
    pdf(filename,height=plot_height,width=plot_width)
    draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T)
    dev.off()    
  }
  draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T)
  return(list("heatmap" = hm, "legend" = hm_legend))
  
}


# ******************************************************************************************
# Taken from : https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
# ******************************************************************************************

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
calculate_alpha_diversity_significance <- function(mydata, variable){
  # Assumes there are Shannon, Chao1 and Simpson columns
  # This is run assuming unpaired data.
  results.df <- data.frame("Variable" = character(),
                           "Group_1" = character(),
                           "Group_2" = character(),
                           "N_samples_group_1" = character(),
                           "N_samples_group_2" = character(),
                           "Shannon_MannW_pvalue" = character(),
                           "Simpson_MannW_pvalue" = character(),
                           "Chao1_MannW_pvalue" = character(),
                           # "Shannon_KrusW_pvalue" = character(),
                           # "Simpson_KrusW_pvalue" = character(),
                           # "Chao1_KrusW_pvalue" = character(),
                           "Shannon_MannW_padj" = character(),
                           "Simpson_MannW_padj" = character(),
                           "Chao1_MannW_padj" = character()
                           # "Shannon_KrusW_padj" = character(),
                           # "Simpson_KrusW_padj" = character(),
                           # "Chao1_KrusW_padj" = character()
                           )
  
  if (length(as.character(unique(mydata[,variable]))) < 2){
    return(results.df)
  }
  group_combinations <- combn(as.character(unique(mydata[,variable])), 2)
  
  for (i in 1:ncol(group_combinations)) {
    group_1 <- group_combinations[1,i]
    group_2 <- group_combinations[2,i]
    group_1_meta <- subset(mydata, get(variable) == group_1)
    group_2_meta <- subset(mydata, get(variable) == group_2)
    if (is.na(group_1) | is.na(group_2)) {next}
    if (dim(group_1_meta)[1] == 0 | dim(group_2_meta)[1] == 0) {next}
    N_samples_group_1 <- dim(group_1_meta)[1]
    N_samples_group_2 <- dim(group_2_meta)[1]
    
    # Mann-Whitney test on the Shannon diversity
    wilcox_shannon_test <- wilcox.test(group_1_meta$Shannon, group_2_meta$Shannon, exact = F)
    # Mann-Whitney test on the Simpson diversity
    wilcox_simpson_test <- wilcox.test(group_1_meta$Simpson, group_2_meta$Simpson, exact = F)
    # Mann-Whitney test on the Chao1 diversity
    wilcox_chao1_test <- wilcox.test(group_1_meta$Chao1, group_2_meta$Chao1, exact = F)
    
    # Kruskal-Wallis (pairwise) test on the Shannon diversity
    kruskal_shannon_test <- kruskal.test(Shannon~get(variable), data = subset(mydata, get(variable) %in% c(group_1, group_2)))
    # Kruskal-Wallis (pairwise) test on the Simpson diversity
    kruskal_simpson_test <- kruskal.test(Simpson~get(variable), data = subset(mydata, get(variable) %in% c(group_1, group_2)))
    # Kruskal-Wallis (pairwise) test on the Chao1 diversity
    kruskal_chao1_test <- kruskal.test(Chao1~get(variable), data = subset(mydata, get(variable) %in% c(group_1, group_2)))
    
    results.df <- rbind(results.df, data.frame("Variable" = variable,
                                               "Group_1" = group_1, 
                                               "Group_2" = group_2, 
                                               "N_samples_group_1" = N_samples_group_1,
                                               "N_samples_group_2" = N_samples_group_2,
                                               "Shannon_MannW_pvalue" = round(wilcox_shannon_test$p.value,6),
                                               "Simpson_MannW_pvalue" = round(wilcox_simpson_test$p.value,6),
                                               "Chao1_MannW_pvalue" = round(wilcox_chao1_test$p.value,6),
                                               # "Shannon_KrusW_pvalue" = round(kruskal_shannon_test$p.value,6),
                                               # "Simpson_KrusW_pvalue" = round(kruskal_simpson_test$p.value,6),
                                               # "Chao1_KrusW_pvalue" = round(kruskal_chao1_test$p.value,6),
                                               
                                               # "Shannon_MannW_padj" = round(wilcox_shannon_test$p.value,6),
                                               # "Simpson_MannW_padj" = round(wilcox_simpson_test$p.value,6),
                                               # "Chao1_MannW_padj" = round(wilcox_chao1_test$p.value,6),
                                               # "Shannon_KrusW_padj" = round(kruskal_shannon_test$p.value,6),
                                               # "Simpson_KrusW_padj" = round(kruskal_simpson_test$p.value,6),
                                               # "Chao1_KrusW_padj" = round(kruskal_chao1_test$p.value,6)),
                                               "Shannon_MannW_padj" = NA,
                                               "Simpson_MannW_padj" = NA,
                                               "Chao1_MannW_padj" = NA
                                               # "Shannon_KrusW_padj" = NA,
                                               # "Simpson_KrusW_padj" = NA,
                                               # "Chao1_KrusW_padj" = NA)
                        ))
                        # "Shannon_MannW_padj" = "x",
                        # "Simpson_MannW_padj" = "x",
                        # "Chao1_MannW_padj" = "x",
                        # "Shannon_KrusW_padj" = "x",
                        # "Simpson_KrusW_padj" = "x",
                        # "Chao1_KrusW_padj" = "x"))
  }

  results.df$Shannon_MannW_padj <- round(p.adjust(results.df$Shannon_MannW_pvalue,method = "BH"),6)
  results.df$Simpson_MannW_padj <- round(p.adjust(results.df$Simpson_MannW_pvalue,method = "BH"),6)
  results.df$Chao1_MannW_padj <- round(p.adjust(results.df$Chao1_MannW_pvalue,method = "BH"),6)
  # results.df$Shannon_KrusW_padj <- round(p.adjust(results.df$Shannon_KrusW_pvalue,method = "BH"),6)
  # results.df$Simpson_KrusW_padj <- round(p.adjust(results.df$Simpson_KrusW_pvalue,method = "BH"),6)
  # results.df$Chao1_KrusW_padj <- round(p.adjust(results.df$Chao1_KrusW_pvalue,method = "BH"),6)
  results.df
}


calculate_alpha_diversity_significance_multiple <- function(mydata, variable){
  # Assumes there are Shannon, Chao1 and Simpson columns
  mydata <- mydata[!is.na(mydata[,variable]),]
  n_groups = length(as.character(unique(mydata[,variable])))
  if (any(is.na(mydata[,variable]))){
    return()
  }
  if (n_groups > 2){
    
    kw_shannon <- kruskal.test(Shannon~get(variable), data = mydata)
    kw_simpson <- kruskal.test(Simpson~get(variable), data = mydata)
    kw_chao1 <- kruskal.test(Shannon~get(variable), data = mydata)
    dunn_shannon <- dunnTest(x = Shannon~get(variable), data = mydata, method = "bh", alpha = 0.05)
    dunn_simpson <- dunnTest(x = Simpson~get(variable), data = mydata, method = "bh", alpha = 0.05)
    dunn_chao1 <- dunnTest(x = Chao1~get(variable), data = mydata, method = "bh", alpha = 0.05)
    
    dunn_shannon <- separate(dunn_shannon$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
    names(dunn_shannon) <- c("Group_1","Group_2","Shannon_Dunn_pvalue","Shannon_Dunn_padj")
    
    dunn_simpson <- separate(dunn_simpson$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
    names(dunn_simpson) <- c("Group_1","Group_2","Simpson_Dunn_pvalue","Simpson_Dunn_padj")
    
    dunn_chao1 <- separate(dunn_chao1$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
    names(dunn_chao1) <- c("Group_1","Group_2","Chao1_Dunn_pvalue","Chao1_Dunn_padj")
    
    multiple_group_comparison.df <- merge(merge(x = dunn_shannon, y = dunn_simpson, by = c("Group_1", "Group_2")), y = dunn_chao1,  by = c("Group_1", "Group_2"))
    multiple_group_comparison.df$Shannon_KrusW_pvalue <- kw_shannon$p.value
    multiple_group_comparison.df$Simpson_KrusW_pvalue <- kw_simpson$p.value
    multiple_group_comparison.df$Chao1_KrusW_pvalue <- kw_chao1$p.value
    multiple_group_comparison.df$Variable = variable
    
    multiple_group_comparison.df$N_samples_group_1 <- NA
    multiple_group_comparison.df$N_samples_group_2 <- NA
    # print(multiple_group_comparison.df)
    for (row in 1:nrow(multiple_group_comparison.df)){
      group_1 <- multiple_group_comparison.df[row,"Group_1"]
      group_2 <- multiple_group_comparison.df[row,"Group_2"]
      group_1_meta <- subset(mydata, get(variable) == group_1)
      group_2_meta <- subset(mydata, get(variable) == group_2)
      N_samples_group_1 <- dim(group_1_meta)[1]
      N_samples_group_2 <- dim(group_2_meta)[1]
      multiple_group_comparison.df[row,"N_samples_group_1"] <- N_samples_group_1
      multiple_group_comparison.df[row,"N_samples_group_2"] <- N_samples_group_2
    }
    
    return(multiple_group_comparison.df)
  }
  multiple_group_comparison.df <- data.frame("Variable" = character(),
                                             "Group_1" = character(),
                                             "Group_2" = character(),
                                             "N_samples_group_1" = character(),
                                             "N_samples_group_2" = character(),
                                             "Shannon_Dunn_pvalue" = character(),
                                             "Shannon_Dunn_padj" = character(),
                                             "Simpson_Dunn_pvalue" = character(),
                                             "Simpson_Dunn_padj" = character(),
                                             "Chao1_Dunn_pvalue" = character(),
                                             "Chao1_Dunn_padj" = character(),
                                             "Shannon_KrusW_pvalue" = character(),
                                             "Simpson_KrusW_pvalue" = character(),
                                             "Chao1_KrusW_pvalue" = character())
  return(multiple_group_comparison.df)
  
}


summarise_alpha_diversities <- function(mydata, group_by_columns){
  summary.df <- mydata %>% 
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::summarise(Shannon_Mean =mean(Shannon),
                     Shannon_Stdev=sd(Shannon),
                     Shannon_Max=max(Shannon), 
                     Shannon_Min=min(Shannon), 
                     Shannon_Median=median(Shannon), 
                     
                     Simpson_Mean=mean(Simpson), 
                     Simpson_Stdev=sd(Simpson),
                     Simpson_Max=max(Simpson), 
                     Simpson_Min=min(Simpson), 
                     Simpson_Median=median(Simpson), 
                     
                     Chao1_Mean=mean(Chao1), 
                     Chao1_Stdev=sd(Chao1),
                     Chao1_Max=max(Chao1), 
                     Chao1_Min=min(Chao1), 
                     Chao1_Median=median(Chao1),
                     
                     N_samples = n_distinct(Index)
    ) %>% as.data.frame()
  summary.df
}


summarise_diversities_each_variable <- function(mydata, variables){
  mydata_summary.df <- NULL
  for (myvar in variables){
    mydata_subset.df <- mydata[which(!is.na(mydata[,myvar])),,drop =F]
    if (is.null(mydata_summary.df)){
      mydata_summary.df <- summarise_alpha_diversities(mydata_subset.df, myvar)
      # Melt because we are combining variables into a single table
      mydata_summary.df <- melt(mydata_summary.df, measure.vars = myvar, 
                                value.name = "Group", 
                                variable.name = "Variable")
    } else{
      mydata_summary.df <- rbind(mydata_summary.df, melt(summarise_alpha_diversities(mydata_subset.df, myvar),
                                                         measure.vars = myvar, 
                                                         value.name = "Group", 
                                                         variable.name = "Variable"))      
    }
    
  }
  mydata_summary.df
}

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


run_permanova_custom <- function(my_metadata, my_formula, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  result <- adonis(my_formula,data = my_metadata, permu=permutations,method= my_method)
  # result <- adonis(my_formula,data = my_metadata, permu=999,method="bray")
  for (r in rownames(result$aov.tab)){
    variable <- r
    Degress_of_freedom <- result$aov.tab[r,]$Df[1]
    SumOfSqs <- round(result$aov.tab[r,]$SumsOfSqs[1], 3)
    meanSqs <- round(result$aov.tab[r,]$MeanSqs[1], 3)
    F.model <- round(result$aov.tab[r,]$F.Model[1], 3)
    R2 <- round(result$aov.tab[r,]$R2[1], 3)
    p_value <- round(result$aov.tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degress_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       R2,
                                                       p_value))
  }
  # my_formula_string <- paste0(as.character(my_formula)[2], as.character(my_formula)[1], as.character(my_formula)[3])
  my_formula_string <- paste0(as.character(my_formula)[1], as.character(my_formula)[3])
  print(paste0("FORMULA: ", my_formula_string))
  print(result)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Formula <- my_formula_string
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}



run_permdisp_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  dist_matrix <- vegdist(t(my_data), method = my_method)
  betadisper_object <- with(my_metadata, betadisper(dist_matrix, group = get(my_group)))
  permutest_results <- permutest(betadisper_object, permutations = permutations, parallel = 2)
  
  for (r in rownames(permutest_results$tab)){
    variable <- r
    Degrees_of_freedom <- permutest_results$tab[r,]$Df[1]
    SumOfSqs <- round(permutest_results$tab[r,]$`Sum Sq`[1],3)
    meanSqs <- round(permutest_results$tab[r,]$`Mean Sq`[1], 3)
    F.model <- round(permutest_results$tab[r,]$F[1], 3)
    N_permutations <- permutest_results$tab[r,]$N.Perm[1]
    p_value <- round(permutest_results$tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degrees_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       N_permutations,
                                                       p_value))
  }
  print(permutest_results)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","Permutations","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Group <- my_group
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}

run_anosim_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  anosim_object <- with(my_metadata, anosim(x = t(my_data), 
                                            grouping = get(my_group),
                                            permutations = permutations,
                                            distance = my_method,
                                            parallel = 2))
  
  stat_sig_table <- rbind(stat_sig_table, data.frame(my_group,
                                                     anosim_object$statistic,
                                                     anosim_object$signif,
                                                     anosim_object$permutations))
  names(stat_sig_table) <- c("Variable","R_statistic", "Significance","Permutations")
  stat_sig_table$Method <- my_method
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Correlation analysis

# calculate_correlation_matrix_stats <- function(correlation_matrix, method = "pearson", adjust = "none"){
#   cor_result <- corr.test(correlation_matrix, adjust = adjust, method = method)
#   cor_result$p
# }

calculate_correlation_matrix <- function(mydata, method = "pearson", adjust = "BH"){
  
  # Remove row entries that don't vary across all samples
  internal_data.m <- mydata
  zv <- apply(internal_data.m, 1, function(x) length(unique(x)) == 1)
  internal_data.m <- internal_data.m[!zv, ]
  
  # Take a two vectors and perform a significance/correlation test
  calculate_stats <- function(x, y, dist.name) {
    k <- cor.test(x, y, method=dist.name)
    c(k$estimate, k$stat, k$p.value)
  }
  
  # cor_result <- corr.test(t(internal_data.m), adjust = adjust, method = method)
  r <- apply(t(internal_data.m), 2,function(col) 
    t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,1,drop =F])
  p <- apply(t(internal_data.m), 2,function(col) 
    t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,3,drop =F])
  padj <- apply(t(internal_data.m), 2,function(col) 
    p.adjust(t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,3,drop =F],
             method = adjust))
  
  rownames(r) <- colnames(r)
  rownames(p) <- colnames(p)
  rownames(padj) <- colnames(padj)
  r <- signif(r, digits = 5)
  p <- signif(p, digits = 5)
  padj <- signif(padj, digits = 5)
  
  # r	The matrix of correlations
  # n	Number of cases per correlation
  # t	value of t-test for each correlation
  # p	two tailed probability of t for each correlation. For symmetric matrices, p values adjusted for multiple tests are reported above the diagonal.
  # se	standard error of the correlation
  # ci	the alpha/2 lower and upper values, as well as the (Holm or Bonferroni) adjusted confidence intervals.
  # list(cor_matrix = cor_result$r, cor_pval_matrix = cor_result$p)
  list(cor_matrix = r, cor_pval_matrix = p, cor_padj_matrix = padj)
}

# Calculates correlations between a feature and all other features. Works by row.
# Assumes feature is in the row, generally sample will be the column
calculate_feature_correlations <- function(mydata, feature, method = "pearson", adjust = "BH", filename = NULL){
  
  # Remove entries that don't vary across all samples
  internal_data.m <- mydata
  zv <- apply(internal_data.m, 1, function(x) length(unique(x)) == 1)
  internal_data.m <- internal_data.m[!zv, ]
  if (!feature %in% rownames(internal_data.m)){
    print(paste0("feature ", feature, " not found in input data. Possibly does not vary across samples?"))
    sig.df <- NULL
    return(sig.df)
  }
  
  # Take a two vectors and perform a signficance/correlation test
  calculate_stats <- function(x, y, dist.name) {
    k<-cor.test(x, y, method=dist.name);
    c(k$estimate, k$stat, k$p.value)
  }
  
  feature_data <- t(internal_data.m)[,feature]
  correlation_results.m <- apply(t(internal_data.m), 2, calculate_stats, feature_data, method)
  correlation_results.m <- t(correlation_results.m)
  pval_adj.col <- p.adjust(correlation_results.m[,3], adjust)
  correlation_results.m <- cbind(correlation_results.m, pval_adj.col)
  colnames(correlation_results.m) <- c("correlation", "t-stat", "p-value", "pval_adj");
  ord.inx <- order(correlation_results.m[,3]) # order by p-value
  sig.m <- signif(correlation_results.m[ord.inx,],5); # Round values
  
  # Remove self-comparison
  sig.m <- sig.m[!rownames(sig.m) == feature, ]
  
  # Write to file
  if (!is.null(filename)){
    write.csv(sig.m, file = filename)
  }
  
  # Remove Inf values from table
  sig.df <- as.data.frame(sig.m)
  is.na(sig.df) <- sapply(sig.df, is.infinite);
  sig.df[is.na(sig.df)] <- 0
  
  sig.df
}

# Given data matrix, plot correlations between a feature/variable and other features/variables
plot_feature_correlations <- function(mydata, feature, method = "pearson", adjust = "BH", top_n = 25,
                                      filename = NULL, plot_width = 10, plot_height = 10, format = "pdf"){
  # internal_data.m <- mydata
  # internal_data.m <- internal_data.m[apply(internal_data.m, 1, function(x) {length(which(x > 0))}) /length(internal_data.m[1,]) > 0.1,]
  
  correlation_result.df <- calculate_feature_correlations(mydata, feature, method = method, adjust = adjust)
  if(is.null(correlation_result.df)){
    return()
  }
  
  # First get most signficant correlations (p-value)
  ord.inx <- order(correlation_result.df[,3]);
  correlation_result.df <- correlation_result.df[ord.inx,] # Should be ordered already, though just ensure
  
  if(nrow(correlation_result.df) > top_n){
    correlation_result.df <- correlation_result.df[1:top_n, ];
  }
  
  # Then order by correlation direction
  ord.inx <- order(correlation_result.df[,1]);
  
  if(sum(correlation_result.df[,1] > 0) == 0){ # all negative correlation
    ord.inx <- rev(ord.inx);
  }
  correlation_result.df <- correlation_result.df[ord.inx,]
  if (!is.null(filename)){
    Cairo::Cairo(file = filename, unit="cm", dpi=100, width=plot_width, height=plot_height, type=format, bg="white",units = "cm");
  }
  plot_title <- paste("Top",nrow(correlation_result.df), "features correlated with", feature);
  
  cols <- ifelse(correlation_result.df[,1] > 0, "mistyrose","lightblue");
  
  par(xaxt = "n")
  xticks <- seq(-1,1,by = .2)
  dotchart(correlation_result.df[,1], 
           labels = rownames(correlation_result.df),
           pch="", 
           xlim=c(-1.1,1.1), 
           xlab="Correlation Coefficients", 
           # main=title,
           # yaxs = "n"
           cex = 0.7)
  par(xaxt = "s")
  axis(1, at=xticks, cex.axis=0.7) 
  title(main = plot_title, cex.main = 0.7)
  # rownames(correlation_result.m) <- NULL;
  barplot(correlation_result.df[,1], 
          space=c(0.5, rep(0, nrow(correlation_result.df)-1)),
          xlim=c(-1,1),
          xaxt="n", 
          col = cols, 
          add=T,
          horiz=T);
  
  for (row in 1:nrow(correlation_result.df)){
    offset <- ifelse(correlation_result.df[row,"correlation"] < 0, -0.1, 0.1)
    if (correlation_result.df[row,"pval_adj"] <= 0.05 & correlation_result.df[row,"pval_adj"] > 0.01){
      text(correlation_result.df[row,"correlation"] + offset, row, labels = "*")
    } else if (correlation_result.df[row,"pval_adj"] <= 0.01 & correlation_result.df[row,"pval_adj"] > 0.001){
      text(correlation_result.df[row,"correlation"] + offset, row, labels = "**")
    } else if (correlation_result.df[row,"pval_adj"] <= 0.001){
      text(correlation_result.df[row,"correlation"] + offset, row, labels = "***")
    }
  }
  if (!is.null(filename)){
    dev.off()
  }
}

# Given correlation matrix and (optional) p-value matrix, 
# plot correlations between a feature/variable and other features/variables
plot_feature_correlations_external <- function(cor_matrix, feature, p_value_matrix = NULL, top_n = 25,
                                               y_label_size = 0.6,
                                               filename = NULL, plot_width = 50, plot_height = 50, format = "pdf",
                                               include_self = T, include_title = T){
  
  cor_feature.m <- cor_matrix[,feature, drop =F]

  if (!is.null(p_value_matrix)){
    p_value_feature.m <- p_value_matrix[,feature,drop =F] 
    correlation_result.df <- as.data.frame(cbind(cor_feature.m,p_value_feature.m))
    names(correlation_result.df) <- c("correlation", "p_value")
  } else{
    correlation_result.df <- as.data.frame(cor_feature.m)
    names(correlation_result.df) <- c("correlation")
  }
  
  # First get most significant correlations (p-value)
  if (!is.null(p_value_matrix)){
    ord.inx <- order(correlation_result.df[,2])
    correlation_result.df <- correlation_result.df[ord.inx,,drop =F]
    if(nrow(correlation_result.df) > top_n){
      correlation_result.df <- correlation_result.df[1:top_n,,drop =F]
    }  
  }
  
  # Order by degree of correlation if no p-value matrix provided
  if (is.null(p_value_matrix)){
    ord.inx <- order(abs(correlation_result.df[,1]),decreasing = T)
    correlation_result.df <- correlation_result.df[ord.inx,,drop =F]
    if(nrow(correlation_result.df) > top_n){
      correlation_result.df <- correlation_result.df[1:top_n,,drop =F]
    }
  }
  
  # Order by correlation direction
  ord.inx <- order(correlation_result.df[,1])
  if(sum(correlation_result.df[,1] > 0) == 0){ # all negative correlation
    ord.inx <- rev(ord.inx);
  }
  correlation_result.df <- correlation_result.df[ord.inx,,drop =F]
 
  # Save plot to file
  if (!is.null(filename)){
    if (format == "pdf"){
      cairo_pdf(filename = filename,height = plot_height, width = plot_width)
      # dev.off()      
    } else if (format == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      # svg(filename = filename,height = plot_height, width = plot_width)
      svglite(file = filename,height = plot_height, width = plot_width)
      # dev.off()
    }
  }
  
  # if (!is.null(filename)){
  #   Cairo::Cairo(file = filename, unit="cm", dpi=200, width=plot_width, height=plot_height, type=format, bg="white",units = "cm")
  # }
  
  
  cols <- ifelse(correlation_result.df[,1] > 0, "mistyrose","lightblue");
  if (!include_self){
    correlation_result.df <- correlation_result.df[rownames(correlation_result.df) != feature,]  
  }

  par(xaxt = "n")
  xticks <- seq(-1,1,by = .2)
  dotchart(correlation_result.df[,1], 
           labels = rownames(correlation_result.df),
           pch="", 
           xlim=c(-1.1,1.1), 
           xlab="Correlation Coefficients", 
           # mgp = c(4,1,0), # margin of the axis title, labels and line
           # main=title,
           # yaxs = "n"
           cex = 0.7)
  par(xaxt = "s")
  axis(1, at=xticks, cex.axis=0.7,pos = c(0,0)) 
  
  if (include_title){
    plot_title <- paste("Top",nrow(correlation_result.df), "features correlated with", feature)
    title(main = plot_title, cex.main = 0.7)  
  }

  barplot(correlation_result.df[,1], 
          space=c(0.5, rep(0, nrow(correlation_result.df)-1)), 
          xlim=c(-1,1), 
          xaxt="n", 
          col = cols, 
          add=T,
          horiz=T)
  
  if (!is.null(p_value_matrix)){
    for (row in 1:nrow(correlation_result.df)){
      offset <- ifelse(correlation_result.df[row,"correlation"] < 0, -0.1, 0.1)
      if (correlation_result.df[row,"p_value"] <= 0.05 & correlation_result.df[row,"p_value"] > 0.01){
        text(correlation_result.df[row,"correlation"] + offset, row, labels = "*")
      } else if (correlation_result.df[row,"p_value"] <= 0.01 & correlation_result.df[row,"p_value"] > 0.001){
        text(correlation_result.df[row,"correlation"] + offset, row, labels = "**")
      } else if (correlation_result.df[row,"p_value"] <= 0.001){
        text(correlation_result.df[row,"correlation"] + offset, row, labels = "***")
      }
    }
  }
  
  if (!is.null(filename)){
    dev.off()
  }
}

# ------------------------------------------------
generate_correlation_network_from_count_data <- function(mydata, p_value_threshold = 0.05, cor_threshold = 0.5, method = "pearson", adjust = "BH"){
  correlation_results <- calculate_correlation_matrix(mydata, adjust = adjust, method = method)
  cor.m <- correlation_results$cor_matrix
  cor_pval.m <- correlation_results$cor_pval_matrix
  generate_correlation_network(cor.m, cor_pval.m,p_value_threshold = 0.05, cor_threshold = 0.5)
}

# Generate a correlation network
generate_correlation_network <- function(cor_matrix, p_matrix = NULL, p_value_threshold = 0.05, cor_threshold = 0.5,
                                         filename = NULL, relabeller_function = NULL,
                                         node_size = 4, node_colour = "grey20", node_fill = "grey20", node_shape = 21,
                                         node_label_segment_colour = "grey",
                                         label_colour = "black",label_size = 1,
                                         plot_height = 10, plot_width = 10, plot_title = "",
                                         edge_width_min = .3, edge_width_max = 1,edge_alpha =1,
                                         network_layout = "fr", exclude_to_from_df = NULL,
                                         myseed=NULL, edgetype = "link", show_p_label=F,
                                         file_type = "pdf"){
  if (!is.null(myseed)){
    set.seed(myseed)  
  }
  
  cor.m <- as.matrix(cor_matrix)
  
  if (max(cor.m[row(cor.m)!=col(cor.m)]) < cor_threshold){
    print("No entries left after filtering at correlation threshold")
    return()
  }
  if (!is.null(relabeller_function)){
    colnames(cor.m) <- lapply(colnames(cor.m), relabeller_function)
    rownames(cor.m) <- lapply(rownames(cor.m), relabeller_function)
  }
  
  if (!is.null(p_matrix)){
    if (min(p_matrix[row(p_matrix)!=col(p_matrix)]) > p_value_threshold){
      print("No entries left after filtering at p-value threshold")
      return()
    }
    cor_pval.m <- as.matrix(p_matrix)  
    if (!is.null(relabeller_function)){
      # colnames(cor_pval.m) <- relabeller_function(colnames(cor_pval.m))
      # rownames(cor_pval.m) <- relabeller_function(rownames(cor_pval.m))
      colnames(cor_pval.m) <- lapply(colnames(cor_pval.m), relabeller_function)
      rownames(cor_pval.m) <- lapply(rownames(cor_pval.m), relabeller_function)
    }
  }

  # Melt correlation matrix
  graph.df <- melt(cor.m,value.name = "Correlation",varnames = c("Variable_1", "Variable_2"))
  if (!is.null(p_matrix)){
    graph.df$P_value <- melt(cor_pval.m)$value
    graph.df <- graph.df[graph.df$P_value <= p_value_threshold,]  
    graph.df$P_value_label <- NA
    if (show_p_label==T){
      graph.df$P_value_label <- lapply(graph.df$P_value, function(x) ifelse(x <= 0.001, "***", ifelse(x <= 0.01, "**", ifelse(x <= 0.05, "*", ""))))  
    }
    # graph.df$P_value_label <- round(graph.df$Correlation,2) # For sanity testing
    
  } else{
    graph.df$P_value_label <- NA
    # graph.df$P_value_label <- round(graph.df$Correlation,2) # For sanity testing
    names(graph.df) <- c("Variable_1", "Variable_2", "Correlation", "P_value_label")
  }
  
  # Remove bidirectional links
  graph.df$ordered_label <- apply(graph.df, 1, function(x) {paste0(sort(c(as.character(x[1]), as.character(x[2]))),collapse = ":")})
  graph.df <- graph.df[!duplicated(graph.df$ordered_label),]
  
  # Remove edges with a correlation of 0, regardless what threshold
  graph.df <- graph.df[graph.df$Correlation != 0,]

  # Remove edges that are below the correlation threshold
  graph.df <- graph.df[abs(graph.df$Correlation) >= cor_threshold,]
  
  # Remove edges specified in supplied dataframe.
  if (!is.null(exclude_to_from_df)){
    graph.df <- graph.df[!paste0(graph.df$Variable_1, "-", graph.df$Variable_2) %in% paste0(edges_to_remove.df[,1], "-", edges_to_remove.df[,2]),]
    graph.df <- graph.df[!paste0(graph.df$Variable_1, "-", graph.df$Variable_2) %in% paste0(edges_to_remove.df[,2], "-", edges_to_remove.df[,1]),]
  }
  
  # Generate graph object and remove looped edges and isolated nodes
  graph.df <- 
    tidygraph::as_tbl_graph(graph.df) %>%
    # Remove loops
    tidygraph::activate(edges) %>%
    filter(!edge_is_loop()) %>%
    # filter(!edge_is_mutual()) %>%
    # filter(!edge_is_multiple()) %>%
    # Remove isolated nodes
    tidygraph::activate(nodes) %>%
    filter(!node_is_isolated())
  
  # Build plot
  set_graph_style(plot_margin = margin(1,1,1,1))
  correlation_graph_plot <- ggraph(graph.df, layout = network_layout)
  if (edgetype == "link"){
    correlation_graph_plot <- correlation_graph_plot + 
      geom_edge_link(aes(colour = Correlation, width=abs(Correlation), label = P_value_label),  show.legend = T, alpha = edge_alpha,
                     angle_calc = "along", label_colour = "black",
                     label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "fan"){
    correlation_graph_plot <- correlation_graph_plot +
    geom_edge_fan(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha,
                  angle_calc = "along", label_colour = "black",
                  label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "elbow"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_elbow(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha, strength = 1,
                      angle_calc = "along", label_colour = "black",
                      label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "bend"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_bend(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha, strength = 1,
                     angle_calc = "along", label_colour = "black",
                     label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "hive"){
    correlation_graph_plot <- correlation_graph_plot + 
      geom_edge_hive(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha, strength = 1,
                     angle_calc = "along", label_colour = "black",
                     label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  }
  correlation_graph_plot <- correlation_graph_plot +
    scale_edge_width_continuous(name="Correlation", range = c(edge_width_min,edge_width_max),
                                breaks = seq(-1,1,.2)) +
    scale_edge_colour_gradientn(colours = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                                     "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                     "#4393C3", "#2166AC", "#053061")))(11),
                          limits = c(-1,1), # limit colours to full corrrelation range
                          breaks = seq(-1,1,.2), # Break colours from -1 to 1 in steps of 0.2
                          guide = guide_edge_colourbar(barwidth = 0.5, barheight = 10)) +
    geom_node_point(colour = node_colour, fill = node_fill, pch = node_shape, size =node_size)  +
    geom_node_text(aes(label = name), colour = label_colour, 
                   size = label_size, nudge_y = -.01, #fontface = "bold",
                   point.padding = unit(0.3, "lines"),
                   segment.size = 0.3,
                   segment.colour = node_label_segment_colour,
                   # family = "serif",
                   repel = T) +
    ggtitle(label = plot_title) +
    theme_graph(background = "white")
    # Override edge colour to show coloured width, entries for width need to match the
    # breaks defined in the colour gradient (at least in length)
    # guides(edge_color = guide_legend(override.aes = list(edge_width = abs(seq(-1,1,.2)*1.2))),
    break_length <- length(seq(-1,1,.2))
    edge_widths <- abs(c(rev(rev(seq(-edge_width_max, -edge_width_min, length.out = break_length/2))[-1]),
                         rev(seq(edge_width_max, edge_width_min, length.out = break_length/2))))
    # edge_widths <- abs(seq(-1,1,.2))
    # edge_widths <- abs(edge_widths *seq(-1,1,.2))
    correlation_graph_plot <- correlation_graph_plot +
    guides(edge_color = guide_legend(override.aes = list(edge_width = edge_widths)),
           edge_width = F)
  
  # Save plot to file
  if (!is.null(filename)){
    if (file_type == "pdf"){
      cairo_pdf(filename = filename,height = plot_height, width = plot_width)
      plot(correlation_graph_plot)
      dev.off()      
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      # svg(filename = filename,height = plot_height, width = plot_width)
      svglite(file = filename,height = plot_height, width = plot_width)
      plot(correlation_graph_plot)
      dev.off()
    }
  }
  list(network_data = graph.df, network_plot = correlation_graph_plot)
}

# Generate a dot correlation plot generated by corrplot
plot_corrplot <- function(correlation_matrix, 
                          p_value_matrix = NULL, 
                          plot_title = "", 
                          plot_title_size = .6,
                          plot_height = 10, plot_width = 10,
                          p_value_threshold = 0.05, 
                          label_size = 1,
                          relabeller_function = NULL, 
                          insig = "blank", insig_pch = 4, insig_pch_cex = 1, insig_pch_col = "black",
                          make_insig_na = F,
                          to_exclude = NULL, 
                          method = "circle", 
                          outline = T,
                          label_colour = "black", 
                          colour_label_size = 1,
                          grid_colour = "black",
                          pairs_to_na = NULL, # Should be two column dataframe (row column / column row)
                          order = "hclust", 
                          col = NULL, 
                          filename = NULL,
                          file_type = "pdf"){
  
  cor.m <- correlation_matrix
  cor_pval.m <- p_value_matrix
  
  # Entries to remove completely
  if (!is.null(to_exclude)){ 
    cor.m <- cor.m[!rownames(cor.m) %in% to_exclude, !colnames(cor.m) %in% to_exclude]
    if (!is.null(p_value_matrix)){
      cor_pval.m <- cor_pval.m[!rownames(cor_pval.m) %in% to_exclude, !colnames(cor_pval.m) %in% to_exclude]  
    }
    
  }
  # Entries to make NA. Should be a two column dataframe
  if (!is.null(pairs_to_na)){
    if (order == "hclust"){
      # stop("Cannot use hclust ordering with NA values in matrix")
      # So before inserting NA values, first order the matrix
      # print(summary(hclust(dist(cor.m),method = "average")))
      ord <- hclust(dist(cor.m),method = "average")$order
      cor.m <- cor.m[ord,ord]
      if (!is.null(p_value_matrix)){
        cor_pval.m <- cor_pval.m[ord,ord]
      }
      order <- "original"
    }
    for (row in 1:nrow(pairs_to_na)){
      a <- as.character(pairs_to_na[row,1])
      b <- as.character(pairs_to_na[row,2])
      if (a %in% rownames(cor.m) & b %in% rownames(cor.m)){
        cor.m[a, b] <- NA
        cor.m[b, a] <- NA
        if (!is.null(p_value_matrix)){
          cor_pval.m[a, b] <- NA
          cor_pval.m[b, a] <- NA
        }
      }
    }
  }
  
  if (make_insig_na == T){
    cor.m[which(cor_pval.m >= p_value_threshold)] <- NA
    cor_pval.m[cor_pval.m >= p_value_threshold] <- NA
  }
  
  if (!is.null(relabeller_function)){
    colnames(cor.m) <- lapply(colnames(cor.m), relabeller_function)
    rownames(cor.m) <- lapply(colnames(cor.m), relabeller_function)
    if (!is.null(p_value_matrix)){
      colnames(cor_pval.m) <- lapply(colnames(cor_pval.m), relabeller_function)
      rownames(cor_pval.m) <- lapply(rownames(cor_pval.m), relabeller_function)
      # colnames(cor_pval.m) <- relabeller_function(colnames(cor_pval.m))
      # rownames(cor_pval.m) <- relabeller_function(rownames(cor_pval.m))
    }
  }
  if (is.null(col)){
    col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                  "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061")))(200)
  }
  
  if (!is.null(filename)){
    # pdf(file = filename,height = plot_height, width = plot_width)
    if (file_type == "pdf"){
      pdf(file = filename, height = plot_height, width = plot_width)
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      # svg(filename = filename,height = plot_height, width = plot_width)
      svglite(file = filename,height = plot_height, width = plot_width)
    }
  }
  
  corrplot(corr = cor.m,
           method = method,
           outline = outline,
           tl.col = label_colour,
           tl.cex = label_size,
           addgrid.col = grid_colour,
           # tl.srt = 45,
           # title = plot_title,
           col = col,
           # col = brewer.pal(n = 8, name = "RdYlBu"),
           type = "lower",
           diag = F,
           na.label = "square",
           na.label.col = "grey",
           order = order,
           hclust.method = "average",
           p.mat = cor_pval.m,
           sig.level = p_value_threshold,
           # insig = "blank",
           insig = insig,
           pch = insig_pch,
           pch.cex = insig_pch_cex,
           pch.col = insig_pch_col,
           cl.pos = 'r',
           cl.cex = colour_label_size,
           
           mar=c(1,0,3,1))
  title(main = plot_title,cex.main = plot_title_size)
  if (!is.null(filename)){
    dev.off()
  }
}

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

library(ggsignif)
generate_significance_boxplots <- function(mydata.df, # Main dataframe
                                           variable_column, # Variable to group by
                                           value_column, # Values to use in boxplots
                                           variable_colours_available = F, # requires a column named "variable_column_colour"
                                           significances.df, # Dataframe containing significance values for each group comparison
                                           p_value_column, # Column in significances.df containing p-values
                                           fill_palette = NULL, # Custom fill palette.
                                           sig_threshold = 0.05, # p-values less than this will be displayed
                                           sig_line_scaling_percentage = 0.08, # percent of max y value for spacing significance lines
                                           sig_vjust = 0.03, # amount to vertically adjust the significance annotations
                                           sig_tip_length = 0.01, # length of tips on significance lines
                                           sig_linetype = 1, # linetype of significance lines
                                           sig_colour = "grey20", # colour of significance lines
                                           sig_line_starting_scale = 1.05
){
  # Requires ggplot and ggsignif packages
  
  # mydata.df must have the following columns: "variable_column" , "value_column" and (optionally) "variable_column_colour"
  # significances.df must have the following columns : "Variable", "Group_1", "Group_2"
  # p_value_column must be defined by the user and also be in the dataframe
  variable_values <- levels(mydata.df[[variable_column]])
  if (variable_colours_available == T){
    color_col_name <- paste0(variable_column, "_colour")
    variable_colours <- setNames(as.character(unique(mydata.df[[color_col_name]])), as.character(unique(mydata.df[[variable_column]])))
  } else{
    if (is.null(fill_palette)){
      # my_colour_palette_30_distinct
      internal_colour_palette <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4",
                                   "#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8",
                                   "#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d",
                                   "#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a",
                                   "#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
    } else{
      internal_colour_palette <- fill_palette
    }
    variable_colours <- setNames(internal_colour_palette[1:length(variable_values)], variable_values)  
  }
  
  # Subset significances to the entries for variable of interest variable
  sig_subset.df <- subset(significances.df, Variable == variable_column) %>% dplyr::select(Variable, Group_1, Group_2,p_value_column)
  
  # Filter to siginificant values
  sig_subset.df <- sig_subset.df[which(sig_subset.df[,p_value_column] < sig_threshold),]
  
  # Generate plot
  myplot <- ggplot(mydata.df, aes(x = get(variable_column), fill = get(variable_column), y = get(value_column))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=0.5, width = 0.10, height=0) +
    guides(fill=FALSE) +
    # scale_y_continuous(breaks = seq(0, max(sig_subset.df$y_position)+.1, 1)) +
    scale_fill_manual(values = variable_colours, name = variable_column) +
    xlab(variable_column) +
    ylab(value_column) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = 0.5),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "white", colour = "white", size = 1),
          strip.text = element_text(size = 6),
          legend.key=element_blank(),
          legend.direction="vertical",
          legend.background = element_rect(colour ="white", size = .3),
          legend.text.align = 0,
          legend.title = element_text(size=10, face="bold"),
          legend.title.align = 0.5,
          legend.margin = margin(c(2,2,2,2)),
          legend.key.height=unit(.4,"cm"),
          legend.text = element_text(size = 8),
          axis.text = element_text(size = 8, colour = "black"),
          axis.text.x = element_text(angle = 0, vjust = .5),
          axis.title = element_text(size = 10,face = "bold"),
          complete = F,
          plot.title = element_text(size = 6,hjust = 0.5))
  
  # Return the plot without siginificances if there are no significant comparisons
  if (dim(sig_subset.df)[1] == 0){
    print("No sigificant comparisons")
    return(myplot)
  }
  
  # Create p-value label
  sig_subset.df$P_value_label <- as.character(lapply(sig_subset.df[,p_value_column], function(x) ifelse(x <= 0.001, "***", 
                                                                                                        ifelse(x <= 0.01, "**", 
                                                                                                               ifelse(x <= 0.05, "*", "ns")))))
  # Determine the maximum diversity value for the pair of groups being compared
  for (row in 1:nrow(sig_subset.df)){
    group_1 <- as.character(sig_subset.df[row,"Group_1"])
    group_2 <- as.character(sig_subset.df[row,"Group_2"])
    y_max <- max(mydata.df[which(mydata.df[,variable_column] == group_1),][,value_column],
                 mydata.df[which(mydata.df[,variable_column] == group_2),][,value_column])
    sig_subset.df[row,"y_max"] <- y_max
    sig_subset.df[row,"level_index_group_1"] <- which(levels(mydata.df[,variable_column]) == group_1)
    sig_subset.df[row,"level_index_group_2"] <- which(levels(mydata.df[,variable_column]) == group_2)
    sig_subset.df[row,"level_distance"] <- abs(sig_subset.df[row,"level_index_group_2"] - sig_subset.df[row,"level_index_group_1"])
  }
  
  # Define the y position of the significance lines to ensure no overlaps.
  # This is based on the maximum value multiplied by a small scaling value
  # that gradually increases.
  
  # sig_subset.df <- sig_subset.df[order(sig_subset.df$y_max),]
  # sig_subset.df <- sig_subset.df[order(sig_subset.df$level_index_group_1),]
  # sig_subset.df <- sig_subset.df[order(sig_subset.df$level_index_group_2),]
  sig_subset.df <- sig_subset.df[order(sig_subset.df$level_distance),]
  scale <- sig_line_starting_scale # starting scale
  for (row in 1:nrow(sig_subset.df)){
    sig_subset.df[row, "y_position"] <- max(sig_subset.df[, "y_max"]) * scale
    scale <- scale + sig_line_scaling_percentage # increase scale value
  }
  
  myplot <- myplot + ggsignif::geom_signif(data = sig_subset.df,
                                           manual = T,
                                           inherit.aes = F,
                                           aes(xmin = Group_1, xmax = Group_2, annotations = P_value_label, y_position = y_position),
                                           linetype = sig_linetype,
                                           color = sig_colour,
                                           size = .5,
                                           tip_length = sig_tip_length, vjust = sig_vjust)
  myplot
}