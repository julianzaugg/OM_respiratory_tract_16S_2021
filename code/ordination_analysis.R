detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

library(vegan)
library(ggplot2)
# library(ggfortify)

############################################################
# Various colour palettes
my_colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
my_colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
my_colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
my_colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
my_colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
my_colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
my_colour_palette_32_distinct <- c("#ea7e00","#ca0074","#d1c69b","#474007","#bb00ad","#9c80ff","#be3300","#542e72","#00b9f5","#09436b","#8b0036","#9ac8e6","#ff1059","#959eff","#154a11","#0290f4","#ff7762","#7dbf00","#ff8194","#834c00","#006e73","#f9bb5d","#d6c943","#017229","#00d3a8","#732427","#36e191","#6a8200","#efb3ea","#3227bb","#ff90e1","#e92a12")
# lesion_palette_7 <- c("#8558d6","#6ee268","#d247ad","#c9d743","#d7453e","#59a237","#d78f2a")
# patient_palette_45 <- c("#d64530","#585fb1","#795d97","#9e4773","#3f6921","#71692c","#a2b93c","#d571cc","#9b3e97","#33947a","#98ad66","#448a4e","#869ae0","#5ce7af","#e085a3","#dfdc87","#d19be2","#5cb735","#e38269","#3db6c0","#50b565","#50902c","#a98a2c","#dde84a","#db3d76","#5fe485","#7c8329","#b3e791","#6fe965","#5ebce9","#3c86c1","#2a6a45","#65b688","#6651d1","#af4ed3","#df872f","#56e4db","#737cea","#ac464b","#dd37b5","#995b2b","#daac6f","#92e2be","#a2e24b","#e0be3a")
patient_palette_270 <- c("#456c00","#77bbff","#75fb3f","#273300","#f5006f","#ac008b","#125700","#ffef77","#00278e","#3d005e","#d84100","#015686","#01289f","#ff8c4c","#0070b5","#8015cd","#feffd6","#02d8aa","#019cf4","#4f2f00","#bbffe9","#c52900","#1b002c","#a3ac00","#5d9200","#f29fff","#231500","#934cfc","#988a00","#002cbb","#ffeb36","#ffa758","#f1f9ff","#000045","#b4004b","#602900","#390048","#e6c400","#00ce3c","#ff7bd0","#8cff56","#e60051","#b89aff","#00474b","#d5fbff","#ff79c2","#1d0016","#00635d","#ff8e33","#992300","#ff6e91","#ffa081","#534a00","#61002d","#ffe1c1","#8c0072","#00405d","#89ffc6","#607500","#64ff6f","#002e52","#9b97ff","#b1ccff","#02c5cd","#5dffba","#beff45","#00112b","#b8ff74","#7f0027","#0074cd","#005c6f","#3f00af","#dd7900","#cced00","#77ffd6","#ffc5b5","#99ffb1","#01ea72","#f0ff88","#007f63","#abff9d","#391200","#003a74","#114b00","#0a0100","#ff5fff","#ffccfb","#00d6b7","#c7ff93","#1efa5a","#005437","#f6af00","#a60024","#ffb7e6","#ea0043","#c7ffbc","#72ab00","#789300","#585500","#c3ff14","#00f398","#ab4a00","#9b7600","#85e5ff","#006235","#130020","#006825","#ff735c","#007a7f","#02a3a4","#4856ff","#bf52ff","#00edbc","#a31bd6","#009642","#e93bee","#e400ae","#ffbdd2","#00cfc7","#f1ffaa","#009b7a","#dd00c9","#ff697d","#004a14","#ff72ac","#ff3d1f","#fffaa3","#5d0400","#027ba4","#01c774","#002655","#00941f","#0a32d7","#82acff","#ff8af3","#ff4165","#001104","#ffd6f2","#efebff","#aebc00","#3e0030","#c5abff","#00402e","#ff4bae","#0275f1","#be89ff","#ffd899","#00c765","#01b208","#97ffd4","#7e9fff","#00fde1","#0050c9","#ff8eb5","#c800cd","#005173","#ff2b95","#76ff7a","#ea0074","#001d70","#009856","#f100a8","#ba6b00","#0293df","#00462d","#ff6862","#f6ff65","#02bbda","#2c2200","#01a876","#e35a00","#e3000f","#ff819e","#5a0039","#a558ff","#e2ffb2","#784800","#016def","#b400a2","#00143c","#00212d","#403d00","#ff75fe","#975300","#166c00","#260008","#917fff","#ff8d89","#01bf7a","#ffa6bf","#800086","#90a100","#cce4ff","#dad800","#52c900","#46a700","#0c0039","#0b0052","#79009d","#003c85","#bb0034","#59e7ff","#af0064","#64001e","#c0007e","#000897","#bd8400","#2b007f","#318400","#31f1ff","#7c8600","#807300","#ffc072","#6f005f","#770040","#e62c00","#2e0019","#005599","#6535e1","#5b0099","#006bd5","#0142a1","#baaf00","#00ab2d","#ffcc40","#edffec","#ef0031","#153100","#abe9ff","#6bbd00","#e5ff4e","#ffdb43","#ffa5ef","#01c4f3","#ffbd8f","#84d100","#bbff84","#9fcdff","#7b0059","#ffe897","#ff8711","#ffa869","#febaff","#20003a","#94002b","#5387ff","#756dff","#fff494","#a5c1ff","#e0ffcf","#002417","#530076","#ff8459","#ffe4ec","#00b650","#0119b7","#c963ff","#a2ff64","#9c6800","#03b6f8","#00a0c2","#00240b","#6297ff","#bd0010","#fff7af","#7d2d00","#cf7aff","#af5600","#322c00","#500028")
my_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")


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

######################## Functions #########################

############################################################

source("code/helper_functions.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Convert variables to factors
discrete_variables <- c("Nose","Tympanic_membrane", "Otitis_Status",
                        "Season","Community","Gold_Star", "No_peop_res_discrete",
                        "H.influenzae_culture","M.catarrhalis_culture","S.pneumoniae_culture",
                        "Otitis_Status__Gold_Star", "Tympanic_membrane__Gold_Star",
                        "Community__Season","Community__Gold_Star","Community__Otitis_Status",
                        "H.Influenzae_qPCR", "M.catarrhalis_qPCR", "S.pneumoniae_qPCR",
                        "Corynebacterium_pseudodiphtheriticum","Dolosigranulum_pigrum","N_HRV")

metadata.df$Tympanic_membrane[metadata.df$Tympanic_membrane == "Unable to visualise/Not examined"] <- NA
# Remove AOM, just make values NA
metadata.df[metadata.df$Otitis_Status == "Acute Otitis Media","Otitis_Status"] <- NA

metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)



# Need to factorise the colour columns as well
colour_columns <- names(metadata.df)[grepl("colour", names(metadata.df))]
metadata.df[colour_columns] <- lapply(metadata.df[colour_columns], factor)

# Set order of Otitis Status
metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status , levels = c("Never OM", "HxOM", "Effusion","Perforation"))
# Set order of Season
metadata.df$Season <- factor(metadata.df$Season, levels = c("Spring", "Winter", "Autumn"))

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the counts
otu.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts.csv", header =T, row.names = 1))
genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

# Remove samples from the OTU table that are not in the filtered metadata
otu.m <- otu.m[,colnames(otu.m) %in% c("OTU.ID", as.character(metadata.df$Index))]
genus.m <- genus.m[,colnames(genus.m) %in% c("OTU.ID", as.character(metadata.df$Index))]

# Order the metadata.df by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]

# Order the matrices and metadata to be the same order
metadata.df <- metadata.df[order(rownames(metadata.df)),]
otu.m <- otu.m[,order(rownames(metadata.df))]
genus.m <- genus.m[,order(rownames(metadata.df))]

# CLR transform the otu matrix.
otu_clr.m <- clr(otu.m)
genus_clr.m <- clr(genus.m)


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Ordination analysis
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
                  paste(phylostring[5], phylostring[6], sep = ";")
                  # paste(phylostring[6], sep = ";")
                }))
}


combined_otu_labeller <- function(x){
  # print(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
  first_resolved_taxonomy(as.character(otu_taxonomy_map.df[otu_taxonomy_map.df$OTU.ID == x,]$taxonomy_species))
}


# --------------------
# Generate PCA plots.

### Generate ordination objects
# All samples
otu_pca <- rda(t(otu_clr.m), data = metadata.df)
genus_pca <- rda(t(genus_clr.m), data = metadata.df)

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ----------------------------------------- TESTING ----------------------------------------------
# remotes::install_github("gavinsimpson/ggvegan")
# library(tidyverse)
# library(vegan)
# library(ggvegan)
# PCAscores <- scores(genus_pca, display = "sites") %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Index") %>% 
#   full_join(metadata.df, by = "Index")
# 
# PCAvect <- scores(genus_pca, display = "species") %>% 
#   as.data.frame()
# PCAvect <- PCAvect[grepl("g__Staph|g__Haemo|g__Orni|g__Strep|g__Coryn|g__Dolosi", rownames(PCAvect)),]
# 
# clean_background <- theme(plot.background = element_rect("white"),
#                           panel.background = element_rect("white"),
#                           panel.grid = element_line("white"),
#                           axis.line = element_line("gray25"),
#                           axis.text = element_text(size = 12, color = "gray25"),
#                           axis.title = element_text(color = "gray25"),
#                           legend.text = element_text(size = 12),
#                           legend.key = element_rect("white"))
# plot_PCA <- ggplot() +
#   geom_point(data = PCAscores, aes(x = PC1, y = PC2, color = Otitis_Status)) +
#   # scale_color_manual(values = pal) +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   geom_segment(data = PCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
#                arrow = arrow(length = unit(0.2, "cm")), colour = "grey") +
#   geom_text(data = PCAvect, aes(x = PC1, y = PC2, label = rownames(PCAvect)),size = 2) +
#   clean_background +
#   labs(x = "PC1 (23.57%)",
#        y = "PC2 (12.23%)",
#        title = "Principal Components Analysis") 
# plot_PCA
# 
# 
# PCA_fortify <- fortify(genus_pca)
# PCA_fort_sites <- PCA_fortify %>% 
#   filter(Score == "sites") %>% 
#   full_join(., metadata.df, by = c("Label" = "Index"))
# PCA_fort_species <- PCA_fortify %>% 
#   filter(Score == "species")
# 
# PCA_fortify_plot <- ggplot() +
#   geom_point(data = PCA_fort_sites, aes(x = PC1, y = PC2, col = Otitis_Status)) +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   # scale_color_manual(values = c("coral", "lightgreen", "darkblue")) +
#   geom_segment(data = PCA_fort_species, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
#                arrow = arrow(length = unit(0.2, "cm"))) +
#   geom_text(data = PCA_fort_species, aes(x = PC1, y = PC2, label = Label)) +
#   clean_background +
#   labs(x = "PC1 (23.57%)",
#        y = "PC2 (12.23%)",
#        title = "Principal Components Analysis - using fortify()")
# PCA_fortify_plot


# temp <- metadata.df[,c("Otitis_Status", "Nose", "Season","No_peop_res_discrete")]
# temp <- temp[!rownames(temp) %in% get_samples_missing_data(temp,c("Otitis_Status", "Nose", "Season","No_peop_res_discrete")),]

# genus_NMDS <- metaMDS(t(genus_clr.m[,rownames(temp)]),distance = "euclidean",)
# stressplot(genus_NMDS)

# fit <- envfit(genus_NMDS, temp[,c("Otitis_Status", "Nose", "Season","No_peop_res_discrete")], perm = 999) 
# fit$vectors
# extract p-values for each species
# fit_pvals <- fit$vectors$pvals %>% 
#   as.data.frame() %>% 
#   rownames_to_column("species") %>% 
#   dplyr::rename("pvals" = ".")
# 
# # extract coordinates for species, only keep species with p-val = 0.001
# fit_spp <- fit %>% 
#   scores(., display = "vectors") %>% 
#   as.data.frame() %>% 
#   rownames_to_column("species") %>% 
#   full_join(., fit_pvals, by = "species") %>% 
#   filter(pvals == 0.001)

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

otu_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/OTU_counts_abundances_and_metadata.csv",header = T)
genus_data.df <- read.csv("Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv",header = T)

temp <- calculate_PC_abundance_correlations(genus_pca, mydata.df = genus_data.df,
                                            taxa_column = "taxonomy_genus",
                                            variables = discrete_variables)
write.csv(x = temp,
          file = "Result_tables/ordination_genus_PC_correlations.csv",
          quote =F,
          row.names = F)
# ------------------------------------------------------------------------------------------------
#           PUBLICATION

#### Within community
rural_metadata.df <- subset(metadata.df, Community == "Rural")
remote_metadata.df <- subset(metadata.df, Community == "Remote")
genus_rural_pca <- rda(t(genus_clr.m[,rownames(rural_metadata.df)]), data = rural_metadata.df)
genus_remote_pca <- rda(t(genus_clr.m[,rownames(remote_metadata.df)]), data = remote_metadata.df)

source("code/helper_functions.R")
file_type <- "svg"

# Genus, Otitis_Status, rural
generate_pca_plot(pca_object = genus_rural_pca,
                  my_metadata.df = rural_metadata.df,
                  plot_height = 7,
                  plot_width = 7,
                  variable_colours_available = T,
                  variable_to_plot = "Otitis_Status",
                  point_line_thickness = 0.7,
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = T,
                  hull_alpha = .7,
                  axis_limits = c(-6,3.5,-3,4.5),
                  
                  include_legend = T,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -6,
                  legend_y = 4,
                  legend_title = "Otitis status",
                  
                  plot_arrows = T,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  # arrow_label_colour = "midnightblue",
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =F,
                  file_type = file_type,
                  filename = paste0("Result_figures/pca_plots/genus_within_community/genus_Otitis_Status__Community_rural.", file_type)
)

# Genus, Otitis_Status, remote
generate_pca_plot(pca_object = genus_remote_pca,
                  my_metadata.df = remote_metadata.df,
                  plot_height = 7,
                  plot_width = 7,

                  variable_colours_available = T,
                  variable_to_plot = "Otitis_Status",
                  point_line_thickness = 0.7,
                  
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = T,
                  hull_alpha = .7,
                  axis_limits = c(-4.5,6,-4,3),
                  
                  include_legend = T,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -4.5,
                  legend_y = 3,
                  legend_title = "Otitis status",
                  
                  plot_arrows = T,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  # arrow_label_colour = "midnightblue",
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =F,
                  file_type = file_type,
                  filename = paste0("Result_figures/pca_plots/genus_within_community/genus_Otitis_Status__Community_remote.", file_type)
)


# ARROWS
generate_pca_plot(pca_object = genus_pca,
                  my_metadata.df = metadata.df,
                  plot_height = 7,
                  plot_width = 7,
                  
                  variable_colours_available = T,
                  variable_to_plot = "No_peop_res_discrete",
                  point_line_thickness = 0.7,
                  point_size = 0,
                  
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = F,
                  hull_alpha = .7,
                  axis_limits = c(-4.5,5,-4.5,3),
                  
                  include_legend = F,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -4.5,
                  legend_y = 3,
                  legend_title = "",
                  
                  plot_arrows = T,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =T,
                  file_type = "svg",
                  filename = paste0("Result_figures/pca_plots/genus/genus_arrows.svg")
)


# House hold size
generate_pca_plot(pca_object = genus_pca,
                  my_metadata.df = metadata.df,
                  plot_height = 7,
                  plot_width = 7,
                  
                  variable_colours_available = T,
                  variable_to_plot = "No_peop_res_discrete",
                  point_line_thickness = 0.7,
                  
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = T,
                  hull_alpha = .7,
                  axis_limits = c(-4.5,5,-4.5,3),
                  
                  include_legend = T,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -4.5,
                  legend_y = 3,
                  legend_title = "Household size",
                  
                  plot_arrows = F,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =F,
                  file_type = file_type,
                  filename = paste0("Result_figures/pca_plots/genus/genus_No_people_in_household.", file_type)
)

# Genus, Community
generate_pca_plot(pca_object = genus_pca,
                  my_metadata.df = metadata.df,
                  plot_height = 7,
                  plot_width = 7,
                  
                  variable_colours_available = T,
                  variable_to_plot = "Community",
                  point_line_thickness = 0.7,
                  
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = T,
                  hull_alpha = .7,
                  axis_limits = c(-4.5,5,-4.5,3),
                  
                  include_legend = T,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -4.5,
                  legend_y = 3,
                  legend_title = "Community",
                  
                  plot_arrows = F,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =F,
                  file_type = file_type,
                  filename = paste0("Result_figures/pca_plots/genus/genus_Community.", file_type)
)

# Genus, Season
generate_pca_plot(pca_object = genus_pca,
                  my_metadata.df = metadata.df,
                  plot_height = 7,
                  plot_width = 7,
                  
                  variable_colours_available = T,
                  variable_to_plot = "Season",
                  point_line_thickness = 0.7,
                  
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = T,
                  hull_alpha = .7,
                  axis_limits = c(-4.5,5,-4.5,3),
                  
                  include_legend = T,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -4.5,
                  legend_y = 3,
                  legend_title = "Season",
                  
                  plot_arrows = F,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =F,
                  file_type = file_type,
                  filename = paste0("Result_figures/pca_plots/genus/genus_Season.", file_type)
)

# Genus, Otitis_Status
generate_pca_plot(pca_object = genus_pca,
                  my_metadata.df = metadata.df,
                  plot_height = 7,
                  plot_width = 7,
                  
                  variable_colours_available = T,
                  variable_to_plot = "Otitis_Status",
                  point_line_thickness = 0.7,
                  
                  use_shapes = T,
                  plot_spiders = F,
                  label_sites = F,
                  # label_spider = T,
                  plot_hulls = T,
                  hull_alpha = .7,
                  axis_limits = c(-4.5,5,-4.5,3),
                  
                  include_legend = T,
                  legend_columns = 1,
                  legend_cex = .8,
                  legend_x = -4.5,
                  legend_y = 3,
                  legend_title = "Otitis Status",
                  
                  plot_arrows = F,
                  num_top_species = 3,
                  arrow_colour = "grey20",
                  arrow_scalar = 1.5,
                  arrow_label_colour = "royalblue4",
                  arrow_label_font_type = 1,
                  arrow_label_size = .7,
                  arrow_thickness = 1,
                  arrow_alpha = .7,
                  specie_labeller_function = first_resolved_taxonomy,
                  hide_grid =F,
                  file_type = file_type,
                  filename = paste0("Result_figures/pca_plots/genus/genus_Otitis_Status.", file_type)
)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# PERMANOVA

# Permutational Multivariate Analysis of Variance (PERMANOVA) can be used to 
# determine if the structure of the microbial communities is significantly different between
# environmental variables. This is done using the adonis function from Vegan with a distance metric, e.g. Bray-Curtis, and a
# specified number of permutations, e.g. 10,000.
# The analysis measures the degree each environmental variable affects the community composition and indicates 
# the significance of that effect on beta diversity (described by p-values and R2 values). 
# The R2 value corresponds to the proportion of variability observed in the dissimilarity.

# Genus, clr euclidean
print("Centred-log ratio transformed counts - Euclidean distance")
otu_permanova_results <- data.frame()
genus_permanova_results <- data.frame()

otu_within_community_permanova_results <- data.frame()
genus_within_community_permanova_results <- data.frame()

discrete_variables <- c("Nose", "Gold_Star", "Season", "Community", "N_HRV", "Dolosigranulum_pigrum", "Otitis_Status", "No_peop_res_discrete")

for (myvar in discrete_variables){
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  
  otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  
  otu_permanova_results <- rbind(otu_permanova_results,run_permanova_custom(my_metadata = metadata_subset.df, 
                                                      my_formula = as.formula(paste0("t(otu_clr_subset.m)~", myvar)),
                                                      my_method = "euclidean",label = "CLR",permutations = 9999))
  
  genus_permanova_results <- rbind(genus_permanova_results,run_permanova_custom(my_metadata = metadata_subset.df, 
                                                      my_formula = as.formula(paste0("t(genus_clr_subset.m)~", myvar)),
                                                      my_method = "euclidean",label = "CLR",permutations = 9999))
  
  if (myvar == "Community") {next}
  for (community in unique(metadata.df$Community)){
    
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Community = community)
    
    otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
    
    temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                 my_formula = as.formula(paste0("t(otu_clr_subset.m)~", myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 9999)
    temp$Community <- community
    otu_within_community_permanova_results <- rbind(otu_within_community_permanova_results, temp)
    
    temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                         my_formula = as.formula(paste0("t(genus_clr_subset.m)~", myvar)),
                         my_method = "euclidean",label = "CLR",permutations = 9999)
    temp$Community <- community
    genus_within_community_permanova_results <- rbind(genus_within_community_permanova_results,temp)

  }
}

write.csv(otu_permanova_results, file = "Result_tables/stats_various/otu_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_permanova_results, file = "Result_tables/stats_various/genus_PERMANOVA.csv", row.names = F, quote = F)

write.csv(otu_within_community_permanova_results, file = "Result_tables/stats_various/otu_within_community_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_within_community_permanova_results, file = "Result_tables/stats_various/genus_within_community_PERMANOVA.csv", row.names = F, quote = F)



# ---------------------------------------------
# PERMDISP (betadisper)
# See: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
# "Test the homogeneity of within-group multivariate dispersions on the basis of any resemblance measure."

otu_beta_diversity_dist.m <- as.matrix(vegdist(t(otu_clr.m), method = "euclidean",upper = F))
genus_beta_diversity_dist.m <- as.matrix(vegdist(t(genus_clr.m), method = "euclidean",upper = F))
write.csv(otu_beta_diversity_dist.m, file = "Result_tables/stats_various/otu_betadiversities.csv",row.names = T, quote = F)
write.csv(genus_beta_diversity_dist.m, file = "Result_tables/stats_various/genus_betadiversities.csv",row.names = T, quote = F)

otu_permdisp_results <- data.frame()
genus_permdisp_results <- data.frame()

otu_within_community_permdisp_results <- data.frame()
genus_within_community_permdisp_results <- data.frame()

for (myvar in discrete_variables){
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  
  otu_permdisp_results <- rbind(otu_permdisp_results, run_permdisp_custom(my_metadata = metadata_subset.df, 
                                                                          my_data = otu_clr_subset.m,
                                                                          my_group = myvar,
                                                                          my_method = "euclidean",
                                                                          permutations = 9999,
                                                                          label = "CLR"))
  
  genus_permdisp_results <- rbind(genus_permdisp_results, run_permdisp_custom(my_metadata = metadata_subset.df, 
                                                                              my_data = genus_clr_subset.m,
                                                                              my_group = myvar,
                                                                              my_method = "euclidean",
                                                                              permutations = 9999,
                                                                              label = "CLR"))
  
  if (myvar == "Community") {next}
  for (community in unique(metadata.df$Community)){
    
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Community == community)
    
    otu_clr_subset.m <- otu_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]

    
    temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                                my_data = otu_clr_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 9999,
                                label = "CLR")
    temp$Remote_Community <- community
    otu_within_community_permdisp_results <- rbind(otu_within_community_permdisp_results, temp)
    
    temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                                my_data = genus_clr_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 9999,
                                label = "CLR")
    temp$Remote_Community <- community
    genus_within_community_permdisp_results <- rbind(genus_within_community_permdisp_results, temp)
    
  }
  
}

write.csv(otu_permdisp_results, file = "Result_tables/stats_various/otu_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_permdisp_results, file = "Result_tables/stats_various/genus_PERMDISP.csv", row.names = F, quote = F)

write.csv(otu_within_community_permdisp_results, file = "Result_tables/stats_various/otu_within_community_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_within_community_permdisp_results, file = "Result_tables/stats_various/genus_within_community_PERMDISP.csv", row.names = F, quote = F)



