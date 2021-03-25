# Diversity calculations Alpha (Shannon, Chao1, Simpson) diversity for each sample.

# "Alpha diversity measures the diversity within a single sample and is generally 
# based on the number and relative abundance of taxa at some rank (e.g. species or OTUs).
# Beta diversity also uses the number of relative abundance of taxa at some rank, but measures 
# variation between samples. In other words, an alpha diversity statistic describes a single 
# sample and a beta diversity statistic describes how two samples compare."

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()
library(vegan)
library(reshape2)
library(ggplot2)
library(dplyr)

library(FSA)
library(phyloseq)
# library(nlme)

common_theme <- theme(
  panel.border = element_blank(), 
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
  axis.text = element_text(size = 9, colour = "black"),
  axis.title = element_text(size = 10,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))


generate_p_labels <- function(sig_table){
  for (sig_column in c("Chao1_Dunn_padj", "Shannon_Dunn_padj", "Simpson_Dunn_padj")){
    metric = strsplit(sig_column, "_")[[1]][1]
    sig_table[,paste0(metric, "_p_label")] <-
      as.character(lapply(sig_table[,sig_column], 
                          function(x) ifelse(x <= 0.001, "***", 
                                             ifelse(x <= 0.01, "**",
                                                    ifelse(x <= 0.05, "*", "ns")))))
  }
  sig_table
}


source("code/helper_functions.R")



# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Remove AOM, just make values NA
metadata.df[metadata.df$Otitis_Status == "Acute Otitis Media","Otitis_Status"] <- NA
metadata.df <- metadata.df[!is.na(metadata.df$Otitis_Status),]

# Set the Index to be the row name
rownames(metadata.df) <- metadata.df$Index

discrete_variables <- c("Community", "Nose", "Otitis_Status", "Season", "No_peop_res_discrete")

metadata.df$Tympanic_membrane[metadata.df$Tympanic_membrane == "Unable to visualise/Not examined"] <- NA

# Combine each variable with Community
for (x in discrete_variables){
  if (!grepl("Community", x) & !grepl("__", x)){
    metadata.df[,paste0(x,'__Community')] <-  with(metadata.df, paste0(get(x), "__",Community))
  }
}
discrete_variables <- c(discrete_variables, grep(".*__Community", colnames(metadata.df), value =T))

# Load the counts
otu.m <- as.matrix(read.csv("Result_tables/count_tables/OTU_counts.csv", header =T, row.names = 1))
genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

# Order the matrices and metadata to be the same order
otu.m <- otu.m[,rownames(metadata.df)]
genus.m <- genus.m[,rownames(metadata.df)]

# Create the rarefied matrix
rarefy_threshold <- 10000
otu_rare.m <- t(rrarefy(t(otu.m[,colSums(otu.m) >= rarefy_threshold]), rarefy_threshold))
genus_rare.m <- t(rrarefy(t(genus.m[,colSums(genus.m) >= rarefy_threshold]), rarefy_threshold))

# Create phyloseq object
otu_rare_phyloseq <- otu_table(otu_rare.m, taxa_are_rows=TRUE)
genus_rare_phyloseq <- otu_table(genus_rare.m, taxa_are_rows=TRUE)

# Estimate alpha diversities
otu_rare_alpha.df <- estimate_richness(otu_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
otu_rare_alpha.df <- otu_rare_alpha.df[rownames(otu_rare_alpha.df) %in% rownames(metadata.df),]

genus_rare_alpha.df <- estimate_richness(genus_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
genus_rare_alpha.df <- genus_rare_alpha.df[rownames(genus_rare_alpha.df) %in% rownames(metadata.df),]

# ---------------------------
# Combine with metadata
otu_rare_alpha.df <- right_join(metadata.df[c("Index", 
                                             discrete_variables, 
                                             grep("_colour", names(metadata.df), value = T))],m2df(otu_rare_alpha.df, "Index"), by = "Index")
genus_rare_alpha.df <- right_join(metadata.df[c("Index",
                                                   discrete_variables, 
                                                   grep("_colour", names(metadata.df), value = T))],m2df(genus_rare_alpha.df, "Index"), by = "Index")

# Remove colour columns
otu_rare_alpha.df <- otu_rare_alpha.df[,!grepl("_colour", names(otu_rare_alpha.df))]
genus_rare_alpha.df <- genus_rare_alpha.df[,!grepl("_colour", names(genus_rare_alpha.df))]
# ---------------------------

# Calculate summary for each variable
otu_alpha_diversity_summary.df <- summarise_diversities_each_variable(otu_rare_alpha.df, 
                                                                        variables = discrete_variables)

genus_alpha_diversity_summary.df <- summarise_diversities_each_variable(genus_rare_alpha.df, 
                                                                        variables = discrete_variables)


# Calculate summary for each variable for each COMMUNITY
# genus_alpha_diversity_summary_per_community.df <- data.frame()
# for (community in as.character(unique(metadata.df$Community))){
#   discrete_variables_subset <- discrete_variables[!discrete_variables %in% c("Community")]
#   genus_rare_alpha_subset.df <- subset(genus_rare_alpha.df, Community == "Remote")
#   temp <- summarise_diversities_each_variable(genus_rare_alpha_subset.df, 
#                                               variables = discrete_variables_subset)
#   temp$Community <- community
#   genus_alpha_diversity_summary_per_community.df <- rbind(genus_alpha_diversity_summary_per_community.df,temp)
# }

# Calculate significances for each variable (and for each COMMUNITY)
otu_alpha_mann_significances.df <- data.frame()
otu_alpha_dunn_significances_multiple.df <- data.frame()

genus_alpha_mann_significances.df <- data.frame()
genus_alpha_dunn_significances_multiple.df <- data.frame()
# genus_alpha_dunn_significances_per_community.df <- data.frame()
# genus_alpha_dunn_significances_multiple_per_community.df <- data.frame()
for (variable in discrete_variables){
  print(variable)
  n_groups <- unique(genus_rare_alpha.df[[variable]])
  if (length(n_groups) < 3){
    print(paste0("Pair: ", variable))
    otu_alpha_mann_significances.df <- rbind(otu_alpha_mann_significances.df, 
                                               calculate_alpha_diversity_significance(otu_rare_alpha.df,
                                                                                      variable = variable))
    
    genus_alpha_mann_significances.df <- rbind(genus_alpha_mann_significances.df, 
                                               calculate_alpha_diversity_significance(genus_rare_alpha.df,
                                                                                      variable = variable))
  } else{
    print(paste0("Multiple: ", variable))
    otu_alpha_dunn_significances_multiple.df <- rbind(otu_alpha_dunn_significances_multiple.df, 
                                                        calculate_alpha_diversity_significance_multiple(otu_rare_alpha.df,
                                                                                                        variable))
    genus_alpha_dunn_significances_multiple.df <- rbind(genus_alpha_dunn_significances_multiple.df, 
                                                        calculate_alpha_diversity_significance_multiple(genus_rare_alpha.df,
                                                                                                        variable))
  }
}

# Remove entries that are not (near) significant
threshold <- 0.05
otu_alpha_mann_significances.df <- otu_alpha_mann_significances.df[apply(otu_alpha_mann_significances.df[,c("Shannon_MannW_padj", "Simpson_MannW_padj", "Chao1_MannW_padj")],1,min) <= threshold,]
genus_alpha_mann_significances.df <- genus_alpha_mann_significances.df[apply(genus_alpha_mann_significances.df[,c("Shannon_MannW_padj", "Simpson_MannW_padj", "Chao1_MannW_padj")],1,min) <= threshold,]

otu_alpha_dunn_significances_multiple.df <- otu_alpha_dunn_significances_multiple.df[apply(otu_alpha_dunn_significances_multiple.df[,c("Shannon_Dunn_padj", "Simpson_Dunn_padj", "Chao1_Dunn_padj")], 1, min) <= threshold,]
genus_alpha_dunn_significances_multiple.df <- genus_alpha_dunn_significances_multiple.df[apply(genus_alpha_dunn_significances_multiple.df[,c("Shannon_Dunn_padj", "Simpson_Dunn_padj", "Chao1_Dunn_padj")], 1, min) <= threshold,]

# Write per-sample diversities to file
write.csv(otu_rare_alpha.df,
          "Result_tables/diversity_analysis/otu/sample_otu_alpha_diversities.csv", quote = F, row.names = F)
write.csv(genus_rare_alpha.df,
          "Result_tables/diversity_analysis/genus/sample_genus_alpha_diversities.csv", quote = F, row.names = F
)
# Write variable summaries to file 
write.csv(otu_alpha_diversity_summary.df,
          "Result_tables/diversity_analysis/otu/otu_alpha_diversities_summary.csv", quote = F, row.names = F
)
write.csv(genus_alpha_diversity_summary.df,
          "Result_tables/diversity_analysis/genus/genus_alpha_diversities_summary.csv", quote = F, row.names = F
)
# Write pair significances
write.csv(otu_alpha_mann_significances.df,
          "Result_tables/diversity_analysis/otu/otu_alpha_diversities_signficances_mannwhitney.csv", quote = F, row.names = F
)
write.csv(genus_alpha_mann_significances.df,
          "Result_tables/diversity_analysis/genus/genus_alpha_diversities_signficances_mannwhitney.csv", quote = F, row.names = F
)
# Write Dunn significances
write.csv(otu_alpha_dunn_significances_multiple.df,
          "Result_tables/diversity_analysis/otu/otu_alpha_diversities_signficances_dunn.csv", quote = F, row.names = F
)
write.csv(genus_alpha_dunn_significances_multiple.df,
          "Result_tables/diversity_analysis/genus/genus_alpha_diversities_signficances_dunn.csv", quote = F, row.names = F
)




