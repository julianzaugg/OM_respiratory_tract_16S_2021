library(DESeq2)
# library(BiocParallel)



############################################################
######################## Functions #########################

# For each rowname (OTU), get the corresponding taxonomy_species
# Assumes "OTU.ID" and "taxonomy_species" columns in the provided map dataframe
assign_taxonomy_to_otu <- function(otutable, taxon_map){
  taxonomies <- c()
  for (otuid in rownames(otutable)){
    taxonomies <- c(taxonomies, as.character(taxon_map[taxon_map$OTU.ID == otuid,]$taxonomy_species))
  }
  return(taxonomies)
}

# Filter and sort DESeq result tables
filter_and_sort_dds_results <- function(x, p_value_threshold = 0.05){
  filtered_table <- x
  filtered_table <- filtered_table[!is.na(filtered_table$padj),]
  filtered_table <- filtered_table[filtered_table$padj <= p_value_threshold,]
  filtered_table <- filtered_table[order(filtered_table$padj),]
  return(filtered_table)
}


filter_matrix_rows <- function(my_matrix, row_max){
  rows_before <- dim(my_matrix)[1]
  rows_after <- dim(my_matrix[apply(my_matrix,1,max) >= row_max,])[1]
  print(paste0("Rows before = ", rows_before))
  print(paste0("Rows after = ", rows_after))
  print(paste0("Lost % = ", round((rows_before-rows_after)/rows_before*100, 2), "%"))
  return(my_matrix[apply(my_matrix,1,max) >= row_max,])
}

############################################################

# Set the working directory
source("code/helper_functions.R")

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Define the discrete variables
discrete_variables <- c("Community", "Nose", "Otitis_Status", "Season", "No_peop_res_discrete")

metadata.df$Tympanic_membrane[metadata.df$Tympanic_membrane == "Unable to visualise/Not examined"] <- NA
metadata.df[metadata.df$Otitis_Status == "Acute Otitis Media","Otitis_Status"] <- NA

# Load count tables
otu.m <- as.matrix(read.table("Result_tables/count_tables/OTU_counts.csv", sep =",", header =T, row.names = 1))
genus.m <- as.matrix(read.table("Result_tables/count_tables/Genus_counts.csv", sep =",", header =T, row.names = 1))

# Filter out features/taxa that do not have at # reads in at least one sample
dim(otu.m)
dim(genus.m)
# head(melt(sort(apply(genus.m, 1, max))), 10)
otu.m <- filter_matrix_rows(otu.m,20)
genus.m <- filter_matrix_rows(genus.m,20)
dim(otu.m)
dim(genus.m)
head(melt(sort(apply(genus.m, 1, max))),10)
# Prevalence
# prevalence = 0.1
# rownames(genus.m)[apply(genus.m, 1, function(x) {length(which(x > 0))}) /length(colnames(genus.m)) >= prevalence]


# Only keep columns (samples) in the metadata
otu.m <- otu.m[,colnames(otu.m) %in% as.character(metadata.df$Index)]
genus.m <- genus.m[,colnames(genus.m) %in% as.character(metadata.df$Index)]

# Order the metadata by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]

# Rownames should match the sample columns in the otu table
rownames(metadata.df) <- metadata.df$Index

# Order the otu_tables the same order as the metadata
otu.m <- otu.m[,rownames(metadata.df)]
genus.m <- genus.m[,rownames(metadata.df)]

# Ensure names of the otu / genus count matrices match the order of the metadata.df!
# Assumes number of samples in metadata.df and count data are the same
all(colnames(otu.m) == metadata.df$Index) # Should be 'True'
all(colnames(otu.m) == rownames(metadata.df)) # Should be 'True'

# Convert variables to factors
metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)

# ---------------------------------------------------------------------------------------------------------
# Perform differential abundance calculations at the OTU level and genus level, 
# comparing between the groups within variables of interest

# DESeq requires a count matrix ('countData'), a corresponding metadata_decontaminated.df ('colData') and a 'design' formula. The formula expresses
# how the counts for each OTU/genus depend on the variables defined in the 'colData'. See help(DESeqDataSetFromMatrix) for more information.
# The first column of the metadata_decontaminated.df ('colData') must match the ordering of the columns of the countData

compare_groups_deseq <- function(mydata.m, mymetadata.df, myvariables, assign_taxonomy = T){
  # Compare groups for all variables
  combined_results_ordered.df <- data.frame()
  
  for (myvar in myvariables){
    print(paste0("Processing ", myvar))
    
    # Get all non-NA entries in the metadata
    mymetadata_filtered.df <- mymetadata.df[!is.na(mymetadata.df[,myvar]),]
    
    # Ensure factored variable
    mymetadata_filtered.df[,myvar] <- factor(mymetadata_filtered.df[,myvar])
    
    # Extract corresponding entries from data
    mydata_filtered.m <- mydata.m[,rownames(mymetadata_filtered.df)]
    
    # If the number of samples is 1 or there is only one unique variable
    if (dim(mymetadata_filtered.df)[2] == 1 | length(unique(mymetadata_filtered.df[,myvar])) == 1){
      print("Only one sample or only one unique group")
      next
    }
    if (dim(mymetadata_filtered.df)[1] == 0 | dim(mydata_filtered.m)[2] == 0){
      print("No samples after filtering")
      next
    }
    
    # If the column and rownames do not match, entries are missing
    if (!all(rownames(mymetadata_filtered.df) == colnames(mydata_filtered.m))){
      print("Colnames and metadata names don't match!!!")
      next
    }
    
    # Run DESeq
    dds <- DESeqDataSetFromMatrix(countData = mydata_filtered.m, colData = mymetadata_filtered.df, design = as.formula(paste0("~", myvar)))
    
    geoMeans <- apply(counts(dds), 1, gm_mean)
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = F))
    if(inherits(dds, "try-error")) {
      next
    }
    group_combinations <- combn(sort(unique(mymetadata_filtered.df[,myvar])),2)
    
    for (i in 1:ncol(group_combinations)){
      group_1 <- as.character(group_combinations[1,i])
      group_2 <- as.character(group_combinations[2,i])
      
      # group_1_meta <- subset(full, get(myvar) == group_1)
      # group_2_meta <- subset(full, get(myvar) == group_2)
      # n_group_1 <- dim(group_1_meta)[1]
      # n_group_2 <- dim(group_2_meta)[1]
      
      n_group_1 <- dim(subset(mymetadata_filtered.df, get(myvar) == group_1))[1]
      n_group_2 <- dim(subset(mymetadata_filtered.df, get(myvar) == group_2))[1]
      
      # Extract results for contrasted groups
      print(paste0(myvar, ": ", group_1, " vs ", group_2))
      resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.05, independentFiltering = F, cooksCutoff = F,parallel = T)
      # print(resMFSource)
      resMFSource$Group_1 <- group_1
      resMFSource$Group_2 <- group_2
      resMFSource$Variable <- myvar
      resMFSource$N_Group_1 <- n_group_1
      resMFSource$N_Group_2 <- n_group_2
      
      # Assign the taxonomy to the results. Assumes feature.
      if (assign_taxonomy == T){
        resMFSource$Taxonomy <- assign_taxonomy_to_otu(resMFSource, otu_taxonomy_map.df)   
        # Convert to dataframe
        resMFSource <- m2df(resMFSource, "OTU")
      } else{
        # Convert to dataframe
        resMFSource <- m2df(resMFSource, "Taxonomy")
      }
      # print(resMFSource)
      resMFSource <- filter_and_sort_dds_results(resMFSource, 0.05)
      combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
    }
  }
  combined_results_ordered.df
}

compare_groups_deseq_within_group <- function(mydata.m, mymetadata.df, myvariables, within_group_variable, assign_taxonomy = F){
  combined_results.df <- data.frame()
  myvariables <- myvariables[which(!myvariables == within_group_variable)]
  for (myvar_value in unique(mymetadata.df[,within_group_variable])){
    print(paste0("Processing ", myvar_value))
    temp <- compare_groups_deseq(mydata.m = mydata.m, 
                                 mymetadata.df = subset(mymetadata.df, get(within_group_variable) == myvar_value), 
                                 myvariables = myvariables, 
                                 assign_taxonomy = assign_taxonomy)
    if (dim(temp)[1] == 0){
      next
    }
    temp[,within_group_variable] <- myvar_value
    combined_results.df <- rbind(combined_results.df, temp)
  }
  combined_results.df
}

otu_group_comparison.df <- compare_groups_deseq(mydata.m = otu.m, mymetadata.df = metadata.df, myvariables = discrete_variables, assign_taxonomy = T)
genus_group_comparison.df <- compare_groups_deseq(mydata.m = genus.m, mymetadata.df = metadata.df, myvariables = discrete_variables, assign_taxonomy = F)

write.csv(x = otu_group_comparison.df,file ="Result_tables/DESeq_results/OTU_deseq.csv",quote = F, row.names =F)
write.csv(x = genus_group_comparison.df,file ="Result_tables/DESeq_results/Genus_deseq.csv",quote = F, row.names =F)

# temp <- read.csv("Result_tables/DESeq_results/OTU_deseq.csv", header =T )
# temp <- left_join(temp, otu_taxonomy_map.df[,c("OTU.ID", "ASV_ID")], by = c("OTU" = "OTU.ID"))
# temp <- temp %>% relocate(c("OTU", "ASV_ID")) %>% mutate(ASV_ID = gsub("_", "", ASV_ID))
# write.csv(x = temp,file ="Result_tables/DESeq_results/OTU_deseq.csv",quote = F, row.names =F)

# ----------------------------------------------------------------------------------------------
reduced_variables <- discrete_variables[which(!discrete_variables == "Community")]
otu_group_comparison_within_community.df <- compare_groups_deseq_within_group(mydata.m = otu.m, 
                                                                              mymetadata.df = metadata.df, 
                                                                              myvariables = reduced_variables, 
                                                                              within_group_variable = "Community", 
                                                                              assign_taxonomy = T)

genus_group_comparison_within_community.df <- compare_groups_deseq_within_group(mydata.m = genus.m, 
                                                                                mymetadata.df = metadata.df, 
                                                                                myvariables = reduced_variables, 
                                                                                within_group_variable = "Community", 
                                                                                assign_taxonomy = F)

write.csv(x =otu_group_comparison_within_community.df,file ="Result_tables/DESeq_results/OTU_deseq_within_community.csv",quote = F, row.names =F)
write.csv(x =genus_group_comparison_within_community.df,file ="Result_tables/DESeq_results/Genus_deseq_within_community.csv",quote = F, row.names =F)

temp <- read.csv("Result_tables/DESeq_results/OTU_deseq_within_community.csv", header =T )
temp <- left_join(temp, otu_taxonomy_map.df[,c("OTU.ID", "ASV_ID")], by = c("OTU" = "OTU.ID"))
temp <- temp %>% relocate(c("OTU", "ASV_ID")) %>% mutate(ASV_ID = gsub("_", "", ASV_ID))
write.csv(x = temp,file ="Result_tables/DESeq_results/OTU_deseq_within_community.csv",quote = F, row.names =F)
# ----------------------------------------------------------------------------------------------
reduced_variables <- discrete_variables[which(!discrete_variables == "Otitis_Status")]

otu_group_comparison_within_otitis_status.df <- compare_groups_deseq_within_group(mydata.m = otu.m, 
                                                                              mymetadata.df = metadata.df, 
                                                                              myvariables = reduced_variables, 
                                                                              within_group_variable = "Otitis_Status", 
                                                                              assign_taxonomy = T)

genus_group_comparison_within_otitis_status.df <- compare_groups_deseq_within_group(mydata.m = genus.m, 
                                                                                mymetadata.df = metadata.df, 
                                                                                myvariables = reduced_variables, 
                                                                                within_group_variable = "Otitis_Status", 
                                                                                assign_taxonomy = F)

write.csv(x =otu_group_comparison_within_otitis_status.df,file ="Result_tables/DESeq_results/OTU_deseq_within_Otitis_Status.csv",quote = F, row.names =F)
write.csv(x =genus_group_comparison_within_otitis_status.df,file ="Result_tables/DESeq_results/Genus_deseq_within_Otitis_Status.csv",quote = F, row.names =F)

# ----------------------------------------------------------------------------------------------
reduced_variables <- discrete_variables[which(!discrete_variables == "Nose")]
otu_group_comparison_within_nose.df <- compare_groups_deseq_within_group(mydata.m = otu.m, 
                                                                              mymetadata.df = metadata.df, 
                                                                              myvariables = reduced_variables, 
                                                                              within_group_variable = "Nose", 
                                                                              assign_taxonomy = T)

genus_group_comparison_within_nose.df <- compare_groups_deseq_within_group(mydata.m = genus.m, 
                                                                                mymetadata.df = metadata.df, 
                                                                                myvariables = reduced_variables, 
                                                                                within_group_variable = "Nose", 
                                                                                assign_taxonomy = F)

write.csv(x =otu_group_comparison_within_nose.df,file ="Result_tables/DESeq_results/OTU_deseq_within_Nose.csv",quote = F, row.names =F)
write.csv(x =genus_group_comparison_within_nose.df,file ="Result_tables/DESeq_results/Genus_deseq_within_Nose.csv",quote = F, row.names =F)

# ----------------------------------------------------------------------------------------------
reduced_variables <- discrete_variables[which(!discrete_variables == "Community__Otitis_Status")]
otu_group_comparison_within_community_otitis_status.df <- compare_groups_deseq_within_group(mydata.m = otu.m, 
                                                                              mymetadata.df = metadata.df, 
                                                                              myvariables = reduced_variables, 
                                                                              within_group_variable = "Community__Otitis_Status", 
                                                                              assign_taxonomy = T)

genus_group_comparison_within_community_otitis_status.df <- compare_groups_deseq_within_group(mydata.m = genus.m, 
                                                                                mymetadata.df = metadata.df, 
                                                                                myvariables = reduced_variables, 
                                                                                within_group_variable = "Community__Otitis_Status", 
                                                                                assign_taxonomy = F)

write.csv(x =otu_group_comparison_within_community_otitis_status.df,file ="Result_tables/DESeq_results/OTU_deseq_within_Community__Otitis_Status.csv",quote = F, row.names =F)
write.csv(x =genus_group_comparison_within_community_otitis_status.df,file ="Result_tables/DESeq_results/Genus_deseq_within_Community__Otitis_Status.csv",quote = F, row.names =F)


# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------