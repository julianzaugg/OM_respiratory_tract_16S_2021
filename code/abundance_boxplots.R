library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)


# Calculate the significance values for taxa between multiple groups
calculate_taxa_significances_multiple <- function(mydata, variable_column, value_column, taxonomy_column){
  results.df <- data.frame("Taxonomy" = character(),
                           "Variable" = character(),
                           "Group_1" = character(),
                           "Group_2" = character(),
                           "Dunn_pvalue" = character(),
                           "Dunn_padj" = character(),
                           "KrusW_pvalue" = character()
  )
  for (taxa in unique(mydata[,taxonomy_column])){ # For each taxa in the taxonomy column
    taxa_data <- subset(mydata, get(taxonomy_column) == taxa)
    n_groups = length(as.character(unique(taxa_data[,variable_column])))
    if (any(is.na(taxa_data[,variable_column]))){
      return()
    }
    if (n_groups > 2){
      kw <- kruskal.test(get(value_column)~get(variable_column), data = taxa_data)
      dunn <- dunnTest(x = get(value_column)~get(variable_column), data = taxa_data, method = "bh", alpha = 0.05)
      dunn <- separate(dunn$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
      names(dunn) <- c("Group_1","Group_2","Dunn_pvalue","Dunn_padj")
      dunn$Taxonomy <- taxa
      dunn$KrusW_pvalue <- kw$p.value
      dunn$Variable <- variable_column
      results.df <- rbind(results.df, dunn)
    }
  }
  results.df$P_value_label <- as.character(lapply(results.df[,"Dunn_padj"], function(x) ifelse(x <= 0.001, "***", 
                                                                                                 ifelse(x <= 0.01, "**", 
                                                                                                        ifelse(x <= 0.05, "*", "ns")))))
  
  results.df[,c("Taxonomy", "Variable", "Group_1","Group_2", "Dunn_pvalue", "Dunn_padj", "KrusW_pvalue")]
}


generate_p_labels <- function(sig_table){
  for (sig_column in c("Dunn_padj")){
    metric = strsplit(sig_column, "_")[[1]][1]
    sig_table[,paste0(metric, "_p_label")] <-
      as.character(lapply(sig_table[,sig_column], 
                          function(x) ifelse(x <= 0.001, "***", 
                                             ifelse(x <= 0.01, "**",
                                                    ifelse(x <= 0.05, "*", "ns")))))
  }
  sig_table
}

source("Code/helper_functions.R")

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# ---------------------------------------------------------------
# Remove AOM, just make values NA
metadata.df[metadata.df$Otitis_Status == "Acute Otitis Media","Otitis_Status"] <- NA
metadata.df <- metadata.df[!is.na(metadata.df$Otitis_Status),]

# Set order of Otitis Status
metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status , levels = c("Never OM", "HxOM", "Effusion","Perforation"))
# Set order of Season
metadata.df$Season <- factor(metadata.df$Season, levels = c("Spring", "Winter", "Autumn"))
# Set order of Nose
metadata.df$Nose <- factor(metadata.df$Nose, levels = c("Normal", "Serous", "Purulent"))
# Set order of No_peop_res_discrete
metadata.df$No_peop_res_discrete <- factor(metadata.df$No_peop_res_discrete, levels = c("2 to 3", "4 to 6", "7 to 12", "Unknown"))
# ---------------------------------------------------------------
variables_of_interest <- c("Community", "Nose", "Otitis_Status", "Season", "No_peop_res_discrete")

# Ensure rownames are the Index
rownames(metadata.df) <- metadata.df$Index

# Load abundance data
otu_rel.df <- read.csv("Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", header = T)
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)

# Remove entries not in metadata
genus_rel.df <- genus_rel.df[,names(genus_rel.df) %in% c("taxonomy_genus", metadata.df$Index)]

# Get top genus for samples
top_genus_all_samples <- abundance_summary(genus_rel.df,20)$taxonomy_genus

# Create palette
palette_20c <- c("#a84b54","#a2b432","#c057c5","#6f64cf","#5aae36","#527c2e","#826729","#6690cf","#7fb875","#d68c62","#a8a352","#e484a0","#d2408e","#d49731","#ca5729","#3a8862","#49c1b5","#ab6daa","#dc3f53","#4bc06d")
# genus_palette <- setNames(c(palette_20c[1:length(top_genus_all_samples)], "grey"), c(top_genus_all_samples, "Other"))

# Get counts for top genus and combine with metadata
counts.df <- read.csv("Result_tables/count_tables/Genus_counts.csv")
counts.df <- counts.df[counts.df$taxonomy_genus %in% top_genus_all_samples,]
counts.df <- collapse_abundance_dataframe(counts.df)
counts.df$Label <- as.character(lapply(counts.df$taxonomy_genus,first_resolved_taxonomy))
counts.df <- left_join(counts.df, metadata.df, by = c("Sample" = "Index"))
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# All samples

# Filter abundance dataframe to top genus
genus_rel_filtered.df <- genus_rel.df[genus_rel.df$taxonomy_genus %in% top_genus_all_samples,]

# Collapse abundance dataframe to taxa and relative abundance
genus_rel_filtered.df <- collapse_abundance_dataframe(genus_rel_filtered.df)

# Combine with metadata
genus_rel_filtered.df <- left_join(genus_rel_filtered.df, metadata.df, by = c("Sample"= "Index"))

# Create label for the taxonomy
genus_rel_filtered.df$Label <- as.character(lapply(genus_rel_filtered.df$taxonomy_genus,first_resolved_taxonomy))

# ------------------------------------------------------------------------------------------------
# Calculate the significances for taxa for Otitis status variable
otitis_status_genus_significances.df <- calculate_taxa_significances_multiple(mydata = genus_rel_filtered.df, 
                                                                              variable_column = "Otitis_Status",
                                                                              value_column = "Relative_abundance",
                                                                              taxonomy_column = "taxonomy_genus")
otitis_status_genus_significances.df <- subset(otitis_status_genus_significances.df, Dunn_padj <= 0.05)
otitis_status_genus_significances.df <- otitis_status_genus_significances.df[order(otitis_status_genus_significances.df$Taxonomy),]
write.csv(file = "Result_tables/abundance_analysis_tables/otitis_status_genus_significances.csv",
          x = otitis_status_genus_significances.df, row.names = F, quote = F)

# Also for each community
otitis_status_genus_significances_remote.df <- calculate_taxa_significances_multiple(mydata = subset(genus_rel_filtered.df, Community == "Remote"), 
                                                                              variable_column = "Otitis_Status",
                                                                              value_column = "Relative_abundance",
                                                                              taxonomy_column = "taxonomy_genus")
otitis_status_genus_significances_remote.df <- subset(otitis_status_genus_significances_remote.df, Dunn_padj <= 0.05)
otitis_status_genus_significances_remote.df <- otitis_status_genus_significances_remote.df[order(otitis_status_genus_significances_remote.df$Taxonomy),]
otitis_status_genus_significances_remote.df$Community <- "Remote"

otitis_status_genus_significances_rural.df <- calculate_taxa_significances_multiple(mydata = subset(genus_rel_filtered.df, Community == "Rural"), 
                                                                                    variable_column = "Otitis_Status",
                                                                                    value_column = "Relative_abundance",
                                                                                    taxonomy_column = "taxonomy_genus")
otitis_status_genus_significances_rural.df <- subset(otitis_status_genus_significances_rural.df, Dunn_padj <= 0.05)
otitis_status_genus_significances_rural.df <- otitis_status_genus_significances_rural.df[order(otitis_status_genus_significances_rural.df$Taxonomy),]
otitis_status_genus_significances_rural.df$Community <- "Rural"

otitis_status_genus_significances_community.df <- rbind(otitis_status_genus_significances_remote.df, otitis_status_genus_significances_rural.df)

write.csv(file = "Result_tables/abundance_analysis_tables/otitis_status_community_genus_significances.csv",
          x = otitis_status_genus_significances_community.df, row.names = F, quote = F)
# ------------------------------------------------------------------------------------------------

# Get abundances for just Ornithobacterium
orni_abundances.df <- genus_rel_filtered.df[grepl("Ornithobacterium", genus_rel_filtered.df$taxonomy_genus),c("Sample",
                                                                                                              "Label",
                                                                                                           "taxonomy_genus",
                                                                                                           "Relative_abundance",
                                                                                                           "Otitis_Status",
                                                                                                           "Otitis_Status_colour",
                                                                                                           "Community",
                                                                                                           "Community__Otitis_Status")]

# Create colour palette for Otitis_Status
colours.l <- unique(genus_rel_filtered.df[,c("Otitis_Status", "Otitis_Status_colour")])
colours.l <- setNames(as.character(colours.l[,"Otitis_Status_colour"]), 
                      colours.l[,"Otitis_Status"])

generate_p_labels(otitis_status_genus_significances_remote.df)
generate_p_labels(otitis_status_genus_significances_rural.df)

myplot <- 
  ggplot(orni_abundances.df, aes(x = Otitis_Status, y = Relative_abundance*100,fill = Otitis_Status, shape = Otitis_Status)) +
  geom_boxplot(position = position_dodge(width =.75), 
               outlier.shape = NA, width=.5,lwd =.2) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = 0.75,
                                              dodge.width = .75)) +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = .8,
               position = position_dodge(width = .75),show.legend = F) +
  xlab("Otitis status") +
  ylab("Relative abundance %") +
  ggtitle(expression(paste(italic(Ornithobacterium)," abundance"))) +
  scale_shape_manual(values = c(25,24,23,22,21),name = "Otitis status") +
  scale_fill_manual(values = colours.l, name = "Otitis status") +
  scale_y_continuous(limits = c(0,30), breaks = seq(0,100, by = 5)) +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black"),#face = "italic"
    axis.text.y = element_text(size = 6, colour = "black"),#face = "italic"
    axis.title = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_blank()
  )
ggsave(filename = "Result_figures/abundance_analysis_plots/ornithobacterium_abundance.pdf",
       units = "cm",
       width = 8,
       height = 6,
       device = "pdf")

ggsave(filename = "Result_figures/abundance_analysis_plots/ornithobacterium_abundance.svg",
       units = "cm",
       width = 8,
       height = 6,
       device = "svg")

myplot <- 
  ggplot(orni_abundances.df, aes(x = Otitis_Status, y = Relative_abundance*100,fill = Otitis_Status, shape = Otitis_Status)) +
  geom_boxplot(position = position_dodge(width =.75), 
               outlier.shape = NA, width=.5,lwd =.2) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = 0.75,
                                              dodge.width = .75)) +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = .8,
               position = position_dodge(width = .75),show.legend = F) +
  facet_wrap(~Community) +
  xlab("Otitis status") +
  ylab("Relative abundance %") +
  ggtitle(expression(paste(italic(Ornithobacterium)," abundance"))) +
  scale_shape_manual(values = c(25,24,23,22,21),name = "Otitis status") +
  scale_fill_manual(values = colours.l, name = "Otitis status") +
  scale_y_continuous(limits = c(0,30), breaks = seq(0,100, by = 5)) +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black"),#face = "italic"
    axis.text.y = element_text(size = 6, colour = "black"),#face = "italic"
    axis.title = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_blank()
  )
ggsave(filename = "Result_figures/abundance_analysis_plots/ornithobacterium_abundance_community.pdf",
       units = "cm",
       width = 12,
       height = 6,
       device = "pdf")

ggsave(filename = "Result_figures/abundance_analysis_plots/ornithobacterium_abundance_community.svg",
       units = "cm",
       width = 12,
       height = 6,
       device = "svg")

otitis_top_genus_boxplot <- 
  ggplot(genus_rel_filtered.df, aes(x = Label, 
                                  y = Relative_abundance*100,
                                  fill = Otitis_Status,
                                  shape = Otitis_Status)) +
  geom_boxplot(position = position_dodge(width =.75), 
               outlier.shape = NA, width=.5,lwd =.2) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = 0,
                                              dodge.width = .75)) +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = .8,
               position = position_dodge(width = .75),show.legend = F) +
  xlab("Taxonomy") +
  ylab("Relative abundance %") +
  scale_shape_manual(values = c(25,24,23,22,21),name = "Otitis status") +
  scale_fill_manual(values = colours.l, name = "Otitis status") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 10)) + 
  # facet_wrap(~Community) +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black", angle = 45,vjust = 1,hjust = 1),#face = "italic"
    axis.title = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 5),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_blank()
  )
ggsave(filename = "Result_figures/abundance_analysis_plots/otitis_top_genus_boxplot.pdf",
       plot = otitis_top_genus_boxplot,
       width = 25,
       height = 7,
       units = "cm")


otitis_top_genus_mean_sd_plot <-  ggplot(genus_rel_filtered.df, aes(x = Label, 
                                  y = Relative_abundance*100,
                                  colour = Otitis_Status,
                                  fill = Otitis_Status,
                                  shape = Otitis_Status)) +
  
  # scale_y_log10() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),colour="black",
               geom="errorbar", width=0.5,lwd=.8, position = position_dodge(width = 0.75)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0.4,lwd=.3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),colour = "grey10",
               geom="point", position = position_dodge(width = 0.75)) +
  xlab("Taxonomy") +
  ylab("Relative abundance %") +
  scale_shape_manual(values = c(25,24,23,22,21),name = "Otitis status") +
  scale_fill_manual(values = colours.l, name = "Otitis status") +
  scale_colour_manual(values = colours.l, name = "Otitis status") +
  scale_y_continuous(breaks = seq(0,70, by = 10),expand = c(0,0)) +
  coord_cartesian(ylim = c(0,70)) +
  theme(
    # axis.text.x = element_text(angle = 90, size = 6,hjust = 0, vjust = 0.5),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major.y = element_line(colour = "grey50",size = .3,linetype = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black"),#face = "italic"
    axis.title = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 5),
    legend.key.size = unit(.3,"cm"),
    legend.key = element_blank()
  )

ggsave(filename = "Result_figures/abundance_analysis_plots/otitis_top_genus_mean_sd_plot.pdf",
       plot = otitis_top_genus_mean_sd_plot,
       width = 25,
       height = 7,
       units = "cm")

