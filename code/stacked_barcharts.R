library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)

source("code/helper_functions.R")

# Load metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv")

# Remove AOM, just make values NA
# metadata.df[metadata.df$Otitis_Status == "Acute Otitis Media","Otitis_Status"] <- NA

# Set order of Otitis Status
metadata.df$Otitis_Status <- factor(metadata.df$Otitis_Status, levels = c("Never OM", "HxOM", "Acute Otitis Media", "Effusion","Perforation"))
# Set order of Season
metadata.df$Season <- factor(metadata.df$Season, levels = c("Spring", "Winter", "Autumn"))
# Set order of Nose
metadata.df$Nose <- factor(metadata.df$Nose, levels = c("Normal", "Serous", "Purulent"))
# Set order of No_peop_res_discrete
metadata.df$No_peop_res_discrete <- factor(metadata.df$No_peop_res_discrete, levels = c("2 to 3", "4 to 6", "7 to 12", "Unknown"))

# Ensure rownames are the Index
rownames(metadata.df) <- metadata.df$Index

# Define variables of interest
variables_of_interest <- c("Community", "Nose", "Otitis_Status", "Season", "No_peop_res_discrete")

# Load the OTU - taxonomy mapping file
otu_taxonomy_map.df <- read.csv("Result_tables/other/otu_taxonomy_map.csv", header = T)

# Load abundance data
otu_rel.df <- read.csv("Result_tables/relative_abundance_tables/OTU_relative_abundances.csv", header = T)
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)

# First create a global palette of the top taxa (more than we need though not everything)
palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
palette_20b <- c("#d05c2a","#a9b635","#b3a95e","#3fbcc1","#397f50","#e19073","#68752b","#d790ce","#7663cd","#d99c34","#97518d","#cf4248","#63c34e","#5fbf85","#d94983","#59993a","#658acc","#ad5262","#c653bd","#95642d")
palette_20c <- c("#a84b54","#a2b432","#c057c5","#6f64cf","#5aae36","#527c2e","#826729","#6690cf","#7fb875","#d68c62","#a8a352","#e484a0","#d2408e","#d49731","#ca5729","#3a8862","#49c1b5","#ab6daa","#dc3f53","#4bc06d")
palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
palette_100_distinct <- c("#b2e54d","#acb5ff","#7ce5e2","#ab3d00","#a4e75a","#699164","#00b9e2","#5a3e6b","#dbcf00","#26511d","#00b979","#f7d152","#008885","#0072b9","#493d85","#6a67fe","#01cdbb","#f8a7ff","#ff54ac","#a8b800","#e675ff","#3080ff","#edd37c","#15ca42","#d1db7b","#014cb1","#01d179","#800489","#ffa27c","#b70094","#ffcc6a","#c1db03","#00ac42","#ffb971","#e00097","#5aef87","#fac0ff","#ff74b9","#f1053b","#1f5a00","#898ac0","#7fbf00","#016cc8","#514804","#00578b","#abc994","#4ba300","#ffad4d","#9720ba","#cd27bf","#ff68d7","#01b0b1","#bd6700","#87a1ff","#a20019","#842708","#ea0090","#006e04","#facc9c","#8ae6c3","#b28400","#4d8cff","#d86400","#84223c","#918b55","#234198","#ff8242","#6d3260","#e31db5","#7c7f00","#3a8500","#70ec9f","#4a5700","#a065fd","#2f41d0","#f0a300","#bc0042","#c184ff","#93e973","#6d4400","#58699b","#ff8d1c","#00997f","#ff6655","#88066c","#b6e190","#ff2899","#ce2a00","#ff2b91","#a20062","#7e2d23","#ff7c88","#783403","#009f49","#a29600","#017a67","#ffaeb6","#4e329e","#00722f","#5e23a7")

# Determine the top genus by mean abundance
top_genus <- abundance_summary(genus_rel.df,top_n = 20)$taxonomy_genus

# Create palette for the genus
genus_palette <- setNames(c(palette_20c[1:length(top_genus)], "grey"), c(top_genus, "Other"))
genus_palette


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# All samples (grouped by Otitis_Status)

# otitis_status_samples.v <- metadata.df[!is.na(metadata.df$Otitis_Status),]$Index

# Collapse abundance dataframe
genus_rel_collapsed.df <- collapse_abundance_dataframe(genus_rel.df, top_genus)
genus_rel_collapsed.df <- genus_rel_collapsed.df[!genus_rel_collapsed.df$Sample %in% metadata.df[is.na(metadata.df$Otitis_Status),]$Index,]

# Combine with metadata
genus_data.df <- left_join(genus_rel_collapsed.df, metadata.df[,c("Index", 
                                                                  "Otitis_Status","Otitis_Status_colour",
                                                                  "Community", "Community_colour",
                                                                  "Nose", "Nose_colour",
                                                                  "Season", "Season_colour",
                                                                  "No_peop_res_discrete", "No_peop_res_discrete_colour"),
                                                               ], by = c("Sample"= "Index"))

# Get order of top taxa by abundance
temp <- genus_data.df[,c("taxonomy_genus","Sample", "Relative_abundance")]
temp <- df2matrix(dcast(temp,taxonomy_genus~Sample))
ordered_taxa <- names(sort(apply(temp,1,mean),decreasing = T))
ordered_taxa <- c(ordered_taxa[!ordered_taxa == "Other"],"Other")

# Order by abundance
genus_data.df <- arrange(genus_data.df, dplyr::desc(Relative_abundance),taxonomy_genus)
sample_order <- genus_data.df %>% filter(grepl("g__Dolosigranulum", taxonomy_genus)) %>% arrange(Relative_abundance) %>% pull(Sample) %>% as.character()
genus_data.df <- genus_data.df %>% arrange(match(Sample,sample_order))

# Order by Otitis status and Community
genus_data.df <- 
  genus_data.df %>% 
  arrange(
    Otitis_Status,
    Community
    )

# Factorise sample
genus_data.df$Sample <- factor(genus_data.df$Sample, 
                                    levels = unique(as.character(genus_data.df$Sample)))
# Factorise taxonomy
genus_data.df$taxonomy_genus <- factor(genus_data.df$taxonomy_genus, levels = ordered_taxa)

# Create label for the taxonomy
genus_data.df$Label <- unlist(lapply(genus_data.df$taxonomy_genus,
                                          function(x) ifelse(grepl(";", x),
                                                             first_resolved_taxonomy(as.character(x)),
                                                             as.character(x)))
)

labels <- unlist(lapply(as.character(unique(genus_data.df$taxonomy_genus)),function(x) ifelse(grepl(";", x),
                                                                                                   first_resolved_taxonomy(as.character(x)),
                                                                                                   as.character(x))))

label_list <- setNames(labels,as.character(unique(genus_data.df$taxonomy_genus)))

# Create palette lists
otitis_status_palette <- unique(genus_data.df[,c("Otitis_Status", "Otitis_Status_colour")])
otitis_status_palette <- setNames(otitis_status_palette$Otitis_Status_colour,otitis_status_palette$Otitis_Status)

community_palette <- unique(genus_data.df[,c("Community", "Community_colour")])
community_palette <- setNames(community_palette$Community_colour,community_palette$Community)

season_palette <- unique(genus_data.df[,c("Season", "Season_colour")])
season_palette <- setNames(season_palette$Season_colour,season_palette$Season)

nose_palette <- unique(genus_data.df[,c("Nose", "Nose_colour")])
nose_palette <- setNames(nose_palette$Nose_colour,nose_palette$Nose)

no_peop_res_discrete_palette <- unique(genus_data.df[,c("No_peop_res_discrete", "No_peop_res_discrete_colour")])
no_peop_res_discrete_palette <- setNames(no_peop_res_discrete_palette$No_peop_res_discrete_colour,no_peop_res_discrete_palette$No_peop_res_discrete)

# temp <- genus_data.df[grepl("g__Rothia",genus_data.df$taxonomy_genus),]
# temp <- temp[temp$Relative_abundance > 0,]

myplot <- ggplot(genus_data.df, aes(x = Sample, y = Relative_abundance * 100,fill = taxonomy_genus)) +
  # coord_flip() +
  # geom_vline(data = temp,aes(xintercept = Sample)) +
  geom_bar(stat = "identity", colour = "black", lwd = 0.1, width = 0.75) +
  xlab("Sample") + 
  ylab("Relative abundance %") + 
  
  scale_fill_manual(values = genus_palette, 
                    labels = label_list,
                    name = "Genus\n(or lowest resolved taxonomy level)") +
  
  new_scale_fill() +
  geom_rect(aes(xmin = as.numeric(Sample)-.5,xmax = as.numeric(Sample)+.5, fill = No_peop_res_discrete),
            ymin = 101,ymax = 104, colour = "black", lwd = .1)+
  scale_fill_manual(values = no_peop_res_discrete_palette,
                    name = "Household size") +
  
  new_scale_fill() +
  geom_rect(aes(xmin = as.numeric(Sample)-.5,xmax = as.numeric(Sample)+.5, fill = Season),
            ymin = 105,ymax = 108, colour = "black", lwd = .1)+
  scale_fill_manual(values = season_palette,
                    name = "Season") +
  
  new_scale_fill() +
  geom_rect(aes(xmin = as.numeric(Sample)-.5,xmax = as.numeric(Sample)+.5, fill = Nose),
            ymin = 109,ymax = 112, colour = "black", lwd = .1)+
  scale_fill_manual(values = nose_palette,
                    name = "Nose") +
  
  new_scale_fill() +
  geom_rect(aes(xmin = as.numeric(Sample)-.5,xmax = as.numeric(Sample)+.5, fill = Community),
            ymin = 113,ymax = 116, colour = "black", lwd = .1)+
  scale_fill_manual(values = community_palette,
                    name = "Community") +
  new_scale_fill() +
  geom_rect(aes(xmin = as.numeric(Sample)-.5,xmax = as.numeric(Sample)+.5, fill = Otitis_Status),
            ymin = 117,ymax = 120, colour = "black", lwd = .1)+
  scale_fill_manual(values = otitis_status_palette,
                    name = "Otitis Status") +
  
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10), limits = c(0,121)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white", size = 1),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 10,face = "bold"),
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1),
    legend.title = element_text(hjust = .5, face = "bold"),
    # legend.key.size = unit(.1, "cm")
    legend.key.width = unit(.5, "cm"),
    legend.key.height = unit(.5, "cm"),
    legend.position = "bottom",
    legend.direction = "vertical"
  )
myplot
ggsave(
  filename = "Result_figures/abundance_analysis_plots/genus_all_samples_stacked_barchart.pdf",
  plot = myplot,
  height = 25,
  width = 55,
  units = "cm")

ggsave(
  filename = "Result_figures/abundance_analysis_plots/genus_all_samples_stacked_barchart.svg",
  plot = myplot,
  height = 25,
  width = 55,
  units = "cm")
